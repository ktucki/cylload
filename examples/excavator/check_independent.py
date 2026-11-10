#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Independent check of the CylLoad results for the excavator example.

A separate implementation in Python/NumPy that shares no code with CylLoad and uses a
different formulation: absolute (maximal) coordinates (x, y, theta) for every body, the
joints and the actuator lengths as algebraic constraints, and Lagrange multipliers for the
joint reactions and the actuator forces - the formulation of general-purpose multibody
codes. CylLoad uses relative coordinates and virtual work. The modelling assumptions are
the same (rigid links, half of each cylinder mass lumped at each of its pins).

Checked over the stroke of each cylinder (CylLoad sweep results):
  * actuator force          F   from  Q_g + Phi_q^T lambda + sum F_k dL_k/dq = 0
  * reduced mass            m_r = qdot^T M qdot for a unit velocity of the actuator
                                  (the other actuators locked)
  * hydraulic stiffness and natural frequency f = sqrt(k_h / m_r) / (2 pi)
Usage (from the repository root, after examples/excavator/run.sh):
  python3 examples/excavator/check_independent.py
"""
import csv, math, os, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
G = 9.81


def read_model(path):
    bodies, joints, acts, hyd, masses = ['ground'], [], [], {}, None
    for line in open(path):
        t = line.split('#')[0].split()
        if not t: continue
        if t[0] == 'body': bodies += t[1:]
        elif t[0] == 'revolute': joints.append((t[1], t[2], t[3], float(t[4]), float(t[5])))
        elif t[0] == 'actuator':
            kv = dict(x.split('=') for x in t[8:])
            acts.append(dict(name=t[1], a=t[2], pa=(float(t[3]), float(t[4])), b=t[5], pb=(float(t[6]), float(t[7])),
                             bore=float(kv['bore']), rod=float(kv['rod']), lmin=float(kv['min']), lmax=float(kv['max']),
                             mass=float(kv.get('mass', 0))))
        elif t[0] == 'hydraulics': hyd = dict(x.split('=') for x in t[1:])
        elif t[0] == 'masses': masses = os.path.join(os.path.dirname(path), t[1])
    items = []
    for r in csv.DictReader(open(masses)):
        items.append((r['body'], float(r['mass_kg']), float(r['x_mm']) / 1000, float(r['y_mm']) / 1000, float(r['J_kgm2'])))
    for a in acts:  # half of the cylinder mass at each pin
        items.append((a['a'], a['mass'] / 2, a['pa'][0] / 1000, a['pa'][1] / 1000, 0.0))
        items.append((a['b'], a['mass'] / 2, a['pb'][0] / 1000, a['pb'][1] / 1000, 0.0))
    return bodies, joints, acts, hyd, items


class Mech:
    def __init__(self, path):
        self.bodies, self.joints, self.acts, self.hyd, items = read_model(path)
        self.n = len(self.bodies) - 1
        self.idx = {b: i - 1 for i, b in enumerate(self.bodies)}  # ground -> -1
        # body inertia about the body frame origin: mass, first moments, polar moment
        self.m = np.zeros(self.n); self.S = np.zeros((self.n, 2)); self.Jo = np.zeros(self.n)
        for b, m, x, y, J in items:
            i = self.idx[b]
            if i < 0: continue
            self.m[i] += m; self.S[i] += m * np.array([x, y]); self.Jo[i] += J + m * (x * x + y * y)

    @staticmethod
    def A(th): return np.array([[math.cos(th), -math.sin(th)], [math.sin(th), math.cos(th)]])
    @staticmethod
    def dA(th): return np.array([[-math.sin(th), -math.cos(th)], [math.cos(th), -math.sin(th)]])

    def point(self, q, b, s):
        i = self.idx[b]
        s = np.asarray(s) / 1000.0
        if i < 0: return s, None
        x = q[3 * i:3 * i + 2] + self.A(q[3 * i + 2]) @ s
        J = np.zeros((2, 3 * self.n)); J[:, 3 * i:3 * i + 2] = np.eye(2); J[:, 3 * i + 2] = self.dA(q[3 * i + 2]) @ s
        return x, J

    def Pj(self, q, b, s):
        x, J = self.point(q, b, s)
        return x, (J if J is not None else np.zeros((2, 3 * self.n)))

    def joints_phi(self, q):
        phi, Jac = [], []
        for _, a, b, x, y in self.joints:
            pa, Ja = self.Pj(q, a, (x, y)); pb, Jb = self.Pj(q, b, (x, y))
            phi += list(pa - pb); Jac += list(Ja - Jb)
        return np.array(phi), np.array(Jac)

    def length(self, q, k):
        a = self.acts[k]
        pa, Ja = self.Pj(q, a['a'], a['pa']); pb, Jb = self.Pj(q, a['b'], a['pb'])
        d = pb - pa; L = np.linalg.norm(d)
        return L, (d / L) @ (Jb - Ja)

    def solve(self, L, q0):
        q = q0.copy()
        for _ in range(50):
            phi, J = self.joints_phi(q)
            r = list(phi); Jl = list(J)
            for k in range(len(self.acts)):
                Lk, g = self.length(q, k); r.append(Lk - L[k] / 1000); Jl.append(g)
            r, Jl = np.array(r), np.array(Jl)
            dq = np.linalg.solve(Jl, -r); q += dq
            if np.max(np.abs(dq)) < 1e-13: return q
        raise RuntimeError('position solver did not converge')

    def gravity(self, q):
        Q = np.zeros(3 * self.n)
        for i in range(self.n):
            th = q[3 * i + 2]
            Q[3 * i + 1] = -G * self.m[i]
            Q[3 * i + 2] = -G * (self.dA(th) @ self.S[i])[1]   # d/dtheta of -g * (A S)_y
        return Q

    def forces(self, q):
        _, J = self.joints_phi(q)
        cols = [J.T] + [self.length(q, k)[1][:, None] for k in range(len(self.acts))]
        A = np.hstack(cols)                      # [Phi_q^T | dL/dq ...]
        sol = np.linalg.solve(A, -self.gravity(q))
        return sol[-len(self.acts):]             # actuator forces, + pushes the pins apart

    def mass_matrix(self, q):
        M = np.zeros((3 * self.n, 3 * self.n))
        for i in range(self.n):
            th = q[3 * i + 2]; c = self.dA(th) @ self.S[i]
            M[3 * i:3 * i + 2, 3 * i:3 * i + 2] = self.m[i] * np.eye(2)
            M[3 * i:3 * i + 2, 3 * i + 2] = c; M[3 * i + 2, 3 * i:3 * i + 2] = c
            M[3 * i + 2, 3 * i + 2] = self.Jo[i]
        return M

    def reduced_mass(self, q, k):
        _, J = self.joints_phi(q)
        rows = [J] + [self.length(q, j)[1][None, :] for j in range(len(self.acts))]
        C = np.vstack(rows)
        rhs = np.zeros(C.shape[0]); rhs[J.shape[0] + k] = 1.0
        qd = np.linalg.solve(C, rhs)
        return qd @ self.mass_matrix(q) @ qd

    def k_hyd(self, k, L):
        a, h = self.acts[k], self.hyd
        beta = float(h.get('beta', 1e9))
        v0p = float(h.get('dead_volume_piston', h.get('dead_volume', 0.3))) * 1e-3
        v0r = float(h.get('dead_volume_rod', h.get('dead_volume', 0.3))) * 1e-3
        Ap = math.pi / 4 * (a['bore'] / 1000) ** 2; Ar = Ap - math.pi / 4 * (a['rod'] / 1000) ** 2
        x = (L - a['lmin']) / 1000; s = (a['lmax'] - a['lmin']) / 1000
        return beta * (Ap ** 2 / (v0p + Ap * x) + Ar ** 2 / (v0r + Ar * (s - x)))


def compare_mass_properties(R):
    """CylLoad level-1 mass properties (parametric B-spline integration and the 1.1 triangle fans)
    against the exact values written by gen_excavator.py."""
    ref = {r['assembly']: r for r in csv.DictReader(l for l in open(os.path.join(HERE, 'excavator_reference.csv'))
                                                    if not l.startswith('#'))}
    lines = ['assembly,mode,faces_bspline,faces_triangulated,volume_rel_err,centroid_err_mm,Jz_rel_err']
    worst = {}
    for mode, fn in (('parametric', 'step_level1.csv'), ('fan', 'step_level1_fan.csv')):
        for r in csv.DictReader(open(os.path.join(R, fn))):
            q = ref[r['name']]
            V = float(q['volume'])
            dv = (float(r['volume_mm3']) - V) / V
            dc = max(abs(float(r[a]) - float(q[b])) for a, b in (('xg', 'xc'), ('yg', 'yc'), ('zg', 'zc')))
            kz = (float(q['cxx']) + float(q['cyy'])) / V * 1e-6
            dk = (float(r['kz2_m2']) - kz) / kz
            lines.append(f"{r['name']},{mode},{r['faces_bspline_parametric']},{r['faces_triangulated']},{dv:.3e},{dc:.3e},{dk:.3e}")
            w = worst.setdefault(mode, [0.0, 0.0, 0.0])
            w[0], w[1], w[2] = max(w[0], abs(dv)), max(w[1], dc), max(w[2], abs(dk))
    open(os.path.join(R, 'check_mass_properties.csv'), 'w').write('\n'.join(lines) + '\n')
    return worst


def main():
    mech = Mech(os.path.join(HERE, 'excavator.mech'))
    R = os.path.join(HERE, 'results')
    worst = compare_mass_properties(R)
    q_ref = np.zeros(3 * mech.n)
    L_ref = []
    for a in mech.acts:
        L_ref.append(1000 * math.dist(np.array(a['pa']) / 1000, np.array(a['pb']) / 1000))
    out = open(os.path.join(R, 'check_independent.csv'), 'w')
    out.write('actuator,L_mm,F_cylload_N,F_independent_N,m_red_cylload_kg,m_red_independent_kg,f_n_cylload_Hz,f_n_independent_Hz\n')
    summary = []
    for k, a in enumerate(mech.acts):
        rows = list(csv.DictReader(open(os.path.join(R, f"sweep_{a['name']}.csv"))))
        eF = em = ef = 0.0; Fmax = max(abs(float(r['F_N'])) for r in rows)
        # continuation from the reference pose to the first point keeps the assembly branch
        q, Lc = q_ref.copy(), L_ref[k]
        for r in rows:
            L = list(L_ref); L[k] = float(r['L_mm'])
            n = max(1, int(abs(L[k] - Lc) / 5.0))
            for j in range(1, n + 1):
                Lj = list(L_ref); Lj[k] = Lc + (L[k] - Lc) * j / n
                q = mech.solve(Lj, q)
            Lc = L[k]
            F = mech.forces(q)[k]
            mr = mech.reduced_mass(q, k)
            fn = math.sqrt(mech.k_hyd(k, L[k]) / mr) / (2 * math.pi)
            Fc, mc, fc = float(r['F_N']), float(r['m_red_kg']), float(r['f_n_Hz'])
            eF = max(eF, abs(F - Fc) / Fmax); em = max(em, abs(mr - mc) / mc); ef = max(ef, abs(fn - fc) / fc)
            out.write(f"{a['name']},{L[k]:.6g},{Fc:.10g},{F:.10g},{mc:.10g},{mr:.10g},{fc:.10g},{fn:.10g}\n")
        summary.append((a['name'], len(rows), eF, em, ef))
    with open(os.path.join(R, 'check_independent.txt'), 'w') as f:
        f.write('Independent check (absolute coordinates + Lagrange multipliers, NumPy) against CylLoad sweeps\n')
        f.write('actuator  poses  max|dF|/max|F|  max|dm_r|/m_r  max|df_n|/f_n\n')
        for s in summary:
            f.write(f'{s[0]:8s} {s[1]:6d}  {s[2]:14.2e}  {s[3]:13.2e}  {s[4]:13.2e}\n')
        f.write('\nMass properties of the 13 level-1 assemblies against the exact values (largest errors)\n')
        f.write('B-spline faces      |dV|/V      |d centroid|, mm   |dJz|/Jz\n')
        for mode in ('parametric', 'fan'):
            w = worst[mode]
            f.write(f'{mode:18s}  {w[0]:9.2e}   {w[1]:12.2e}      {w[2]:9.2e}\n')
        f.write('(details in check_mass_properties.csv)\n')
    sys.stdout.write(open(os.path.join(R, 'check_independent.txt')).read())


if __name__ == '__main__':
    main()
