# -*- coding: utf-8 -*-
"""Comparison of CylLoad with Exudyn (J. Gerstmayr, Univ. Innsbruck) for the open excavator example.

The mechanism of excavator.mech is rebuilt in Exudyn from planar rigid bodies (NodeRigidBody2D,
ObjectRigidBody2D), revolute joints (ObjectJointRevolute2D) and cylinder lengths as distance
constraints (ObjectConnectorDistance). Gravity acts through LoadMassProportional. Exudyn's static
solver gives the pose and the constraint (cylinder) forces; the reduced mass follows from the
Exudyn pose derivative and mass matrix. Results are compared with the CylLoad sweeps.
Usage (after run.sh): python -m pip install exudyn numpy; python examples/excavator/check_exudyn.py
Writes results/check_exudyn.txt and results/check_exudyn.csv.
"""
import csv, math, os, sys, traceback

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, 'results') if os.path.isdir(os.path.join(HERE, 'results')) else HERE
LOG = open(os.path.join(RES, 'check_exudyn.txt'), 'w', encoding='utf-8')


def log(*a):
    s = ' '.join(str(x) for x in a)
    print(s)
    LOG.write(s + '\n')
    LOG.flush()


try:
    import numpy as np
    import exudyn as exu
    from exudyn.itemInterface import (NodeRigidBody2D, ObjectRigidBody2D, ObjectGround, MarkerBodyPosition,
                                      MarkerBodyMass, ObjectJointRevolute2D, ObjectConnectorDistance,
                                      LoadMassProportional)
    log('exudyn', exu.GetVersionString(), 'numpy', np.__version__, 'python', sys.version.split()[0])
except Exception:
    log('IMPORT ERROR'); log(traceback.format_exc()); sys.exit(1)

# ------------------------------------------------------------------ read the CylLoad model
bodies, joints, acts, hyd = [], [], [], {}
for line in open(os.path.join(HERE, 'excavator.mech'), encoding='utf-8'):
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
items = [(r['body'], float(r['mass_kg']), float(r['x_mm']) / 1000, float(r['y_mm']) / 1000, float(r['J_kgm2']))
         for r in csv.DictReader(open(os.path.join(HERE, 'excavator_masses.csv'), encoding='utf-8'))]
for a in acts:  # half of each cylinder mass at each pin (same modelling assumption as CylLoad)
    items.append((a['a'], a['mass'] / 2, a['pa'][0] / 1000, a['pa'][1] / 1000, 0.0))
    items.append((a['b'], a['mass'] / 2, a['pb'][0] / 1000, a['pb'][1] / 1000, 0.0))
props = {}
for b in bodies:
    its = [i for i in items if i[0] == b]
    m = sum(i[1] for i in its)
    cx = sum(i[1] * i[2] for i in its) / m
    cy = sum(i[1] * i[3] for i in its) / m
    J = sum(i[4] + i[1] * ((i[2] - cx) ** 2 + (i[3] - cy) ** 2) for i in its)
    props[b] = (m, cx, cy, J)
    log(f'body {b:8s} m = {m:10.4f} kg  CoG = ({cx*1000:9.3f}, {cy*1000:9.3f}) mm  J = {J:10.5f} kg m2')

# ------------------------------------------------------------------ Exudyn model
SC = exu.SystemContainer()
mbs = SC.AddSystem()
oGround = mbs.AddObject(ObjectGround(referencePosition=[0, 0, 0]))
node, obj = {}, {}
for b in bodies:
    m, cx, cy, J = props[b]
    node[b] = mbs.AddNode(NodeRigidBody2D(referenceCoordinates=[cx, cy, 0], initialCoordinates=[0, 0, 0],
                                          initialVelocities=[0, 0, 0]))
    obj[b] = mbs.AddObject(ObjectRigidBody2D(physicsMass=m, physicsInertia=J, nodeNumber=node[b]))
    mM = mbs.AddMarker(MarkerBodyMass(bodyNumber=obj[b]))
    mbs.AddLoad(LoadMassProportional(markerNumber=mM, loadVector=[0, -9.81, 0]))


def marker(b, x, y):
    x, y = x / 1000, y / 1000
    if b == 'ground':
        return mbs.AddMarker(MarkerBodyPosition(bodyNumber=oGround, localPosition=[x, y, 0]))
    m, cx, cy, J = props[b]
    return mbs.AddMarker(MarkerBodyPosition(bodyNumber=obj[b], localPosition=[x - cx, y - cy, 0]))


jmark = {}
for name, a, b, x, y in joints:
    jmark[name] = (marker(a, x, y), marker(b, x, y))
    mbs.AddObject(ObjectJointRevolute2D(markerNumbers=list(jmark[name])))
dist = []
for a in acts:
    L0 = math.dist(a['pa'], a['pb']) / 1000
    a['L0'] = L0
    dist.append(mbs.AddObject(ObjectConnectorDistance(markerNumbers=[marker(a['a'], *a['pa']), marker(a['b'], *a['pb'])],
                                                      distance=L0)))
mbs.Assemble()
ss = exu.SimulationSettings()
ss.staticSolver.newton.relativeTolerance = 1e-12
ss.staticSolver.newton.absoluteTolerance = 1e-10
ss.staticSolver.newton.maxIterations = 50
ss.staticSolver.verboseMode = 0
ss.linearSolverType = exu.LinearSolverType.EXUdense
ss.displayStatistics = False
ss.displayComputationTime = False


def solve(L):
    for k, o in enumerate(dist):
        mbs.SetObjectParameter(o, 'distance', L[k])
    ok = exu.SolveStatic(mbs, ss, updateInitialValues=True)
    if not ok:
        raise RuntimeError('Exudyn static solver failed at ' + str(L))
    return np.array(mbs.systemData.GetODE2Coordinates())


def force(k):
    f = mbs.GetObjectOutput(dist[k], exu.OutputVariableType.Force)
    return float(np.atleast_1d(f)[0])


def mass_matrix():
    try:
        return np.array(mbs.ComputeMassMatrix())
    except Exception:
        M = np.zeros((3 * len(bodies),) * 2)
        for i, b in enumerate(bodies):
            m, cx, cy, J = props[b]
            M[3 * i, 3 * i] = M[3 * i + 1, 3 * i + 1] = m
            M[3 * i + 2, 3 * i + 2] = J
        return M


def k_hyd(a, Lmm):
    beta = float(hyd.get('beta', 1e9))
    v0p = float(hyd.get('dead_volume_piston', hyd.get('dead_volume', 0.3))) * 1e-3
    v0r = float(hyd.get('dead_volume_rod', hyd.get('dead_volume', 0.3))) * 1e-3
    Ap = math.pi / 4 * (a['bore'] / 1000) ** 2
    Ar = Ap - math.pi / 4 * (a['rod'] / 1000) ** 2
    x = (Lmm - a['lmin']) / 1000
    s = (a['lmax'] - a['lmin']) / 1000
    return beta * (Ap ** 2 / (v0p + Ap * x) + Ar ** 2 / (v0r + Ar * (s - x)))


# sign convention of the Exudyn distance force, from the reference pose
q = solve([a['L0'] for a in acts])
out = open(os.path.join(RES, 'check_exudyn.csv'), 'w', encoding='utf-8', newline='')
out.write('actuator,L_mm,F_cylload_N,F_exudyn_N,m_red_cylload_kg,m_red_exudyn_kg,f_n_cylload_Hz,f_n_exudyn_Hz\n')
M = mass_matrix()
log('mass matrix from', 'Exudyn ComputeMassMatrix' if hasattr(mbs, 'ComputeMassMatrix') else 'body properties')
summary = []
for k, a in enumerate(acts):
    rows = list(csv.DictReader(open(os.path.join(RES, f"sweep_{a['name']}.csv"), encoding='utf-8')))
    Fmax = max(abs(float(r['F_N'])) for r in rows)
    sgn = None
    eF = em = ef = 0.0
    Lc = a['L0'] * 1000
    # start every sweep from the model pose (otherwise Newton may jump to another assembly branch)
    nq = len(mbs.systemData.GetODE2Coordinates())
    mbs.systemData.SetODE2Coordinates(np.zeros(nq), configuration=exu.ConfigurationType.Initial)
    q = solve([b['L0'] for b in acts])
    for r in rows:
        Lt = float(r['L_mm'])
        n = max(1, int(abs(Lt - Lc) / 5))
        for j in range(1, n + 1):  # continuation in 5 mm steps keeps the assembly branch
            L = [b['L0'] for b in acts]
            L[k] = (Lc + (Lt - Lc) * j / n) / 1000
            q = solve(L)
        Lc = Lt
        Fe = force(k)
        Fc = float(r['F_N'])
        if sgn is None:
            sgn = 1.0 if abs(Fe - Fc) < abs(Fe + Fc) else -1.0
            log(f"{a['name']}: Exudyn distance force sign relative to CylLoad (push positive): {sgn:+.0f}")
        Fe *= sgn
        # pose derivative by central differences of Exudyn static solutions (other cylinders locked)
        h = 1e-5
        L = [b['L0'] for b in acts]
        L[k] = Lt / 1000 + h; qp = solve(L)
        L[k] = Lt / 1000 - h; qm = solve(L)
        L[k] = Lt / 1000; q = solve(L)
        qd = (qp - qm) / (2 * h)
        mr = float(qd @ M @ qd)
        fn = math.sqrt(k_hyd(a, Lt) / mr) / (2 * math.pi)
        mc, fc = float(r['m_red_kg']), float(r['f_n_Hz'])
        eF = max(eF, abs(Fe - Fc) / Fmax)
        em = max(em, abs(mr - mc) / mc)
        ef = max(ef, abs(fn - fc) / fc)
        out.write(f"{a['name']},{Lt:.6g},{Fc:.10g},{Fe:.10g},{mc:.10g},{mr:.10g},{fc:.10g},{fn:.10g}\n")
        # diagnostics: body rotations and joint residuals
        rot = [math.degrees(q[3 * i + 2]) for i in range(len(bodies))]
        rc = [float(r['rot_' + bb + '_deg']) for bb in bodies]
        res = 0.0
        for jn, ja, jb, x, y in joints:
            pa = np.array(mbs.GetMarkerOutput(jmark[jn][0], exu.OutputVariableType.Position))[:2]
            pb = np.array(mbs.GetMarkerOutput(jmark[jn][1], exu.OutputVariableType.Position))[:2]
            res = max(res, float(np.linalg.norm(pa - pb)))
        drot = max(abs(x - y) for x, y in zip(rot, rc))
        if abs(Lt - round(Lt / 50) * 50) < 1e-6:
            log(f"  {a['name']} L={Lt:7.1f}  rot Exudyn " + ' '.join(f'{x:8.3f}' for x in rot) +
                "  CylLoad " + ' '.join(f'{x:8.3f}' for x in rc) + f"  max drot {drot:.2e}  joint residual {res:.1e} m")
    summary.append((a['name'], len(rows), eF, em, ef))
out.close()
log('\nCylLoad vs Exudyn, excavator example')
log('actuator  poses  max|dF|/max|F|  max|dm_r|/m_r  max|df_n|/f_n')
for s in summary:
    log(f'{s[0]:8s} {s[1]:6d}  {s[2]:14.2e}  {s[3]:13.2e}  {s[4]:13.2e}')
log('bodies (rotation order):', ' '.join(bodies))
log('GOTOWE')
