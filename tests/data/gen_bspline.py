#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Generate STEP test solids bounded by B-spline faces, with their exact mass properties.

* bump      - prism over a 80 x 60 mm rectangle whose top is a non-convex bicubic B-spline
              surface with interior knots (a doubly curved, saddle-shaped face);
* fillet    - bar of a 50 x 40 mm section with a convex 15 mm fillet (rational quadratic face);
* bracket   - L-section bar with a concave 12 mm fillet (rational quadratic face).
The fillet arcs of the end caps are rational B-spline curves. Reference values are
computed here independently of CylLoad: by Gauss quadrature over the B-spline parameter
domain (bump) and from a densely sampled section polygon (bars).
Usage: python3 gen_bspline.py OUT_DIR
"""
import math, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gen_revolved import W, rot  # noqa: E402

# ------------------------------------------------------------------ B-spline helpers
def knots_clamped(p, n_ctrl):
    inner = n_ctrl - p - 1
    return [0.0] * (p + 1) + [(i + 1) / (inner + 1) for i in range(inner)] + [1.0] * (p + 1)


def basis(K, p, t):
    n = len(K) - p - 1
    N = [1.0 if (K[i] <= t < K[i + 1]) or (t == K[-1] and K[i] < K[i + 1] == K[-1]) else 0.0 for i in range(len(K) - 1)]
    for d in range(1, p + 1):
        M = []
        for i in range(len(K) - 1 - d):
            a = (t - K[i]) / (K[i + d] - K[i]) * N[i] if K[i + d] > K[i] else 0.0
            b = (K[i + d + 1] - t) / (K[i + d + 1] - K[i + 1]) * N[i + 1] if K[i + d + 1] > K[i + 1] else 0.0
            M.append(a + b)
        N = M
    return N[:n]


def greville(K, p, n):
    return [sum(K[i + 1:i + p + 1]) / p for i in range(n)]


GX = [-0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, -0.1488743389816312, 0.1488743389816312, 0.4333953941292472, 0.6794095682990244, 0.8650633666889845, 0.9739065285171717]
GW = [0.0666713443086881, 0.1494513491505804, 0.2190863625159820, 0.2692667193099965, 0.2955242247147528, 0.2955242247147528, 0.2692667193099965, 0.2190863625159820, 0.1494513491505804, 0.0666713443086881]

# ------------------------------------------------------------------ bump
BUMP = dict(a=80.0, b=60.0, pu=3, pv=2,
            z=[[20, 26, 22, 30], [24, 50, 8, 26], [30, 5, 44, 24], [22, 34, 12, 28], [25, 21, 30, 20]])


def bump_surface():
    d = BUMP
    nu, nv = len(d['z']), len(d['z'][0])
    KU, KV = knots_clamped(d['pu'], nu), knots_clamped(d['pv'], nv)
    gu, gv = greville(KU, d['pu'], nu), greville(KV, d['pv'], nv)
    P = [[(d['a'] * gu[i], d['b'] * gv[j], d['z'][i][j]) for j in range(nv)] for i in range(nu)]
    return P, KU, KV


def bump_props():
    P, KU, KV = bump_surface()
    a, b, pu, pv = BUMP['a'], BUMP['b'], BUMP['pu'], BUMP['pv']
    bu, bv = sorted(set(KU)), sorted(set(KV))
    I = {}
    keys = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (0, 1, 1), (1, 0, 1)]
    for k in keys: I[k] = 0.0
    for i in range(len(bu) - 1):
        for j in range(len(bv) - 1):
            for xu, wu in zip(GX, GW):
                u = bu[i] + (bu[i + 1] - bu[i]) * (xu + 1) / 2
                Nu = basis(KU, pu, u)
                for xv, wv in zip(GX, GW):
                    v = bv[j] + (bv[j + 1] - bv[j]) * (xv + 1) / 2
                    Nv = basis(KV, pv, v)
                    z = sum(Nu[r] * Nv[s] * P[r][s][2] for r in range(len(Nu)) for s in range(len(Nv)))
                    x, y = a * u, b * v
                    wt = wu * wv * (bu[i + 1] - bu[i]) * (bv[j + 1] - bv[j]) / 4 * a * b
                    for (p, q, r) in keys:  # int x^p y^q z^r dV = int x^p y^q z^(r+1)/(r+1) dA
                        I[(p, q, r)] += wt * x ** p * y ** q * z ** (r + 1) / (r + 1)
    return I


# ------------------------------------------------------------------ bars (extruded sections in the (y, z) plane)
S2 = math.sqrt(0.5)
BARS = {
    'fillet': dict(L=100.0, prof=[('L', (0, 0), (50, 0)), ('L', (50, 0), (50, 25)),
                                  ('A', (50, 25), (35, 40), (35, 25)), ('L', (35, 40), (0, 40)), ('L', (0, 40), (0, 0))]),
    'bracket': dict(L=90.0, prof=[('L', (0, 0), (60, 0)), ('L', (60, 0), (60, 10)), ('L', (60, 10), (22, 10)),
                                  ('A', (22, 10), (10, 22), (22, 22)), ('L', (10, 22), (10, 50)), ('L', (10, 50), (0, 50)),
                                  ('L', (0, 50), (0, 0))]),
}


def arc_ctrl(s):
    (y1, z1), (y2, z2), (cy, cz) = s[1], s[2], s[3]
    # 90 degree arc: middle control point at the intersection of the end tangents
    m = (y1 + y2 - cy, z1 + z2 - cz)
    return [(y1, z1), m, (y2, z2)], [1.0, S2, 1.0]


def arc_pt(s, t):
    P, w = arc_ctrl(s)
    B = [(1 - t) ** 2 * w[0], 2 * t * (1 - t) * w[1], t * t * w[2]]
    den = sum(B)
    return (sum(B[i] * P[i][0] for i in range(3)) / den, sum(B[i] * P[i][1] for i in range(3)) / den)


def bar_props(nm):
    d = BARS[nm]
    L = d['L']
    pts = []
    for s in d['prof']:
        if s[0] == 'L':
            pts += [(s[1][0] + (s[2][0] - s[1][0]) * k / 40, s[1][1] + (s[2][1] - s[1][1]) * k / 40) for k in range(40)]
        else:
            pts += [arc_pt(s, k / 40000) for k in range(40000)]

    def A(p, q):  # int y^p z^q dA by Green: oint y^(p+1) z^q/(p+1) dz
        t = 0.0
        for i in range(len(pts)):
            (y1, z1), (y2, z2) = pts[i], pts[(i + 1) % len(pts)]
            for x, w in zip(GX, GW):
                tt = (x + 1) / 2
                y, z = y1 + (y2 - y1) * tt, z1 + (z2 - z1) * tt
                t += 0.5 * w * y ** (p + 1) * z ** q / (p + 1) * (z2 - z1)
        return t
    I = {}
    for (p, q, r) in [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (0, 1, 1), (1, 0, 1)]:
        I[(p, q, r)] = L ** (p + 1) / (p + 1) * A(q, r)  # x in [0, L]
    return I


# ------------------------------------------------------------------ STEP topology
class Topo:
    def __init__(self, w):
        self.w, self.verts, self.edges = w, {}, {}

    def vertex(self, p):
        k = tuple(round(c, 10) for c in p)
        if k not in self.verts:
            self.verts[k] = self.w.add(f"VERTEX_POINT('',#{self.w.pt(p)})")
        return self.verts[k]

    def line(self, a, b):
        key = (tuple(round(c, 10) for c in a), tuple(round(c, 10) for c in b))
        if key[::-1] in self.edges:
            return self.edges[key[::-1]], False
        if key not in self.edges:
            d = [b[i] - a[i] for i in range(3)]
            Ln = math.sqrt(sum(x * x for x in d))
            vec = self.w.add(f"VECTOR('',#{self.w.dr([x / Ln for x in d])},{Ln:.15g})")
            ln = self.w.add(f"LINE('',#{self.w.pt(a)},#{vec})")
            self.edges[key] = self.w.add(f"EDGE_CURVE('',#{self.vertex(a)},#{self.vertex(b)},#{ln},.T.)")
        return self.edges[key], True

    def bspline(self, pts, p, weights=None, tag=None):
        a, b = pts[0], pts[-1]
        key = ('c', tag)
        if key not in self.edges:
            K = knots_clamped(p, len(pts))
            ks = sorted(set(K))
            mult = [K.count(k) for k in ks]
            cp = ','.join('#%d' % self.w.pt(q) for q in pts)
            body = (f"B_SPLINE_CURVE({p},({cp}),.UNSPECIFIED.,.F.,.F.)"
                    f"B_SPLINE_CURVE_WITH_KNOTS(({','.join(map(str, mult))}),({','.join('%.15g' % k for k in ks)}),.UNSPECIFIED.)")
            if weights:
                c = self.w.add(f"(BOUNDED_CURVE(){body}CURVE()GEOMETRIC_REPRESENTATION_ITEM()"
                               f"RATIONAL_B_SPLINE_CURVE(({','.join('%.15g' % x for x in weights)}))REPRESENTATION_ITEM(''))")
            else:
                c = self.w.add(f"B_SPLINE_CURVE_WITH_KNOTS('',{p},({cp}),.UNSPECIFIED.,.F.,.F.,({','.join(map(str, mult))}),"
                               f"({','.join('%.15g' % k for k in ks)}),.UNSPECIFIED.)")
            self.edges[key] = (self.w.add(f"EDGE_CURVE('',#{self.vertex(a)},#{self.vertex(b)},#{c},.T.)"), a)
        e, start = self.edges[key]
        fwd = all(abs(start[i] - pts[0][i]) < 1e-9 for i in range(3))
        return e, fwd

    def face(self, oriented, surf, same):
        ids = [self.w.add(f"ORIENTED_EDGE('',*,*,#{e},{'.T.' if f else '.F.'})") for e, f in oriented]
        el = self.w.add(f"EDGE_LOOP('',({','.join('#%d' % i for i in ids)}))")
        fb = self.w.add(f"FACE_OUTER_BOUND('',#{el},.T.)")
        return self.w.add(f"ADVANCED_FACE('',(#{fb}),#{surf},{'.T.' if same else '.F.'})")

    def plane(self, o, n, x):
        return self.w.add(f"PLANE('',#{self.w.ax(o, n, x)})")


def surface_entity(w, P, KU, KV, pu, pv, Wt=None):
    ksu, ksv = sorted(set(KU)), sorted(set(KV))
    rows = ','.join('(' + ','.join('#%d' % w.pt(q) for q in row) + ')' for row in P)
    mu, mv = [KU.count(k) for k in ksu], [KV.count(k) for k in ksv]
    kn = (f"({','.join(map(str, mu))}),({','.join(map(str, mv))}),({','.join('%.15g' % k for k in ksu)}),"
          f"({','.join('%.15g' % k for k in ksv)}),.UNSPECIFIED.")
    if Wt is None:
        return w.add(f"B_SPLINE_SURFACE_WITH_KNOTS('',{pu},{pv},({rows}),.UNSPECIFIED.,.F.,.F.,.F.,{kn})")
    wt = ','.join('(' + ','.join('%.15g' % x for x in row) + ')' for row in Wt)
    return w.add(f"(BOUNDED_SURFACE()B_SPLINE_SURFACE({pu},{pv},({rows}),.UNSPECIFIED.,.F.,.F.,.F.)"
                 f"B_SPLINE_SURFACE_WITH_KNOTS({kn})GEOMETRIC_REPRESENTATION_ITEM()RATIONAL_B_SPLINE_SURFACE(({wt}))"
                 f"REPRESENTATION_ITEM('')SURFACE())")


def build_bump(w):
    T = Topo(w)
    P, KU, KV = bump_surface()
    a, b = BUMP['a'], BUMP['b']
    nu, nv = len(P), len(P[0])
    B = lambda x, y: (x, y, 0.0)  # noqa: E731
    top = surface_entity(w, P, KU, KV, BUMP['pu'], BUMP['pv'])
    c_v0 = [P[i][0] for i in range(nu)]           # y = 0, along u (x)
    c_v1 = [P[i][nv - 1] for i in range(nu)]      # y = b
    c_u0 = [P[0][j] for j in range(nv)]           # x = 0, along v (y)
    c_u1 = [P[nu - 1][j] for j in range(nv)]      # x = a
    pu, pv = BUMP['pu'], BUMP['pv']
    e_v0, e_v1 = T.bspline(c_v0, pu, tag='v0'), T.bspline(c_v1, pu, tag='v1')
    e_u0, e_u1 = T.bspline(c_u0, pv, tag='u0'), T.bspline(c_u1, pv, tag='u1')
    t00, t10, t11, t01 = P[0][0], P[nu - 1][0], P[nu - 1][nv - 1], P[0][nv - 1]
    inv = lambda e: (e[0], not e[1])  # noqa: E731
    faces = []
    # top (normal +z = Su x Sv): CCW from above t00 -> t10 -> t11 -> t01
    faces.append(T.face([e_v0, e_u1, inv(e_v1), inv(e_u0)], top, True))
    # bottom (normal -z): (0,0) -> (0,b) -> (a,b) -> (a,0)
    faces.append(T.face([T.line(B(0, 0), B(0, b)), T.line(B(0, b), B(a, b)), T.line(B(a, b), B(a, 0)), T.line(B(a, 0), B(0, 0))],
                        T.plane((0, 0, 0), (0, 0, -1), (1, 0, 0)), True))
    # y = 0 (normal -y): B(0,0) -> B(a,0) -> t10 -> (curve back) t00 -> B(0,0)
    faces.append(T.face([T.line(B(0, 0), B(a, 0)), T.line(B(a, 0), t10), inv(e_v0), T.line(t00, B(0, 0))],
                        T.plane((0, 0, 0), (0, -1, 0), (1, 0, 0)), True))
    # y = b (normal +y): B(a,b) -> B(0,b) -> t01 -> t11 -> B(a,b)
    faces.append(T.face([T.line(B(a, b), B(0, b)), T.line(B(0, b), t01), e_v1, T.line(t11, B(a, b))],
                        T.plane((0, b, 0), (0, 1, 0), (1, 0, 0)), True))
    # x = 0 (normal -x): B(0,b) -> B(0,0) -> t00 -> t01 -> B(0,b)
    faces.append(T.face([T.line(B(0, b), B(0, 0)), T.line(B(0, 0), t00), e_u0, T.line(t01, B(0, b))],
                        T.plane((0, 0, 0), (-1, 0, 0), (0, 1, 0)), True))
    # x = a (normal +x): B(a,0) -> B(a,b) -> t11 -> t10 -> B(a,0)
    faces.append(T.face([T.line(B(a, 0), B(a, b)), T.line(B(a, b), t11), inv(e_u1), T.line(t10, B(a, 0))],
                        T.plane((a, 0, 0), (1, 0, 0), (0, 1, 0)), True))
    shell = w.add(f"CLOSED_SHELL('',({','.join('#%d' % f for f in faces)}))")
    return w.add(f"MANIFOLD_SOLID_BREP('',#{shell})")


def build_bar(w, nm):
    T = Topo(w)
    d = BARS[nm]
    L, prof = d['L'], d['prof']
    X = lambda x, q: (x, q[0], q[1])  # noqa: E731
    inv = lambda e: (e[0], not e[1])  # noqa: E731

    def seg_edge(s, x, k):
        if s[0] == 'L':
            return T.line(X(x, s[1]), X(x, s[2]))
        cp, wt = arc_ctrl(s)
        return T.bspline([X(x, q) for q in cp], 2, wt, tag=(nm, k, x))
    faces = []
    for k, s in enumerate(prof):
        e0, e1 = seg_edge(s, 0.0, k), seg_edge(s, L, k)
        loop = [e0, T.line(X(0.0, s[2]), X(L, s[2])), inv(e1), T.line(X(L, s[1]), X(0.0, s[1]))]
        if s[0] == 'L':
            (y1, z1), (y2, z2) = s[1], s[2]
            n = (0.0, z2 - z1, -(y2 - y1))
            ln = math.hypot(n[1], n[2])
            surf = T.plane(X(0.0, s[1]), (0.0, n[1] / ln, n[2] / ln), (1.0, 0.0, 0.0))
            same = True
        else:
            cp, wt = arc_ctrl(s)
            P = [[X(0.0, q) for q in cp], [X(L, q) for q in cp]]
            surf = surface_entity(w, P, [0, 0, 1, 1], [0, 0, 0, 1, 1, 1], 1, 2, [wt, wt])
            # Su x Sv = x x t (t along the profile) is the inward normal of a CCW profile
            same = False
        faces.append(T.face(loop, surf, same))
    # end caps: x = L (normal +x, profile order), x = 0 (normal -x, reversed)
    faces.append(T.face([seg_edge(s, L, k) for k, s in enumerate(prof)], T.plane((L, 0, 0), (1, 0, 0), (0, 1, 0)), True))
    faces.append(T.face([inv(seg_edge(s, 0.0, k)) for k, s in reversed(list(enumerate(prof)))],
                        T.plane((0, 0, 0), (-1, 0, 0), (0, 1, 0)), True))
    shell = w.add(f"CLOSED_SHELL('',({','.join('#%d' % f for f in faces)}))")
    return w.add(f"MANIFOLD_SOLID_BREP('',#{shell})")


PLACEMENT = {'bump': ((0.2, 1, 0.4), 35.0, (300.0, -20.0, 10.0)),
             'fillet': ((1, 0, 0), -20.0, (-150.0, 60.0, 0.0)),
             'bracket': ((0, 0, 1), 60.0, (40.0, 250.0, -30.0))}


def to_global(I, nm):
    V = I[(0, 0, 0)]
    c = [I[(1, 0, 0)] / V, I[(0, 1, 0)] / V, I[(0, 0, 1)] / V]
    S = [[I[(2, 0, 0)], I[(1, 1, 0)], I[(1, 0, 1)]], [I[(1, 1, 0)], I[(0, 2, 0)], I[(0, 1, 1)]], [I[(1, 0, 1)], I[(0, 1, 1)], I[(0, 0, 2)]]]
    C = [[S[i][j] - V * c[i] * c[j] for j in range(3)] for i in range(3)]
    axis, ang, t = PLACEMENT[nm]
    R = rot(axis, ang)
    cg = [sum(R[i][k] * c[k] for k in range(3)) + t[i] for i in range(3)]
    Cg = [[sum(R[i][k] * C[k][l] * R[j][l] for k in range(3) for l in range(3)) for j in range(3)] for i in range(3)]
    return V, cg, Cg


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    w = W()
    ctx = w.add("APPLICATION_CONTEXT('configuration controlled 3D designs')")
    pc = w.add(f"PRODUCT_CONTEXT('',#{ctx},'mechanical')")
    pdc = w.add(f"PRODUCT_DEFINITION_CONTEXT('part definition',#{ctx},'design')")
    unit = w.add("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))")
    rc = w.add(f"(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#{unit}))REPRESENTATION_CONTEXT('',''))")

    def product(name):
        p = w.add(f"PRODUCT('{name}','{name}','',(#{pc}))")
        f = w.add(f"PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE('','',#{p},.NOT_KNOWN.)")
        d = w.add(f"PRODUCT_DEFINITION('design','',#{f},#{pdc})")
        s = w.add(f"PRODUCT_DEFINITION_SHAPE('','',#{d})")
        return d, s
    top_pd, top_pds = product('BSPLINE_TEST')
    top_ax = w.ax((0, 0, 0), (0, 0, 1), (1, 0, 0))
    top_sr = w.add(f"SHAPE_REPRESENTATION('',(#{top_ax}),#{rc})")
    w.add(f"SHAPE_DEFINITION_REPRESENTATION(#{top_pds},#{top_sr})")
    names = ['bump', 'fillet', 'bracket']
    for i, nm in enumerate(names):
        pd, pds = product(nm)
        own = w.ax((0, 0, 0), (0, 0, 1), (1, 0, 0))
        sr = w.add(f"SHAPE_REPRESENTATION('',(#{own}),#{rc})")
        w.add(f"SHAPE_DEFINITION_REPRESENTATION(#{pds},#{sr})")
        solid = build_bump(w) if nm == 'bump' else build_bar(w, nm)
        ab = w.add(f"ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{solid},#{own}),#{rc})")
        w.add(f"SHAPE_REPRESENTATION_RELATIONSHIP('','',#{sr},#{ab})")
        axis, ang, t = PLACEMENT[nm]
        R = rot(axis, ang)
        place = w.ax(t, (R[0][2], R[1][2], R[2][2]), (R[0][0], R[1][0], R[2][0]))
        nauo = w.add(f"NEXT_ASSEMBLY_USAGE_OCCURRENCE('{i}','{nm}:1','',#{top_pd},#{pd},$)")
        npds = w.add(f"PRODUCT_DEFINITION_SHAPE('','',#{nauo})")
        idt = w.add(f"ITEM_DEFINED_TRANSFORMATION('','',#{place},#{own})")
        rr = w.add(f"(REPRESENTATION_RELATIONSHIP('','',#{sr},#{top_sr})REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#{idt})SHAPE_REPRESENTATION_RELATIONSHIP())")
        w.add(f"CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#{rr},#{npds})")
    with open(os.path.join(out, 'bspline_solids.step'), 'w') as f:
        f.write("ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('cylload test solids with B-spline faces'),'2;1');\n"
                "FILE_NAME('bspline_solids.step','2026-11-10',(''),(''),'gen_bspline.py','','');\n"
                "FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));\nENDSEC;\nDATA;\n")
        f.write('\n'.join(w.lines))
        f.write("\nENDSEC;\nEND-ISO-10303-21;\n")
    with open(os.path.join(out, 'bspline_reference.csv'), 'w') as f:
        f.write('name,volume,xc,yc,zc,cxx,cxy,cxz,cyy,cyz,czz\n')
        for nm in names:
            I = bump_props() if nm == 'bump' else bar_props(nm)
            V, cg, Cg = to_global(I, nm)
            f.write(f"{nm},{V:.12g},{cg[0]:.12g},{cg[1]:.12g},{cg[2]:.12g},{Cg[0][0]:.12g},{Cg[0][1]:.12g},{Cg[0][2]:.12g},"
                    f"{Cg[1][1]:.12g},{Cg[1][2]:.12g},{Cg[2][2]:.12g}\n")


if __name__ == '__main__':
    main()
