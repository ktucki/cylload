#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Generate STEP AP203 test solids of revolution and their exact mass properties.

Each solid is a closed profile in the (r, z) half-plane revolved about the z axis and
split into two half faces (u in [0, pi] and [pi, 2 pi]) as CAD systems do. Profile
segments are lines or circular arcs, which produce planes, cylinders, cones (also to
the apex), spheres (to the poles) and tori, including the "lemon" torus with a
negative major radius used for barrel rollers. The solids are placed in a small
assembly with a rotation and a translation so that the placement code is tested too.

Reference values come from Pappus-type integrals over the profile, evaluated with a
densely sampled polygon (independent of the STEP reader).
Usage: python3 gen_revolved.py OUT_DIR
"""
import math, sys, os

# ---------------------------------------------------------------- profiles (CCW in (r, z))
# segment: ('L', (r1, z1), (r2, z2)) or ('A', (r1, z1), (r2, z2), (cr, cz), ccw)
PROFILES = {
    'pin': [  # plane, cylinder, torus fillet, cone chamfer, plane; touches the axis
        ('L', (0, 0), (10, 0)), ('L', (10, 0), (10, 20)),
        ('A', (10, 20), (8, 22), (8, 20), True),
        ('L', (8, 22), (5, 25)), ('L', (5, 25), (0, 25)), ('L', (0, 25), (0, 0))],
    'tube': [  # annulus ends, outer and inner cylinders (inner face reversed), chamfer cone
        ('L', (6, 0), (12, 0)), ('L', (12, 0), (12, 28)), ('L', (12, 28), (9, 31)),
        ('L', (9, 31), (6, 31)), ('L', (6, 31), (6, 0))],
    'cone_tip': [  # cone that ends in its apex
        ('L', (0, 0), (7, 0)), ('L', (7, 0), (0, 9)), ('L', (0, 9), (0, 0))],
    'ball': [  # sphere made of two hemispheres, each face touching one pole
        ('A', (0, -6), (6, 0), (0, 0), True), ('A', (6, 0), (0, 6), (0, 0), True), ('L', (0, 6), (0, -6))],
    'barrel': [  # spherical-roller barrel: lemon torus with a negative major radius
        ('L', (0, -4), (4.5, -4)),
        ('A', (4.5, -4), (4.5, 4), (-30.0, 0.0), True),
        ('L', (4.5, 4), (0, 4)), ('L', (0, 4), (0, -4))],
}
# rigid placement of each solid in the assembly: rotation about (axis) by angle, then translation
PLACEMENT = {
    'pin': ((1, 0, 0), 30.0, (100.0, -50.0, 20.0)),
    'tube': ((0, 1, 0), -90.0, (-200.0, 10.0, 0.0)),
    'cone_tip': ((1, 1, 0), 45.0, (5.0, 5.0, 5.0)),
    'ball': ((0, 0, 1), 10.0, (0.0, 300.0, -40.0)),
    'barrel': ((0.3, 0.2, 0.9), 70.0, (-35.0, 40.0, 7.0)),
}


def arc_pts(s, n):
    (r1, z1), (r2, z2), (cr, cz), ccw = s[1], s[2], s[3], s[4]
    R = math.hypot(r1 - cr, z1 - cz)
    a1, a2 = math.atan2(z1 - cz, r1 - cr), math.atan2(z2 - cz, r2 - cr)
    d = (a2 - a1) % (2 * math.pi) if ccw else -((a1 - a2) % (2 * math.pi))
    return [(cr + R * math.cos(a1 + d * k / n), cz + R * math.sin(a1 + d * k / n)) for k in range(n + 1)], R, a1, d


def polygon(prof, n=40000):
    pts = []
    for s in prof:
        if s[0] == 'L':
            pts += [(s[1][0] + (s[2][0] - s[1][0]) * k / 50, s[1][1] + (s[2][1] - s[1][1]) * k / 50) for k in range(50)]
        else:
            pts += arc_pts(s, n)[0][:-1]
    return pts


def area_int(pts, a, b):
    """int_D r^a z^b dA via Green: oint r^(a+1) z^b / (a+1) dz, Gauss on each polygon edge"""
    gx = [-0.906179845938664, -0.538469310105683, 0.0, 0.538469310105683, 0.906179845938664]
    gw = [0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189]
    tot = 0.0
    for i in range(len(pts)):
        (r1, z1), (r2, z2) = pts[i], pts[(i + 1) % len(pts)]
        for x, w in zip(gx, gw):
            t = 0.5 * (x + 1)
            r, z = r1 + (r2 - r1) * t, z1 + (z2 - z1) * t
            tot += 0.5 * w * r ** (a + 1) * z ** b / (a + 1) * (z2 - z1)
    return tot


def props(prof):
    pts = polygon(prof)
    V = 2 * math.pi * area_int(pts, 1, 0)
    zc = 2 * math.pi * area_int(pts, 1, 1) / V
    Ixx = math.pi * area_int(pts, 3, 0)       # int x^2 dV = int y^2 dV
    Izz = 2 * math.pi * area_int(pts, 1, 2)   # int z^2 dV
    C = [[Ixx, 0, 0], [0, Ixx, 0], [0, 0, Izz - V * zc * zc]]
    return V, [0.0, 0.0, zc], C


def rot(axis, deg):
    ax = [c / math.sqrt(sum(a * a for a in axis)) for c in axis]
    t = math.radians(deg)
    c, s, C = math.cos(t), math.sin(t), 1 - math.cos(t)
    x, y, z = ax
    return [[c + x * x * C, x * y * C - z * s, x * z * C + y * s],
            [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
            [z * x * C - y * s, z * y * C + x * s, c + z * z * C]]


# ---------------------------------------------------------------- STEP writer
class W:
    def __init__(self):
        self.lines, self.n = [], 0

    def add(self, txt):
        self.n += 1
        self.lines.append(f'#{self.n} = {txt};')
        return self.n

    def pt(self, p): return self.add(f"CARTESIAN_POINT('',({p[0]:.15g},{p[1]:.15g},{p[2]:.15g}))")
    def dr(self, d): return self.add(f"DIRECTION('',({d[0]:.15g},{d[1]:.15g},{d[2]:.15g}))")

    def ax(self, o, z, x):
        return self.add(f"AXIS2_PLACEMENT_3D('',#{self.pt(o)},#{self.dr(z)},#{self.dr(x)})")


def build_solid(w, prof):
    """returns MANIFOLD_SOLID_BREP id"""
    verts = {}

    def vertex(p):
        key = tuple(round(c, 12) for c in p)
        if key not in verts:
            verts[key] = w.add(f"VERTEX_POINT('',#{w.pt(p)})")
        return verts[key]

    def P(rz, u): return (rz[0] * math.cos(u), rz[0] * math.sin(u), rz[1])

    # parallels: dict (r, z) -> [edge u 0..pi, edge pi..2pi] (both with same_sense .T.)
    par = {}

    def parallel(rz):
        key = (round(rz[0], 12), round(rz[1], 12))
        if key in par or rz[0] < 1e-12:
            return par.get(key)
        c = w.add(f"CIRCLE('',#{w.ax((0, 0, rz[1]), (0, 0, 1), (1, 0, 0))},{rz[0]:.15g})")
        e0 = w.add(f"EDGE_CURVE('',#{vertex(P(rz, 0))},#{vertex(P(rz, math.pi))},#{c},.T.)")
        e1 = w.add(f"EDGE_CURVE('',#{vertex(P(rz, math.pi))},#{vertex(P(rz, 0))},#{c},.T.)")
        par[key] = [e0, e1]
        return par[key]

    def meridian(s, u):
        a, b = P(s[1], u), P(s[2], u)
        if s[0] == 'L':
            d = [b[i] - a[i] for i in range(3)]
            L = math.sqrt(sum(x * x for x in d))
            line = w.add(f"LINE('',#{w.pt(a)},#{w.add(f'VECTOR({chr(39)}{chr(39)},#{w.dr([x / L for x in d])},{L:.15g})')})")
            return w.add(f"EDGE_CURVE('',#{vertex(a)},#{vertex(b)},#{line},.T.)")
        pts, R, a1, dd = arc_pts(s, 2)
        cr, cz = s[3]
        ctr = (cr * math.cos(u), cr * math.sin(u), cz)
        n = (-math.sin(u), math.cos(u), 0)            # normal of the meridian plane
        X = (math.cos(u), math.sin(u), 0)
        circ = w.add(f"CIRCLE('',#{w.ax(ctr, n, X)},{R:.15g})")
        # curve parameter t: point = ctr + R (cos t X + sin t (n x X)); find the sense of the arc
        Y = (n[1] * X[2] - n[2] * X[1], n[2] * X[0] - n[0] * X[2], n[0] * X[1] - n[1] * X[0])

        def tpar(q):
            d = [q[i] - ctr[i] for i in range(3)]
            return math.atan2(sum(d[i] * Y[i] for i in range(3)), sum(d[i] * X[i] for i in range(3)))
        mid = P(pts[1], u)
        t1, t2, tm = tpar(a), tpar(b), tpar(mid)
        fwd = ((tm - t1) % (2 * math.pi)) < ((t2 - t1) % (2 * math.pi))
        return w.add(f"EDGE_CURVE('',#{vertex(a)},#{vertex(b)},#{circ},{'.T.' if fwd else '.F.'})")

    faces = []
    for s in prof:
        (r1, z1), (r2, z2) = s[1], s[2]
        if r1 < 1e-12 and r2 < 1e-12:
            continue  # on the axis
        # surface and its natural normal in the (r, z) half-plane at the segment middle
        if s[0] == 'L':
            rm, zm = (r1 + r2) / 2, (z1 + z2) / 2
            if abs(z2 - z1) < 1e-12:
                surf = w.add(f"PLANE('',#{w.ax((0, 0, z1), (0, 0, 1), (1, 0, 0))})")
                nat = (0, 1)
            elif abs(r2 - r1) < 1e-12:
                surf = w.add(f"CYLINDRICAL_SURFACE('',#{w.ax((0, 0, 0), (0, 0, 1), (1, 0, 0))},{r1:.15g})")
                nat = (1, 0)
            else:
                up = (r2 - r1) * (z2 - z1) > 0          # radius grows with z
                Z = (0, 0, 1) if up else (0, 0, -1)
                ta = abs((r2 - r1) / (z2 - z1))
                surf = w.add(f"CONICAL_SURFACE('',#{w.ax((0, 0, z1), Z, (1, 0, 0))},{r1:.15g},{math.atan(ta):.15g})")
                nat = (1, -ta * Z[2])
        else:
            pts, R, a1, dd = arc_pts(s, 2)
            rm, zm = pts[1]
            cr, cz = s[3]
            if abs(cr) < 1e-12:
                surf = w.add(f"SPHERICAL_SURFACE('',#{w.ax((0, 0, cz), (0, 0, 1), (1, 0, 0))},{R:.15g})")
            else:
                surf = w.add(f"TOROIDAL_SURFACE('',#{w.ax((0, 0, cz), (0, 0, 1), (1, 0, 0))},{cr:.15g},{R:.15g})")
            nat = (rm - cr, zm - cz)
        # outward normal of a CCW profile: right-hand side of the segment direction
        dr_, dz_ = (r2 - r1, z2 - z1) if s[0] == 'L' else (-(zm - s[3][1]) * (1 if s[4] else -1), (rm - s[3][0]) * (1 if s[4] else -1))
        out = (dz_, -dr_)
        same = nat[0] * out[0] + nat[1] * out[1] > 0
        for h in (0, 1):
            ua, ub = h * math.pi, (h + 1) * math.pi
            m_a, m_b = meridian(s, ua % (2 * math.pi)), meridian(s, ub % (2 * math.pi))
            pa, pb = parallel(s[1]), parallel(s[2])
            oes = []
            if pa: oes.append((pa[h], True))        # along u at P1 from ua to ub
            oes.append((m_b, True))                 # meridian at ub from P1 to P2
            if pb: oes.append((pb[h], False))       # at P2 from ub back to ua
            oes.append((m_a, False))                # meridian at ua from P2 to P1
            # orientation of this loop against the outward normal (vector area of a sampled loop)
            um = (ua + ub) / 2
            nout = (out[0] * math.cos(um), out[0] * math.sin(um), out[1])
            loop = []
            for k in range(200):  # P1 parallel, meridian, P2 parallel, meridian back
                t = k / 200
                loop.append(P(s[1], ua + (ub - ua) * t))
            mp = polygon([s], 50)
            for q in (mp[:51] if s[0] == 'A' else [(r1 + (r2 - r1) * k / 50, z1 + (z2 - z1) * k / 50) for k in range(51)]):
                loop.append(P(q, ub))
            for k in range(200):
                loop.append(P(s[2], ub - (ub - ua) * k / 200))
            for q in reversed((mp[:51] if s[0] == 'A' else [(r1 + (r2 - r1) * k / 50, z1 + (z2 - z1) * k / 50) for k in range(51)])):
                loop.append(P(q, ua))
            A = [0.0, 0.0, 0.0]
            for i in range(len(loop)):
                p, q = loop[i], loop[(i + 1) % len(loop)]
                A[0] += 0.5 * (p[1] * q[2] - p[2] * q[1]); A[1] += 0.5 * (p[2] * q[0] - p[0] * q[2]); A[2] += 0.5 * (p[0] * q[1] - p[1] * q[0])
            ok = sum(A[i] * nout[i] for i in range(3)) > 0
            ids = [w.add(f"ORIENTED_EDGE('',*,*,#{e},{'.T.' if f else '.F.'})") for e, f in oes]
            el = w.add(f"EDGE_LOOP('',({','.join('#%d' % i for i in ids)}))")
            fb = w.add(f"FACE_OUTER_BOUND('',#{el},{'.T.' if ok else '.F.'})")
            faces.append(w.add(f"ADVANCED_FACE('',(#{fb}),#{surf},{'.T.' if same else '.F.'})"))
    shell = w.add(f"CLOSED_SHELL('',({','.join('#%d' % f for f in faces)}))")
    return w.add(f"MANIFOLD_SOLID_BREP('',#{shell})")


def write_assembly(path, names):
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
    top_pd, top_pds = product('REVOLVED_TEST')
    top_ax = w.ax((0, 0, 0), (0, 0, 1), (1, 0, 0))
    top_sr = w.add(f"SHAPE_REPRESENTATION('',(#{top_ax}),#{rc})")
    w.add(f"SHAPE_DEFINITION_REPRESENTATION(#{top_pds},#{top_sr})")
    for i, nm in enumerate(names):
        pd, pds = product(nm)
        own = w.ax((0, 0, 0), (0, 0, 1), (1, 0, 0))
        sr = w.add(f"SHAPE_REPRESENTATION('',(#{own}),#{rc})")
        w.add(f"SHAPE_DEFINITION_REPRESENTATION(#{pds},#{sr})")
        solid = build_solid(w, PROFILES[nm])
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
    with open(path, 'w') as f:
        f.write("ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('cylload test solids of revolution'),'2;1');\n"
                "FILE_NAME('revolved_solids.step','2026-11-01',(''),(''),'gen_revolved.py','','');\n"
                "FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));\nENDSEC;\nDATA;\n")
        f.write('\n'.join(w.lines))
        f.write("\nENDSEC;\nEND-ISO-10303-21;\n")


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else '.'
    names = list(PROFILES)
    write_assembly(os.path.join(out, 'revolved_solids.step'), names)
    with open(os.path.join(out, 'revolved_reference.csv'), 'w') as f:
        f.write('name,volume,xc,yc,zc,cxx,cxy,cxz,cyy,cyz,czz\n')
        for nm in names:
            V, c, C = props(PROFILES[nm])
            axis, ang, t = PLACEMENT[nm]
            R = rot(axis, ang)
            cg = [sum(R[i][k] * c[k] for k in range(3)) + t[i] for i in range(3)]
            Cg = [[sum(R[i][k] * C[k][l] * R[j][l] for k in range(3) for l in range(3)) for j in range(3)] for i in range(3)]
            f.write(f"{nm},{V:.12g},{cg[0]:.12g},{cg[1]:.12g},{cg[2]:.12g},{Cg[0][0]:.12g},{Cg[0][1]:.12g},{Cg[0][2]:.12g},"
                    f"{Cg[1][1]:.12g},{Cg[1][2]:.12g},{Cg[2][2]:.12g}\n")


if __name__ == '__main__':
    main()
