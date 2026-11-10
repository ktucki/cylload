#!/usr/bin/env python3
# SPDX-License-Identifier: MIT
"""Generate the open example: a synthetic backhoe attachment as a STEP AP203 assembly.

The model is invented for this example (it is not a real machine). It has the structure
of a real CAD export:

  EX-100 Backhoe attachment                       (root)
    EX-110 Boom            side plates with rounded ends (rational B-spline faces),
                           pin bosses, cylinder bracket casting (bicubic B-spline face)
    EX-120 Arm             side plates, bosses
    EX-130 Bucket          curved shell (rational B-spline faces), ear plates
    EX-140 Rocker, EX-150 Link        bucket linkage
    EX-160 Valve block     bought-in part, exported only as a simplified envelope
    EX-171..173, 181..183  barrels and piston rods of the boom, arm and bucket cylinders
    EX-105 Upper frame     ground

Machine plane x-y (x forward, y up), z across the machine, millimetres.
The exact mass properties of every level-2 assembly are computed here from the
geometric definitions (independently of CylLoad) and written to
excavator_reference.csv.
Usage: python3 gen_excavator.py [OUT_DIR]
"""
import math, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', '..', 'tests', 'data'))
import gen_revolved as gr          # noqa: E402  STEP writer, solids of revolution
import gen_bspline as gb           # noqa: E402  extruded sections with arcs, B-spline bump

S2 = math.sqrt(0.5)

# ------------------------------------------------------------------ joints, reference pose
J = dict(A=(0.0, 0.0),            # boom foot
         Bg=(450.0, -650.0),      # boom cylinder, ground end
         Bb=(1600.0, 450.0),      # boom cylinder, boom end
         Ab=(1800.0, 1650.0),     # arm cylinder, boom end
         F=(3300.0, 1500.0),      # arm pin (boom tip)
         At=(3150.0, 2050.0),     # arm cylinder, arm end (arm tail)
         Cb=(3420.0, 1250.0),     # bucket cylinder, arm end
         R=(3620.0, 0.0),         # rocker pivot on the arm
         K=(3880.0, 60.0),        # rocker - link - bucket cylinder rod
         G=(3700.0, -300.0),      # bucket pivot
         L=(3960.0, -230.0))      # link - bucket


def sub(a, b): return (a[0] - b[0], a[1] - b[1])
def add(a, b): return (a[0] + b[0], a[1] + b[1])
def mul(a, s): return (a[0] * s, a[1] * s)
def unit(a):
    n = math.hypot(*a)
    return (a[0] / n, a[1] / n)
def perp(a): return (-a[1], a[0])


# ------------------------------------------------------------------ section profiles in the machine plane
def slot_plate(pts, r):
    """r: half width, or a list of half widths at the points (tapered plate)."""
    rs = list(r) if isinstance(r, (list, tuple)) else [r] * len(pts)
    """Plate around a polyline of pin centres pts with half width r: straight sides,
    each end a half circle made of two 90 degree arcs (rational B-spline faces).
    Returns a CCW profile of ('L', p, q) and ('A', p, q, centre) segments."""
    left, right = [], []
    n = len(pts)
    for i in range(n):
        if i == 0: d = unit(sub(pts[1], pts[0]))
        elif i == n - 1: d = unit(sub(pts[-1], pts[-2]))
        else:  # mitre at an interior bend
            d1, d2 = unit(sub(pts[i], pts[i - 1])), unit(sub(pts[i + 1], pts[i]))
            d = unit(add(d1, d2))
        nrm = perp(d)
        if 0 < i < n - 1:
            d1 = unit(sub(pts[i], pts[i - 1]))
            s = rs[i] / (nrm[0] * perp(d1)[0] + nrm[1] * perp(d1)[1])
        else:
            s = rs[i]
        left.append(add(pts[i], mul(nrm, s)))
        right.append(add(pts[i], mul(nrm, -s)))
    prof = []
    # right side forward, end cap around pts[-1], left side backward, cap around pts[0]
    for i in range(n - 1):
        prof.append(('L', right[i], right[i + 1]))
    d = unit(sub(pts[-1], pts[-2]))
    tip = add(pts[-1], mul(d, rs[-1]))
    prof.append(('A', right[-1], tip, pts[-1]))
    prof.append(('A', tip, left[-1], pts[-1]))
    for i in range(n - 1, 0, -1):
        prof.append(('L', left[i], left[i - 1]))
    d = unit(sub(pts[0], pts[1]))
    tail = add(pts[0], mul(d, rs[0]))
    prof.append(('A', left[0], tail, pts[0]))
    prof.append(('A', tail, right[0], pts[0]))
    return prof


def box(x0, y0, x1, y1):
    return [('L', (x0, y0), (x1, y0)), ('L', (x1, y0), (x1, y1)), ('L', (x1, y1), (x0, y1)), ('L', (x0, y1), (x0, y0))]


def bucket_shell():
    """C-shaped bucket section around the pivot G: back wall, rounded bottom (two concentric
    90 degree arcs), floor plate ending in the cutting edge. Wall thickness 14 mm."""
    gx, gy = J['G']
    c = (gx + 250.0, gy - 350.0)    # centre of the bottom radius
    Ro, t = 420.0, 14.0
    Ri = Ro - t
    # outer: from the top of the back wall down, around the bottom, to the cutting edge
    a = (c[0] - Ro, gy + 60.0)
    b = (c[0] - Ro, c[1])
    d = (c[0], c[1] - Ro)
    e = (c[0] + 520.0, c[1] - Ro)
    e2 = (c[0] + 520.0, c[1] - Ri)
    d2 = (c[0], c[1] - Ri)
    b2 = (c[0] - Ri, c[1])
    a2 = (c[0] - Ri, gy + 60.0)
    return [('L', a2, a), ('L', a, b), ('A', b, d, c), ('L', d, e), ('L', e, e2), ('L', e2, d2),
            ('A', d2, b2, c), ('L', b2, a2)]


def bucket_side():
    """Side plate: the outer outline of the shell closed by a straight line from the
    cutting edge to the top of the back wall."""
    sh = bucket_shell()
    a, b, arc, d, e = sh[0][2], sh[1][2], sh[2], sh[3][1], sh[3][2]
    return [('L', a, b), arc, ('L', d, e), ('L', e, a)]


def orient_ccw(prof):
    """Reverse a profile that runs clockwise (shoelace on the segment end points and arc
    midpoints)."""
    pts = []
    for s in prof:
        pts.append(s[1])
        if s[0] == 'A':
            pts.append(gb.arc_pt(('A', s[1], s[2], s[3]), 0.5))
    A = sum(pts[i][0] * pts[(i + 1) % len(pts)][1] - pts[(i + 1) % len(pts)][0] * pts[i][1] for i in range(len(pts)))
    if A > 0:
        return prof
    out = []
    for s in reversed(prof):
        out.append(('L', s[2], s[1]) if s[0] == 'L' else ('A', s[2], s[1], s[3]))
    return out


# ------------------------------------------------------------------ parts
# Plates are sections in the machine plane extruded across the machine (z from z0 to z1).
# gen_bspline's extrusion runs along its local x with the section in local (y, z), so the
# placement maps local (x, y, z) -> global (z, x, y).
EXTRUDE_R = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]  # columns: local x->(0,0,1), y->(1,0,0), z->(0,1,0)

PARTS = {}   # name -> dict(kind, data, R, t)
TREE = []    # (assembly, [part names])


def plate(name, prof, z0, z1):
    gb.BARS[name] = dict(L=z1 - z0, prof=orient_ccw(prof))
    PARTS[name] = dict(kind='bar', R=EXTRUDE_R, t=(0.0, 0.0, z0))


def boss(name, centre, ro, ri, z0, z1):
    """Pin boss: tube along z (solid of revolution, local z = global z)."""
    gr.PROFILES[name] = [('L', (ri, 0), (ro, 0)), ('L', (ro, 0), (ro, z1 - z0)), ('L', (ro, z1 - z0), (ri, z1 - z0)),
                         ('L', (ri, z1 - z0), (ri, 0))]
    PARTS[name] = dict(kind='rev', R=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], t=(centre[0], centre[1], z0))


def rod(name, p, q, ro, ri, L0, L1):
    """Cylinder barrel (ri > 0) or piston rod (ri = 0) on the line p -> q, from L0 to L1 mm from p."""
    d = sub(q, p)
    if ri > 0:
        prof = [('L', (ri, 0), (ro, 0)), ('L', (ro, 0), (ro, L1 - L0)), ('L', (ro, L1 - L0), (ri, L1 - L0)), ('L', (ri, L1 - L0), (ri, 0))]
    else:
        prof = [('L', (0, 0), (ro, 0)), ('L', (ro, 0), (ro, L1 - L0)), ('L', (ro, L1 - L0), (0, L1 - L0)), ('L', (0, L1 - L0), (0, 0))]
    gr.PROFILES[name] = prof
    u = unit(d)
    # local z along the cylinder axis
    # local x in the machine plane, so the seams of the two half faces lie on the outline in side view
    Z = (u[0], u[1], 0.0)
    X = (-u[1], u[0], 0.0)
    Y = (Z[1] * X[2] - Z[2] * X[1], Z[2] * X[0] - Z[0] * X[2], Z[0] * X[1] - Z[1] * X[0])
    R = [[X[i], Y[i], Z[i]] for i in range(3)]
    o = add(p, mul(u, L0))
    PARTS[name] = dict(kind='rev', R=R, t=(o[0], o[1], 0.0))


def casting(name, origin):
    """Cylinder bracket casting: block with a doubly curved bicubic B-spline top
    (gen_bspline 'bump', enlarged to 240 x 180 mm), lying on the boom; local z = global y."""
    PARTS[name] = dict(kind='bump', R=[[1, 0, 0], [0, 0, 1], [0, -1, 0]], t=origin)


def build():
    W2, T = 160.0, 20.0  # half the inner width between side plates, plate thickness
    # boom: bent, tapered box, foot A -> knee -> tip F
    knee = (1900.0, 1250.0)
    for side, z0 in (('L', W2), ('R', -W2 - T)):
        plate(f'EX-111 boom side plate {side}', slot_plate([J['A'], (900.0, 610.0), knee, (2650.0, 1380.0), J['F']],
                                                           [150.0, 200.0, 250.0, 185.0, 140.0]), z0, z0 + T)
    boss('EX-112 boom foot boss', J['A'], 95.0, 50.0, -W2, W2)
    boss('EX-113 boom tip boss', J['F'], 90.0, 45.0, -W2, W2)
    boss('EX-114 boom cylinder boss', J['Bb'], 70.0, 40.0, -W2, W2)
    for side, z0 in (('L', 60.0), ('R', -80.0)):
        plate(f'EX-117 boom cylinder bracket {side}', slot_plate([(1560.0, 900.0), J['Bb']], [130.0, 80.0]), z0, z0 + 20.0)
    plate('EX-115 arm cylinder bracket', slot_plate([(1760.0, 1300.0), J['Ab']], [140.0, 80.0]), -60.0, 60.0)
    casting('EX-116 bracket casting', (2200.0, 1330.0, 90.0))
    TREE.append(('EX-110 Boom', [p for p in PARTS if p.startswith('EX-11')]))
    # arm: tail At -> F -> G
    n0 = len(PARTS)
    for side, z0 in (('L', 90.0), ('R', -110.0)):
        plate(f'EX-121 arm side plate {side}', slot_plate([J['At'], J['F'], (3500.0, 600.0), J['G']],
                                                          [95.0, 165.0, 140.0, 95.0]), z0, z0 + 20.0)
    boss('EX-122 arm pin boss', J['F'], 75.0, 45.0, -90.0, 90.0)
    boss('EX-123 bucket pivot boss', J['G'], 60.0, 35.0, -90.0, 90.0)
    boss('EX-124 rocker pivot boss', J['R'], 50.0, 30.0, -90.0, 90.0)
    boss('EX-125 bucket cylinder boss', J['Cb'], 55.0, 35.0, -90.0, 90.0)
    TREE.append(('EX-120 Arm', list(PARTS)[n0:]))
    # bucket
    n0 = len(PARTS)
    plate('EX-131 bucket shell', bucket_shell(), -400.0, 400.0)
    for side, z0 in (('L', 112.0), ('R', -132.0)):
        plate(f'EX-132 bucket ear {side}', slot_plate([J['G'], J['L']], 60.0), z0, z0 + 20.0)
    for side, z0 in (('L', 400.0), ('R', -416.0)):
        plate(f'EX-133 bucket side {side}', bucket_side(), z0, z0 + 16.0)
    sh = bucket_shell()
    ex, ey = sh[3][2]          # outer point of the cutting edge
    for k, z0 in enumerate((-330.0, -170.0, -10.0, 150.0, 310.0)):
        plate(f'EX-134 tooth {k + 1}', [('L', (ex - 60.0, ey), (ex + 150.0, ey + 6.0)), ('L', (ex + 150.0, ey + 6.0), (ex + 150.0, ey + 26.0)),
                                        ('L', (ex + 150.0, ey + 26.0), (ex - 60.0, ey + 60.0)), ('L', (ex - 60.0, ey + 60.0), (ex - 60.0, ey))],
              z0, z0 + 40.0)
    TREE.append(('EX-130 Bucket', list(PARTS)[n0:]))
    # linkage
    n0 = len(PARTS)
    for side, z0 in (('L', 92.0), ('R', -112.0)):
        plate(f'EX-141 rocker {side}', slot_plate([J['R'], J['K']], 45.0), z0, z0 + 20.0)
    TREE.append(('EX-140 Rocker', list(PARTS)[n0:]))
    n0 = len(PARTS)
    plate('EX-151 link', slot_plate([J['K'], J['L']], 55.0), -60.0, 60.0)
    TREE.append(('EX-150 Link', list(PARTS)[n0:]))
    # bought-in valve block on the boom: the supplier's STEP is a plain envelope
    n0 = len(PARTS)
    plate('EX-161 valve block envelope', box(1150.0, 700.0, 1450.0, 850.0), -90.0, 90.0)
    TREE.append(('EX-160 Valve block', list(PARTS)[n0:]))
    # cylinders: barrel from the base end, rod from the rod end (reference lengths)
    # barrel and rod are each 0.92 of the retracted length, so they overlap over the whole stroke
    for k, (nm, a, b, rb, rr, lmin) in enumerate((('Boom cylinder S1', 'Bg', 'Bb', 75.0, 45.0, 1300.0),
                                                  ('Arm cylinder S2', 'Ab', 'At', 70.0, 40.0, 1200.0),
                                                  ('Bucket cylinder S3', 'Cb', 'K', 60.0, 35.0, 1050.0))):
        n0 = len(PARTS)
        rod(f'EX-17{k + 1}.1 barrel', J[a], J[b], rb, 0.0, 60.0, 0.92 * lmin)
        TREE.append((f'EX-17{k + 1} {nm} barrel', list(PARTS)[n0:]))
        n0 = len(PARTS)
        rod(f'EX-18{k + 1}.1 piston rod', J[b], J[a], rr, 0.0, 40.0, 0.92 * lmin)
        TREE.append((f'EX-18{k + 1} {nm} piston rod', list(PARTS)[n0:]))
    n0 = len(PARTS)
    plate('EX-106 boom mounting bracket', [('L', (-260.0, -860.0), (620.0, -860.0)), ('L', (620.0, -860.0), (620.0, -560.0)),
                                            ('L', (620.0, -560.0), (180.0, -60.0)), ('L', (180.0, -60.0), (-180.0, -60.0)),
                                            ('L', (-180.0, -60.0), (-260.0, -860.0))], -60.0, 60.0)
    plate('EX-107 upper frame', box(-2200.0, -1060.0, 700.0, -860.0), -600.0, 600.0)
    boss('EX-108 boom cylinder pin boss', J['Bg'], 75.0, 40.0, -100.0, 100.0)
    TREE.append(('EX-105 Upper frame', list(PARTS)[n0:]))


# ------------------------------------------------------------------ exact properties
KEYS = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (2, 0, 0), (0, 2, 0), (0, 0, 2), (1, 1, 0), (0, 1, 1), (1, 0, 1)]


def local_props(name):
    k = PARTS[name]['kind']
    if k == 'bar':
        I = gb.bar_props(name)
        V = I[(0, 0, 0)]
        c = [I[(1, 0, 0)] / V, I[(0, 1, 0)] / V, I[(0, 0, 1)] / V]
        S = [[I[(2, 0, 0)], I[(1, 1, 0)], I[(1, 0, 1)]], [I[(1, 1, 0)], I[(0, 2, 0)], I[(0, 1, 1)]], [I[(1, 0, 1)], I[(0, 1, 1)], I[(0, 0, 2)]]]
        C = [[S[i][j] - V * c[i] * c[j] for j in range(3)] for i in range(3)]
        return V, c, C
    if k == 'bump':
        I = gb.bump_props()
        V = I[(0, 0, 0)]
        c = [I[(1, 0, 0)] / V, I[(0, 1, 0)] / V, I[(0, 0, 1)] / V]
        S = [[I[(2, 0, 0)], I[(1, 1, 0)], I[(1, 0, 1)]], [I[(1, 1, 0)], I[(0, 2, 0)], I[(0, 1, 1)]], [I[(1, 0, 1)], I[(0, 1, 1)], I[(0, 0, 2)]]]
        C = [[S[i][j] - V * c[i] * c[j] for j in range(3)] for i in range(3)]
        return V, c, C
    return gr.props(gr.PROFILES[name])


def global_props(name):
    V, c, C = local_props(name)
    R, t = PARTS[name]['R'], PARTS[name]['t']
    cg = [sum(R[i][k] * c[k] for k in range(3)) + t[i] for i in range(3)]
    Cg = [[sum(R[i][k] * C[k][l] * R[j][l] for k in range(3) for l in range(3)) for j in range(3)] for i in range(3)]
    return V, cg, Cg


def union(props):
    V = sum(p[0] for p in props)
    c = [sum(p[0] * p[1][i] for p in props) / V for i in range(3)]
    C = [[sum(p[2][i][j] + p[0] * (p[1][i] - c[i]) * (p[1][j] - c[j]) for p in props) for j in range(3)] for i in range(3)]
    return V, c, C


# ------------------------------------------------------------------ STEP output
def write_step(path):
    w = gr.W()
    ctx = w.add("APPLICATION_CONTEXT('configuration controlled 3D designs')")
    pc = w.add(f"PRODUCT_CONTEXT('',#{ctx},'mechanical')")
    pdc = w.add(f"PRODUCT_DEFINITION_CONTEXT('part definition',#{ctx},'design')")
    unit_ = w.add("(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))")
    rc = w.add(f"(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#{unit_}))REPRESENTATION_CONTEXT('',''))")
    counter = [0]

    def product(name):
        p = w.add(f"PRODUCT('{name}','{name}','',(#{pc}))")
        f = w.add(f"PRODUCT_DEFINITION_FORMATION_WITH_SPECIFIED_SOURCE('','',#{p},.NOT_KNOWN.)")
        d = w.add(f"PRODUCT_DEFINITION('design','',#{f},#{pdc})")
        s = w.add(f"PRODUCT_DEFINITION_SHAPE('','',#{d})")
        own = w.ax((0, 0, 0), (0, 0, 1), (1, 0, 0))
        sr = w.add(f"SHAPE_REPRESENTATION('',(#{own}),#{rc})")
        w.add(f"SHAPE_DEFINITION_REPRESENTATION(#{s},#{sr})")
        return d, sr, own

    def link(parent, child, R, t, label):
        (ppd, psr, _), (cpd, csr, cown) = parent, child
        place = w.ax(t, (R[0][2], R[1][2], R[2][2]), (R[0][0], R[1][0], R[2][0]))
        counter[0] += 1
        nauo = w.add(f"NEXT_ASSEMBLY_USAGE_OCCURRENCE('{counter[0]}','{label}','',#{ppd},#{cpd},$)")
        npds = w.add(f"PRODUCT_DEFINITION_SHAPE('','',#{nauo})")
        idt = w.add(f"ITEM_DEFINED_TRANSFORMATION('','',#{place},#{cown})")
        rr = w.add(f"(REPRESENTATION_RELATIONSHIP('','',#{csr},#{psr})REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#{idt})SHAPE_REPRESENTATION_RELATIONSHIP())")
        w.add(f"CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(#{rr},#{npds})")

    I3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    root = product('EX-100 Backhoe attachment')
    for asm, parts in TREE:
        a = product(asm)
        link(root, a, I3, (0.0, 0.0, 0.0), asm + ':1')
        for nm in parts:
            pp = product(nm)
            k = PARTS[nm]['kind']
            solid = gb.build_bar(w, nm) if k == 'bar' else gb.build_bump(w) if k == 'bump' else gr.build_solid(w, gr.PROFILES[nm])
            ab = w.add(f"ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{solid},#{pp[2]}),#{rc})")
            w.add(f"SHAPE_REPRESENTATION_RELATIONSHIP('','',#{pp[1]},#{ab})")
            link(a, pp, PARTS[nm]['R'], PARTS[nm]['t'], nm + ':1')
    with open(path, 'w') as f:
        f.write("ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('synthetic backhoe attachment, open CylLoad example'),'2;1');\n"
                "FILE_NAME('excavator.step','2026-11-10',('Karol Tucki'),(''),'gen_excavator.py','','');\n"
                "FILE_SCHEMA(('CONFIG_CONTROL_DESIGN'));\nENDSEC;\nDATA;\n")
        f.write('\n'.join(w.lines))
        f.write("\nENDSEC;\nEND-ISO-10303-21;\n")


def main():
    out = sys.argv[1] if len(sys.argv) > 1 else HERE
    # cylinder bracket casting: the B-spline block of gen_bspline enlarged by 3
    gb.BUMP.update(a=240.0, b=180.0, z=[[3 * v for v in row] for row in gb.BUMP['z']])
    build()
    write_step(os.path.join(out, 'excavator.step'))
    with open(os.path.join(out, 'excavator_reference.csv'), 'w') as f:
        f.write('# exact mass properties of the level-1 assemblies (mm, mm3, mm5), computed by gen_excavator.py\n')
        f.write('assembly,parts,volume,xc,yc,zc,cxx,cxy,cxz,cyy,cyz,czz\n')
        for asm, parts in TREE:
            V, c, C = union([global_props(p) for p in parts])
            f.write(f"{asm},{len(parts)},{V:.12g},{c[0]:.12g},{c[1]:.12g},{c[2]:.12g},{C[0][0]:.12g},{C[0][1]:.12g},"
                    f"{C[0][2]:.12g},{C[1][1]:.12g},{C[1][2]:.12g},{C[2][2]:.12g}\n")


if __name__ == '__main__':
    main()
