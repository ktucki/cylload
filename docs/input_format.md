# CylLoad input formats

All files are plain text. `#` starts a comment. Lengths are in mm, masses in kg, moments
of inertia in kg·m² and gravity in m/s².

## Mechanism file (`.mech`)

The model is described in its **reference configuration**, which is the position of the
solid model. Every body pose is measured from this configuration. All coordinates are
global coordinates of the reference configuration, so they can be copied directly from
the CAD model.

```
title  TEXT
gravity GX GY                         # default 0 -9.81
body NAME [NAME ...]                  # moving bodies; "ground" exists implicitly
revolute  NAME BODY_A BODY_B X Y      # pin joint at (X, Y)
prismatic NAME BODY_A BODY_B X Y UX UY  # BODY_B slides along (UX, UY) through (X, Y), no relative rotation
actuator NAME BODY_A XA YA BODY_B XB YB [bore=D rod=d min=Lmin max=Lmax mass=m]
mass NAME BODY X Y M [J] [group=G]    # lumped mass with its centre of gravity
force NAME BODY X Y FX FY             # constant external force [N] at (X, Y), global directions
masses FILE.csv                       # mass table (relative to the .mech file)
hydraulics [beta=Pa] [dead_volume=L] [dead_volume_piston=L] [dead_volume_rod=L]
           [back_pressure=BAR] [reduced_mass_pins=yes|no]
                                      # defaults 1e9 Pa, 0.3 L per chamber, 0 bar, yes
```

* An actuator is a two-point force element whose pin-to-pin length is a driving
  coordinate. With `bore`, `rod`, `min` and `max` it is a hydraulic cylinder. Half of its
  `mass` is lumped at each pin. Without these keys it is a generic actuator (winch, screw
  jack), and only its force is computed.
* A `force` is a constant load such as a hitch load or a digging force. It has the
  potential -f·r, so it enters the virtual work, the energy differences and the link
  equilibrium exactly like a weight (`examples/usarya/drawbar_load.mech`).
* The number of actuators must equal the number of degrees of freedom. The program
  checks this and explains the mismatch.
* Forces are positive when the actuator pushes (extends against the load). Pressures
  refer to the piston side for pushing and to the annulus for pulling. With
  `back_pressure` p_b in the other chamber (counterbalance or meter-out valve), the loaded
  side needs p = F/A_p + p_b A_r/A_p when pushing and p = −F/A_r + p_b A_p/A_r when pulling.
* `beta` is the effective bulk modulus of oil and hoses; `dead_volume` (both chambers),
  `dead_volume_piston` and `dead_volume_rod` are the oil volumes between the piston and the
  valve at the end of the stroke (hoses, pipes, valve passages).
* The cylinder mass is lumped half at each pin, also in the reduced mass;
  `reduced_mass_pins=no` leaves these lumped masses out of the reduced mass (not out of the
  statics), as in models that treat the cylinder as massless.
* The Euler load in `sweep` is that of the piston rod alone, pinned at both ends, with the
  free length equal to the pin-to-pin distance at full extension. This overestimates the
  rod length and ignores the stiffer barrel; it is a screening value, not a replacement
  for a stepped-column check (e.g. ISO/TS 13725) with the guide length of the cylinder.

### Grammar

One statement per line; words are separated by spaces; keys of the form `key=value`
may come in any order after the positional fields. Names are case-sensitive and
must not contain spaces. Bodies must be declared with `body` before they are used, and
joints, actuators and masses may then come in any order. The file is checked on reading:
unknown keywords and keys, undeclared bodies, duplicate names, a cylinder without `bore`,
`rod`, `min` or `max`, a reference length outside the stroke, and a number of actuators
different from the number of degrees of freedom are reported with the line number.

## Mass table (CSV)

```
name,body,group,mass_kg,x_mm,y_mm,J_kgm2,source
upper bracket,AC,AC,86.81,-300.0,592.9,15.12,STEP 366-1.2.0
```

`group` (default: the body name) is used to select items in Monte Carlo studies. `J_kgm2`
is the moment of inertia about the centre of gravity and the axis normal to the plane.
It is used only for the reduced mass. Columns may come in any order, and `;` is
accepted as a separator.

## Assignment file for `cylload masstable`

```
name,body,group,mass_kg,prefix,remainder,x_mm,y_mm,J_kgm2
mast 18,M,M,178.0,366-1.8.0,yes,,,
post,M,post,56.0,,no,-1250,990,74.67
cylinder 8,-,S8,34.0,,yes,,,
```

* For a row with a `prefix`, all assembly-tree groups whose name starts with the prefix
  are combined:
  * groups are the path element at `--level`, where 1 is the children of the root;
  * the centroid and the squared radius of gyration are taken from the solid model
    (`--properties solid`) or from the bounding outline and a uniform prism
    (`--properties outline`);
  * the mass is the title-block value `mass_kg`, so the density is uniform within the
    assembly.
* A row without a prefix takes its `x_mm`, `y_mm` and `J_kgm2` as given.
* A row with a prefix may also give `x_mm`, `y_mm` (and `J_kgm2`); these replace the values
  from the solid model while the mass stays the title-block value. Use it for bought-in
  components (motors, gearboxes, valve blocks, hammers) whose supplier models are hollow
  shells or envelopes, with the centre of gravity from the data sheet
  (`examples/excavator/excavator_assign.csv`, valve block).
* `--remainder KG` distributes an unidentified mass in proportion to the masses of the rows
  with `remainder=yes`: each of these rows is multiplied by the same factor
  1 + remainder / (sum of their masses), so the remainder is placed at their centres of
  gravity, in proportion to their masses. Rows with `body` = `-` take part in this distribution but are not
  written; use them for cylinders, whose mass is set in the `.mech` file. Their scaled
  masses are printed.
* `--plane` selects the machine plane in STEP coordinates (default `xy`).
* `--bspline parametric|fan` (also for `step`): B-spline faces integrated in their parameter
  domain (default) or replaced by triangle fans as in version 1.1, e.g. to measure the effect.
* Lengths in the STEP file are converted from the unit of the file (mm, cm, m, inch, foot)
  to mm; the unit is printed. A file without a length unit is read as mm.

## Monte Carlo specification

```
seed 2026
samples 3000
bootstrap 1000                    # resamples for the confidence intervals (0 = none)
scale  LABEL SELECTORS LO HI      # multiply masses (and cylinder masses) by U(LO, HI)
mass   LABEL SELECTORS LO HI      # set masses to U(LO, HI)
shift  LABEL SELECTORS DXLO DXHI DYLO DYHI   # move centres of gravity (two inputs: LABEL_dx, LABEL_dy)
length LABEL ACTUATOR LO HI       # length of another actuator, e.g. a hammer position
output LABEL ACTUATOR [ACT=L ...] # force and pressure of ACTUATOR in the given pose
```

Selectors are `all`, `group:G`, `item:NAME` and `body:B`, separated by commas.
`group:NAME` and `item:NAME` with an actuator name select the mass of that cylinder.
Inputs are independent and uniformly distributed (bounds without a known distribution
within them). Spearman rank correlations are given between every input and the
pressure of every output (the force for an actuator without cylinder data), together with
the first-order indices S_i (20 equal-probability bins), standardised regression
coefficients and R². The binned estimator of S_i has a positive bias of about
(bins − 1)/N per input, so with many inputs and a small N the sum of S_i can exceed 1.
Bootstrap intervals (95 %) are given for P5, P95 and every S_i; the S_i interval of an
inactive input reaches zero.

## Outputs

* `pose` prints the actuator forces from the three routes, the pressures, the joint
  reactions and the body rotations.
* `sweep` sweeps the length of one actuator. `--report ACT` reports the loads of
  another actuator during the sweep, e.g. the tilting cylinder against the hammer
  position. The sweep writes CSV with columns `L_mm`, `rot_<body>_deg`, `F_N`, `F_fd_N`, `F_eq_N`,
  `p_bar`, `m_red_kg`, `k_hyd_N_per_m`, `f_n_Hz` and `V_J`. A summary goes to the
  standard error stream.
* `ramp` gives the minimum valve-ramp times for extension and retraction and a table of
  peak pressures against the ramp time.
* `mc` writes CSV with all samples. Percentiles with their intervals, rank correlations
  and sensitivity indices go to the standard error stream; `--sensitivity FILE` writes them
  as CSV.
* `step` and `masstable` write CSV with bounding boxes, volumes, centroids and squared
  radii of gyration, and the consistency checks `volume_check` and `moment_check`.

## Plot specification (`cylload plot`)

```
size W_MM H_MM                 # whole figure
columns N                      # panels per row
font_pt 8                      # optional; font "Arial, Helvetica, sans-serif"
palette okabe-ito              # optional: okabe-ito | black | list of colours (#RRGGBB ...)
grid                           # optional light grid at the major ticks
panel (a)                      # starts a panel; the letter is printed below it
title "..."                    # optional title above the plot area
xlabel "Hammer position s_{b} (m)"   # _{..} subscript, ^{..} superscript
ylabel "..."
xrange MIN MAX [STEP]          # otherwise automatic: the axis ends on a labelled tick
yrange MIN MAX [STEP]
xlog / ylog
legend below|inside|none
legend_columns N               # entries per legend row below the plot (default: as many as fit)
legend_rows N                  # rows reserved for the legend, so that panels side by side keep the same plot height
series FILE.csv XCOL YCOL [label=".."] [line=solid|dashed|dotted|dashdot|none]
       [marker=circle|circle-open|square|square-open|triangle|triangle-open|diamond|diamond-open|cross]
       [every=N] [width=0.2] [xscale=1] [yscale=1] [xoffset=0] [yoffset=0] [where=COL:VALUE]
point X Y [label=".."] [marker=circle-open]  # a single published value
hline VALUE [label=".."] [line=dashed]
vline VALUE [label=".."] [line=dashed]
band FILE.csv XCOL LOCOL HICOL [label=".."] [color=#0072B2] [opacity=0.18]
     [xscale=1] [xoffset=0] [yscale=1]           # shaded area between two columns
```

`series`, `point`, `hline` and `vline` also accept `color=#RRGGBB`. Figures have thin
lines (0.2 mm by default). Without `palette` they are black on white; with a palette,
series without an explicit colour take its colours in order. Line styles and markers are
kept in both cases, so the figures can also be printed in black and white. A `#` starts
a comment only as a separate word. CSV paths are relative to the specification file or to
`--data DIR`. `where=` selects the rows of a map, e.g. one hammer position.

## Mechanism drawing (`cylload draw`)

```
cylload draw MODEL [ACT=L ...] [--overlay ACT=L,ACT=L ...] [--color yes] [--labels no] [--items no] [--cg yes]
             [--title TEXT] [--font_pt 7] [--view XMIN,YMIN,XMAX,YMAX]
             [--width MM --height MM | --scale MM_PER_M [--ybottom MM] [--height MM]]
```

Each link is drawn as the rounded convex hull of its joints, cylinder pins and `outline`
points (its centres of mass when it has fewer than two such points); larger links are
drawn first. Two drawing-only keywords of the mechanism file help here:

```
outline BODY X1 Y1 X2 Y2 ...   # extra outline points, e.g. the corners of a mast
label   BODY TEXT              # item label of the link, e.g. an element number
```

Joint and cylinder labels are placed next to their marks where they cover the fewest
other marks. Link labels are item numbers as in an assembly drawing: a leader with a dot
inside the link leads to a column left or right of the drawing, ordered by height so
that leaders do not cross (`--items no` leaves them out). A 1 m scale bar is added; `--cg yes` also marks the centre of gravity of every link. With
`--scale` the drawing size follows from the mechanism; give the same `--scale`,
`--ybottom` (model ordinate at the lower margin) and `--height` to several drawings to
show them in one scale with the ground at the same level. `--view` draws only a region
of the model (for example the linkage of a long mast); links drawn from outlines get
break marks where they leave it. Overlay poses are dashed.

With `--step FILE --parts CSV` the links are drawn from the solid model instead of their
outlines. Every solid of the STEP assembly is projected on the machine plane (`--plane`,
default `xy`, seen from the side of the third axis) and moved rigidly with its link, or
with the barrel or rod of a cylinder; solids nearer to the viewer are drawn later. The
parts table assigns assemblies by the prefix of their name at `--level` (default 2):

```
prefix,body
366-1.1.1,ground                            # frame
366-1.2.0,AC                                # a link
366-1.9.0 Silownik hydraul 100x50x350 Rura,S8.barrel
366-1.9.0 Silownik hydraul 100x50x350 Tłok,S8.rod
366-1.28.0,-                                # left out
```

The first matching row wins; assemblies without a row are left out. Plane faces are
filled with their holes, and curved faces with the hull of their boundary, so that a
tube seen from the side is filled. Item numbers then point at the drawn parts.

A machine without a solid model can be drawn from its 2D documentation with `--sketch FILE`:
each line of the file is `BODY [fill] X1 Y1 X2 Y2 ...` (mm, reference configuration), where
BODY is `ground`, a link or `ACT.barrel` / `ACT.rod`; the polylines are moved with their link
like the solids of `--step`, coloured by link, and `fill` closes and fills a polygon (a
silhouette). Later lines are drawn on top; when the sketch draws the frame, the ground symbol
is left out. `examples/usarya/hopper_sketch.txt` holds the 1978 lines of the hopper system
traced from its assembly drawing:

```
cylload draw examples/usarya/hopper.mech lift=1890.388 tip=1588.435 --sketch examples/usarya/hopper_sketch.txt --scale 10 --color yes --items no --out hopper_tipped.svg
```

## Simulation and maps

* `simulate MODEL ACT [ACT=L ...] --flow LPM [--ramp S] [--dir ext|ret] [--start L]`
  `[--damping ZETA]` writes the time history (`t_s`, `L_mm`, `v_mm_s`, `p1_bar`,
  `p2_bar`, `dF_N`). With `--ramps T1,T2,...` it writes instead a table of simulated and
  formula excess forces and dynamic factors κ against the ramp time.
* `map MODEL ACT_X ACT_Y [ACT=L ...] [--report ACT] [--flow LPM] [--plim BAR]` writes
  one row per grid node with the rotations, force, pressure, reduced mass, natural
  frequency and local minimum ramp times (`inf` where the static load alone exceeds the
  limit).
