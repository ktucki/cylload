# Changelog

## 1.2.0 — 2026-11-10
* STEP: B-spline surface faces (also rational and periodic) are integrated in their
  parameter domain. The boundary is mapped to (u, v) by point inversion on the surface,
  and Green's theorem is applied in (u, v) with the exact surface normal. Faces whose
  inversion fails, or whose loop winds around a closed direction, still fall back to the
  triangle fan and are counted. `--bspline fan` restores the behaviour of 1.1.
* `step` and `masstable` report the number of faces integrated in (u, v), the number of
  triangulated faces and their share of the surface area (CSV columns
  `faces_bspline_parametric`, `faces_triangulated`, `area_mm2`, `area_triangulated_mm2`).
* STEP length units: the unit in force (`GLOBAL_UNIT_ASSIGNED_CONTEXT`) is read; millimetre,
  centimetre, metre, inch and foot are converted to mm, other units are rejected, and a file
  without a unit is read as mm with a note ("mm (assumed)").
* `masstable`: a row with a prefix may give `x_mm`, `y_mm` and `J_kgm2`, which then replace the
  values from the solid model, for bought-in parts exported as shells or envelopes.
* `hydraulics`: `dead_volume_piston`, `dead_volume_rod` (`dead_volume` sets both),
  `back_pressure=BAR` (pressure in the unloaded chamber, added to the loaded side through
  the area ratio) and `reduced_mass_pins=no` (cylinder pin masses left out of the reduced
  mass only).
* `mc`: 95 % bootstrap intervals (`bootstrap B` in the specification or `--bootstrap B`,
  default 1000; 0 switches them off) for P5 and P95 (percentile intervals) and for every
  S_i (basic intervals, which correct the binning bias). The sensitivity CSV has the new
  columns `first_order_ci_lo`, `first_order_ci_hi`, `p5_ci_lo`, `p5_ci_hi`, `p95_ci_lo`
  and `p95_ci_hi`.
* `simulate --ramps`: column `dF_sim_less_damping_drag_N`, the excess force over the
  quasi-steady force including the viscous drag at the commanded velocity.
* New open example `examples/excavator`: a synthetic backhoe attachment from the STEP file
  to forces, with B-spline faces, nested assemblies, a bought-in part with data-sheet
  centre of gravity, a four-bar bucket linkage, exact reference mass properties
  (`gen_excavator.py`), an independent check of forces, reduced masses and natural
  frequencies in absolute coordinates (`check_independent.py`) and the same comparison with
  Exudyn (`check_exudyn.py`).
* `examples/usarya/cases/drawbar_frame_h*.mech`: the drawbar linkage loaded by the frame weight
  at its centre of gravity, for five heights of the centre of gravity.
* Tests: B-spline solids against exact integrals (`tests/data/gen_bspline.py`), length units,
  the excavator example, bootstrap intervals, back-pressure, separate dead volumes and damping.

## 1.1.0 — 2026-11-01
* `simulate`: time-domain start-up of a cylinder (Lagrange equation with the
  configuration-dependent reduced mass, oil compressibility, ideal flow valve with a
  linear opening ramp, optional damping), integrated with RK4. It is compared with the
  ramp formula.
* `map`: grid over two actuators with the loads, natural frequency and local minimum
  ramp time of a reported cylinder, for ramp-time maps in a controller.
* `mc`: first-order variance-based indices (20 equal-probability bins), standardised
  regression coefficients and R²; `--sensitivity` writes them to a CSV file.
* `sweep --report`: the loads of another actuator during a sweep, e.g. tilting pressure
  against the hammer position.
* `plot` and `draw`: publication figures as SVG:
  * thin lines, black by default or the Okabe–Ito colour palette (`palette okabe-ito`),
    with line styles and markers kept, so the figures also read in black and white;
  * optional light grid, shaded bands (e.g. Monte Carlo P5–P95) and `color=` per item;
  * axes ending on a labelled tick and panel letters below the panels;
  * subscripts written as `_{..}` and superscripts as `^{..}`;
  * a mechanism drawing in one or more poses: links as rounded hulls of their joints,
    pins and `outline` points, item numbers (`label`) with leaders, labels placed clear of the
    drawing, a 1 m scale bar, `--scale` and `--ybottom` for several drawings in one
    scale with the ground aligned, `--view` with break marks, `--title`, `--cg yes`
    for centres of gravity and `--color yes`;
  * `draw --step --parts`: the links drawn from the solid model, every solid projected on
    the machine plane and moved with its link or cylinder barrel or rod, painted back to
    front (`examples/horty/tower_parts.csv`).
  * `draw --sketch`: 2D line work traced from drawings, each polyline assigned to a
    link or cylinder part and moved like the solids (`examples/usarya/hopper_sketch.txt`).
* `force NAME BODY X Y FX FY`: constant external forces (hitch load, digging force) in all
  three force routes and in the joint reactions; `examples/usarya/drawbar_load.mech`
  loads the drawbar by the hitch load as in the published analysis.
* `examples/usarya`: four mechanisms of a trailed stone picker checked against a
  closed-form kinematic model and published pressures (regression tests).
* `sweep`: a cylinder that is never pushed is reported as not loaded in compression
  instead of a buckling factor of 0.
* Continuous integration on Linux, Windows and macOS, and a sanitizer build.

## 1.0.0 — 2026-11-01
* First release: STEP assembly reader with B-rep mass properties, `masstable`,
  kinematics, the three routes to the static loads, hydraulic dynamics, valve ramps,
  Monte Carlo, the Horty example and the test suite.
