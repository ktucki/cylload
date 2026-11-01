# Changelog

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
