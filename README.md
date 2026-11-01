# CylLoad — hydraulic cylinder loads from the solid model of a machine

CylLoad computes the forces, pressures and dynamic loads of the hydraulic cylinders of a
planar mechanism (booms, parallelogram linkages, tilting towers, sliding masts) from the
data that a designer already has: joint coordinates, cylinder drawings, title-block
masses and the STEP file of the assembly. It is a small C++17 program and library with
no external dependencies.

What it does:

* **Mass properties from the solid model.** It reads the STEP (AP203/AP214) assembly
  tree with its placements and integrates the volume, centroid and inertia tensor of
  every B-rep solid with the divergence theorem. Planes, cylinders, cones, spheres and
  tori are integrated in their parameter domain (Green's theorem, Gauss–Legendre
  quadrature along the edges); B-spline and other faces are approximated by a fan of
  triangles on their boundary. Each part's result is checked with independent integration fields.
  Assembly properties are scaled to the title-block masses (`cylload masstable`).
* **Kinematics of an arbitrary planar mechanism.** Rigid bodies are connected by
  revolute and prismatic joints and driven by cylinders (or other two-point actuators).
  Positions are solved by Newton's method with continuation, so the assembly branch is
  kept.
* **Static cylinder loads by three independent routes.** It uses the principle of
  virtual work with exact sensitivities, F = ∂V/∂L. It also uses finite differences of
  the potential energy and the force and moment balance of every link, which also gives
  all joint reactions. Agreement of the three routes and the energy balance
  ∫F dL = ΔV are reported for every run.
* **Pressures, buckling and hydraulic dynamics.** It reports the pressure on the loaded
  side and the Euler buckling factor. It computes the reduced mass seen by each cylinder
  (Lagrange), the hydraulic stiffness and natural frequency, and the excess force for a
  linear valve-opening ramp. From these it derives the shortest ramp time that keeps the
  pressure below a limit over the whole stroke.
* **Time-domain verification.** A start-up simulation solves the Lagrange equation with
  the configuration-dependent reduced mass and oil compressibility, with the valve
  opening along a linear ramp. Maps over two actuators give pressures and local ramp
  times for a controller.
* **Uncertainty and sensitivity.** A Monte Carlo study runs over mass scales, masses,
  centroid shifts and positions of other actuators. It reports percentiles, Spearman
  rank correlations, first-order variance-based indices, standardised regression
  coefficients and R². The random stream is reproducible across compilers.
* **Publication figures.** SVG plots and mechanism drawings with thin lines, axes that
  end on a labelled tick and panel letters below the panels. Figures are black by
  default or use the colour-blind-safe Okabe–Ito palette, with an optional grid and
  shaded uncertainty bands. A mechanism can be drawn from its STEP model, with every
  part moved into the computed pose, item numbers and joint labels, or from line work
  traced from its 2D drawings (`--sketch`).

## Build

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

Any C++17 compiler works. Continuous integration builds and tests the code on Linux,
Windows (MSVC) and macOS, and runs a separate build with the address and
undefined-behaviour sanitizers. API documentation: `doxygen Doxyfile`. On Windows, CMake generates a Visual Studio
solution.

## Quick start

```
cylload info  examples/horty/tower.mech
cylload pose  examples/horty/tower.mech S8=973 S10=758
cylload sweep examples/horty/tower.mech S10 S8=973 --n 177 --out tilt.csv
cylload ramp  examples/horty/tower.mech S10 S8=973 --flow 56.1 --plim 140
cylload mc    examples/horty/tower.mech examples/horty/tower_mc.txt --sensitivity sens.csv --out mc.csv
cylload simulate examples/horty/tower.mech S10 S8=973 S10=978 --flow 11.25 --ramp 0.1 --out sim.csv
cylload map   examples/horty/tower.mech S10 winch S8=973 --yfrom 4942.5 --yto 10721.5 --flow 56.1 --out map.csv
cylload plot  examples/horty/figures/fig_static.txt --out fig.svg
cylload draw  examples/horty/tower.mech S8=973 S10=1183 --out mechanism.svg
cylload step  assembly.step --level 2 --out parts.csv
cylload masstable assembly.step examples/horty/tower_assign.csv --remainder 222.7 --out masses.csv
```

All lengths are in mm, masses in kg, forces in N and pressures in bar. Actuators that are
not given on the command line keep their reference length. `sweep` writes one row per
position: body rotations, the force from each of the three routes, pressure, reduced
mass, hydraulic stiffness, natural frequency and potential energy. The summary on the
standard error stream gives the extremes, the agreement of the routes, the energy
balance and the buckling factor. The input formats are described in
[docs/input_format.md](docs/input_format.md).

## Examples

* `examples/boom` is a single boom lifted by one cylinder. The test suite compares it
  with the closed-form solution from the law of cosines; it agrees to 1e-12.
* `examples/horty` is the positioning mechanism of the Horty AutoHammer orchard
  pile driver. It comprises:
  * a lifting parallelogram, a tilting tower and a mast sliding along the tower;
  * the hammer position along the mast;
  * the post-lifting beam on the carriage.

  It reproduces the published analysis of this machine (Tucki and Nowakowski,
  *Applied Sciences*, applsci-4655781). There are two mass-property variants:
  * centres of the bounding outlines with prism inertia (`tower_outline.mech`);
  * properties integrated from the solid model (`tower.mech`, adopted in the paper).

  The simplified input data are released with the permission of Arno Group. The
  design documentation and the STEP files are not included.
* `examples/usarya` contains four mechanisms of a trailed stone picker: folding of a working
  module, the drawbar linkage, hopper lifting and hopper tipping. The kinematics is
  compared with the closed-form triangle model of the machine (Mathcad) and the pressures
  with published values (see [examples/usarya/README.md](examples/usarya/README.md)).

`run_examples.sh` (or `run_examples.bat` on Windows) writes all results of the examples
to `examples/horty/results` and `examples/usarya/results`.

## Tests

`cylload_tests` contains about 3350 checks:
* the time-domain simulation against the ramp formula and the published simulation;
* sensitivity indices against an analytic linear model;
* closed-form solutions: boom, inclined slide, hydraulic formulas;
* STEP solids of revolution covering every analytic face type, a cone apex, sphere
  poles and a negative-major-radius torus, against Pappus integrals
  (`tests/data/gen_revolved.py`);
* the hand-written assembly tree;
* the Monte Carlo helpers;
* regressions against the published results of the Horty example;
* the stone picker example against its closed-form kinematic model and published pressures.

## Licence and citation

MIT licence, see [LICENSE.txt](LICENSE.txt). If you use CylLoad, please cite the SoftwareX paper
(see [CITATION.cff](CITATION.cff)).
