# Open example: synthetic backhoe attachment, from the STEP file to cylinder loads

This example runs the whole CylLoad workflow on data that anyone can inspect and
regenerate. The machine is invented for this purpose: a backhoe attachment with a bent box
boom, an arm, a bucket on a four-bar linkage and three cylinders. Its solid model has the
features of a real CAD export:

* an assembly tree with 13 assemblies and 32 parts, each with its own placement;
* 57 faces that are rational or non-rational B-spline surfaces: the rounded ends of the
  plates, the bucket bottom and a cast bracket with a doubly curved bicubic top;
* a bought-in valve block that the supplier exports only as an envelope, so its centre of
  gravity is taken from the (invented) data sheet;
* title-block masses that differ from the solid-model volumes, an unidentified remainder of
  60 kg and an operating load (bucket fill) without geometry.

| File | Content |
|:--|:--|
| `gen_excavator.py` | generates `excavator.step` and the exact mass properties `excavator_reference.csv` |
| `excavator.step` | STEP AP203 assembly, mm (machine plane x–y, z across the machine) |
| `excavator_assign.csv` | title-block masses of the assemblies, links, remainder flags, data-sheet CoG of the valve block |
| `excavator_masses.csv` | mass table written by `cylload masstable` |
| `excavator.mech` | mechanism: 5 links, 6 revolute joints, cylinders S1 (boom), S2 (arm), S3 (bucket) |
| `excavator_mc.txt` | Monte Carlo specification (masses ±10 %, bucket fill 0–810 kg, remainder 0–60 kg) |
| `excavator_parts.csv` | assignment of the assemblies to links and cylinder barrels and rods, for drawing |
| `run.sh` | all commands below |
| `check_independent.py` | independent check of mass properties, forces, reduced masses and natural frequencies |
| `check_exudyn.py` | the same mechanism in Exudyn (forces, reduced masses, natural frequencies) |

## Commands

From the repository root, after building:

```
sh examples/excavator/run.sh                        # writes examples/excavator/results
python3 examples/excavator/check_independent.py     # needs NumPy
```

`run.sh` runs these steps:

```
cylload step      excavator.step --level 1 --out step_level1.csv              # mass properties of the assemblies
cylload step      excavator.step --level 1 --bspline fan --out step_level1_fan.csv   # the same with the 1.1 triangle fans
cylload masstable excavator.step excavator_assign.csv --level 1 --remainder 60 --out excavator_masses.csv
cylload sweep     excavator.mech S1 --n 91 --out sweep_S1.csv                  # S2, S3 likewise
cylload ramp      excavator.mech S1 --flow 120 --plim 250 --out ramp_S1.csv
cylload simulate  excavator.mech S1 --flow 120 --ramp 0.05 --dir ext --out simulate_S1.csv
cylload map       excavator.mech S1 S2 --report S1 --flow 120 --out map_boom_arm.csv
cylload mc        excavator.mech excavator_mc.txt --sensitivity sensitivity.csv --out monte_carlo.csv
cylload draw      excavator.mech --step excavator.step --parts excavator_parts.csv --level 1 --color yes --out excavator.svg
```

`run.sh` checks that the mass table it writes is identical to `excavator_masses.csv`.

## Expected results

Mass properties of the 13 assemblies against the exact values (`results/check_mass_properties.csv`):

| B-spline faces | largest volume error | largest centroid error | largest error of J_z |
|:--|--:|--:|--:|
| integrated in their parameter domain (1.2, default) | 1.6 × 10⁻⁸ | 5 × 10⁻⁶ mm | 1.1 × 10⁻⁸ |
| triangle fans on the boundary (1.1, `--bspline fan`) | 6.0 % | 13.6 mm | 5.6 % |

The bucket, with 12 B-spline faces, has the largest fan errors; the volume and moment
consistency checks of the fan variant stay small, so they do not reveal these errors.

Loads over the strokes (`results/sweep_S*.txt`; the three force routes agree to 6 × 10⁻¹⁰):

| Cylinder | bore/rod, mm | stroke, mm | largest pressure, bar | reduced mass, t | natural frequency, Hz |
|:--|:--|:--|--:|:--|:--|
| S1 boom | 125/85 | 1300–2200 | 82.9 | 49–124 | 3.6–13.7 |
| S2 arm | 115/80 | 1200–1950 | 72.0 | 23–76 | 4.5–16.2 |
| S3 bucket | 100/70 | 1050–1500 | 101.7 | 7.5–24 | 6.7–14.1 |

Comparison with Exudyn 1.13.6 (`check_exudyn.py`, `results/check_exudyn.txt`; needs
`pip install exudyn numpy`): the mechanism rebuilt from Exudyn planar rigid bodies, revolute
joints and distance constraints gives the same cylinder forces to 5 × 10⁻¹⁰, reduced masses to
1.1 × 10⁻⁸ and natural frequencies to 5 × 10⁻⁹ at all 258 poses.

Independent check (`results/check_independent.txt`): a separate NumPy implementation in
absolute coordinates with Lagrange multipliers gives the same forces, reduced masses and
natural frequencies at all 258 poses of the three sweeps to within 5 × 10⁻¹⁰.

Start-up of the boom cylinder at 120 L/min with a 0.05 s ramp: simulated excess force
219.1 kN against 219.0 kN from the ramp formula (+0.02 %). Monte Carlo (N = 3000, 95 %
bootstrap intervals): P95 of the boom-cylinder pressure in the reference pose 79.1 bar
(78.7–79.4 bar).

The drawings `results/excavator.svg` and `results/excavator_raised.svg` show every part of
the STEP model moved into the computed pose.

## What the example does not show

The masses, dimensions and data-sheet values are invented and have not been compared with
a real machine. The example demonstrates and tests the workflow; it is not a design
reference for an excavator.
