# Stone picker (CleanField / Usarya SCS-100)

Four mechanisms of a trailed stone picker, used to check CylLoad against an independent
model and against published results:

| File | Mechanism | Reference |
|:--|:--|:--|
| `rotor.mech` | folding of the working module (cylinder 125/56 x 570) | parametric Mathcad model of the machine (closed-form triangle solution); pressures from USARYA 2 |
| `drawbar.mech` | drawbar linkage that sets the frame tilt (cylinder 100/50 x 220) | parametric Mathcad model (triangles ABC, BCD, BDE) |
| `cases/lift_m*_b*.mech` | hopper lifting parallelogram, load 0 or 6000 kg, frame tilt 0, 3.1325 or 10.1929 deg | pressures from USARYA 2 |
| `cases/tip_m*_b*.mech` | hopper tipping, same loads and tilts | pressures from USARYA 2 |
| `hopper.mech`, `hopper_sketch.txt` | hopper lifting and tipping together, with 1978 lines of assembly drawing 345-0.0.0 for `draw --sketch` | drawing only |

USARYA 2: K. Tucki, L. Mieszkalski, M. Gołębiewski, Loads on the hydraulic cylinders of a
trailed stone picker determined from the kinematics of its mechanisms, *Agriculture*
(accepted, agriculture-4643982).

* Two cylinders working in parallel are modelled as one cylinder with twice the piston
  area (bore multiplied by the square root of 2), so pressures are those of each cylinder.
* The frame tilt is introduced by rotating the gravity vector (`gravity` keyword).
* The published drawbar loads express the weight carried by the frame by a constant hitch load.
  `drawbar_load.mech` applies that weight as a constant force on the frame (`force` keyword),
  calibrated to the hitch load of 7650 kg at ED = 935 mm, where CylLoad gives 206.79 bar
  against the published 206.8 bar. At ED = 737 mm it gives 186.0 bar against 194.5 bar,
  because with a constant frame weight the hitch reaction falls to 7323 kg over the stroke
  (`results/drawbar_load_737.txt`), while the published analysis keeps it constant.

`run_examples.sh` (or `.bat`) writes the sweeps and the figures to `results/`.
