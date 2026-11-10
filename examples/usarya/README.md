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

* `drawbar_load.mech` contains two forces on the frame. `hitchload` (75 046.5 N downwards at
  the hitch point C) is the hitch load V g = 7650 kg × 9.81 m/s²; only its moment about the
  axle pivot A matters, and it stands for the weight of everything the frame carries.
  `drawbar_on_frame` (3090.15 N upwards at the drawbar centre of gravity) is the drawbar
  weight (315 kg) taken back off the frame, because V already includes the drawbar, whose
  weight is applied separately to the drawbar link.
* `cases/drawbar_frame_h*.mech` replace this equivalent load by the physical one: the
  frame carries its weight (12 400 kg machine mass less the 315 kg drawbar) at its centre
  of gravity, and the axle and hitch reactions follow from the kinematics. The position of
  the centre of gravity along the frame is set so that the hitch load is 7650 kg at
  ED = 935 mm (`make_drawbar_frame.py`); its height h above the frame beam is not known and
  is varied from 0 to 2000 mm. At ED = 737 mm the pressure is then 190.2–201.4 bar and the
  hitch load 7483–7915 kg (h = 0–2000 mm), against 194.5 bar in the published analysis
  (constant hitch load) and 186.0 bar with the constant-force model of `drawbar_load.mech`.
  The published value corresponds to h ≈ 770 mm. The difference between the published
  value and `drawbar_load.mech` is therefore a sensitivity to the load model, not an
  error of either analysis.

`run_examples.sh` (or `.bat`) writes the sweeps and the figures to `results/`.
