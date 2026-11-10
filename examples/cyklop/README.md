# T-214 grapple loader boom (model from manufacturer documentation)

The boom of the T-214/3 "Cyklop-1" trailed grapple loader (POL-MOT Warfama), a machine in series
production, described only from its operator's manual and spare parts catalogue (in Polish), without
a solid model. It shows the typical situation of a user who has drawings but no CAD data.

## Sources of the data

| Data | Source |
|---|---|
| Pins A, F and cylinder pins | read from Fig. 1 of the manual (scale drawing of the working range, 1 m grid); accuracy about ±30–50 mm |
| Cylinders S1 (lift) and S2 (arm) | parts catalogue, table 4 and table 29: CJ2F-12-80/45/630, bore 80 mm, rod 45 mm, stroke 630 mm |
| Rated load 500 kg at the end of the boom, including the grapple | manual, technical data |
| Relief valve 10 MPa, pump 40 l/min | manual, technical data and section 7.8.1 |
| Masses of the boom (200 kg) and arm (100 kg), cylinder masses (25 kg) | **estimated**; not given in the manual |

The reference pose is the pose of Fig. 1 of the manual (boom raised, arm straight), so the lift
cylinder is at its maximum length and the arm cylinder at its minimum length.

## Expected results (`sh run.sh`)

* Lifting height of the grapple bottom (suspension point minus 1.28 m): 4.02 m (manual: 4.0 m).
* Maximum pressures in the whole workspace at the rated load: S1 72.8 bar, S2 80.1 bar, below the
  relief setting of 100 bar. With arm masses 50 % higher: 80.6 and 83.7 bar.
* Extension and retraction times at 40 l/min: 4.75 s and 3.25 s (manual: 4–5 s per movement).

These are checks of consistency with manufacturer data, not a validation of the arm masses.
