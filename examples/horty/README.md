# Example: Horty AutoHammer orchard pile driver

This is the positioning mechanism of the tower of the Horty AutoHammer, a tracked
machine that drives orchard posts. Element numbers follow patent application
P.443257. The example has three parts:

* `tower.mech` — the tower with four actuators:
  * lifting parallelogram A–B–C–D with cylinder S8 (100/50, 698–1048 mm);
  * tower tilted about H by cylinder S10 (80/40, 758–1198 mm);
  * mast sliding along the tower, driven by cylinder S19 (60/40, 1262–2262 mm);
  * hammer position along the mast, set by a winch. The winch length is
    4942.5 mm + s_b, where the hammer position s_b runs from 0 to 5.78 m and
    s_b = 2.90 m is the position in the solid model.
* `post_lift.mech` — the post-lifting beam 2.5.0 on the carriage with cylinder Sp
  (50/30, 630–1030 mm). Variants with 56 kg (intermediate) and 110 kg (end) posts.
* `tower_mc.txt`, `tower_mc_110kg.txt` — parameter uncertainty at fixed operating states.

The mass tables were produced by `cylload masstable` from the STEP models of the tower
and the carriage (`tower_assign.csv`, `post_lift_assign.csv`). They use the title-block
masses, and the 222.7 kg of unidentified mass is distributed in proportion. There are two
mass-property variants:

* `tower.mech`, `post_lift.mech` — centroids and moments of inertia integrated from the
  solid model (`--properties solid`);
* `tower_outline.mech`, `post_lift_outline.mech` — centres of the bounding outlines and
  rectangular-prism moments of inertia (`--properties outline`), the simplified approach
  available when only outlines are known.

The STEP files and the design documentation are not distributed. The simplified input
data are published with the permission of Arno Group.

## Results compared with the published analysis

Published values are from Tucki K., Nowakowski M., *Applied Sciences* (applsci-4655781),
Tables 4, 5 and A2. Static values agree to the published rounding.

| Quantity | outline, paper | outline, CylLoad | solid, paper | solid, CylLoad |
|---|---|---|---|---|
| p10, mast horizontal, s_b = 2.90 m, bar | 101.1 | 101.14 | 108.4 | 108.42 |
| p10, mast vertical (rod side), bar | 36.8 | 36.79 | 41.4 | 41.39 |
| p10, end of stroke (rod side), bar | 47.6 | 47.64 | 53.1 | 53.13 |
| p8, tower lowered, bar | 61.7 | 61.72 | 61.8 | 61.84 |
| pp, beam nearly horizontal, 56 kg post, bar | 64.1 | 64.10 | 58.0 | 57.94 |
| s_b at 140 bar, mast horizontal, m | 4.28 | 4.28 | 4.02 | 4.02 |
| J_H about joint H, kg·m² | 4412 | 4413 | 4612 | 4647 ¹ |
| f_h of cylinder 10, Hz | 3.7–4.5 | 3.70–4.46 | 3.6–4.4 | 3.61–4.34 |
| f_h of cylinder p, Hz | 2.2–6.0 | 2.25–6.04 | 2.2–5.9 | 2.21–5.93 |
| t_r* tilting, 56.1 L/min, 140 bar, s ² | 2.9 | 2.97 | 3.80 | 3.84 |
| t_r* post lifting 56 / 110 kg, s ² | 4.9 / 7.1 | 4.93 / 7.14 | 5.22 / 7.35 | 5.24 / 7.37 |

¹ The whole difference of 35 kg·m² comes from the moment of inertia of the mast about its
own centre of gravity: 767.8 kg·m² here against 740.2 kg·m² in Table A1 of the paper. Both
are scaled by the factor 1.2508 for the remaining mass, so the difference in J_H is
(767.8 − 740.2) × 1.2508 = 34.5 kg·m². The mast term itself is only a small part of J_H;
the difference is a difference in that term, not the term. CylLoad's value is checked by a
second integration field, agreeing to 1e-10, and lies close to the prism estimate of
766.5 kg·m².

² Retraction velocity, 56.1 L/min. The reduced mass in CylLoad also contains the lumped
pins of cylinder 10 (about 0.3 %), which the published Python model neglected. Without them
(`tower_mr_nopins.mech`, `tower_outline_mr_nopins.mech`, `reduced_mass_pins=no`) the natural
frequencies become 3.61–4.35 Hz (solid; published 3.6–4.4) and 3.71–4.46 Hz (outline), and
t_r* becomes 3.83 s (solid; published 3.80) and 2.95 s (outline; published 2.9). The ramp
time is governed by the mast-horizontal end of the stroke (L₁₀ = 758 mm), where
t_r* = 2 m_r v / ((p_lim − p_st) A) magnifies small differences: the remaining 1.8 % of the
outline value corresponds to 0.7 bar of static pressure. The values are therefore given to
within 2 %, not "to the published rounding".

The Monte Carlo percentiles (`results/monte_carlo*.txt`) match Table 6 of the paper
within the sampling error (N = 3000).

Run `../../run_examples.sh` (or `run_examples.bat`) from the repository root to
regenerate `results/`.
