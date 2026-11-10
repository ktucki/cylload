#!/bin/sh
# Open example: complete workflow STEP -> mass table -> forces, pressures, dynamics, Monte Carlo, drawing.
# Run from the repository root: sh examples/excavator/run.sh
set -e
B=${CYLLOAD:-build/cylload}
E=examples/excavator
R=$E/results
mkdir -p $R
# 1. mass properties of the STEP assemblies (level 1 = children of the root), exact B-spline integration
$B step $E/excavator.step --level 1 --out $R/step_level1.csv 2> $R/step_level1.txt
$B step $E/excavator.step --level 1 --bspline fan --out $R/step_level1_fan.csv 2> $R/step_level1_fan.txt
# 2. mass table: title-block masses, 60 kg unidentified remainder, data-sheet CoG of the bought-in valve block
$B masstable $E/excavator.step $E/excavator_assign.csv --level 1 --remainder 60 --out $R/excavator_masses.csv 2> $R/masstable.txt
cmp $R/excavator_masses.csv $E/excavator_masses.csv && echo "mass table reproduced" >> $R/masstable.txt
# 3. statics over the strokes (three force formulations, joint reactions, pressures, reduced mass, f_h)
$B info  $E/excavator.mech > $R/info.txt
$B pose  $E/excavator.mech > $R/pose_reference.txt
$B sweep $E/excavator.mech S1   --n 91 --out $R/sweep_S1.csv   2> $R/sweep_S1.txt
$B sweep $E/excavator.mech S2    --n 76 --out $R/sweep_S2.csv    2> $R/sweep_S2.txt
$B sweep $E/excavator.mech S3 --n 91 --out $R/sweep_S3.csv 2> $R/sweep_S3.txt
# 4. shortest valve ramp time and start-up simulation of the boom cylinder (lifting, 120 L/min)
$B ramp $E/excavator.mech S1 --n 91 --flow 120 --plim 250 --out $R/ramp_S1.csv 2> $R/ramp_S1.txt
$B simulate $E/excavator.mech S1 --flow 120 --ramp 0.05 --dir ext --out $R/simulate_S1.csv 2> $R/simulate_S1.txt
# 5. load map over the boom and arm strokes (bucket cylinder at reference length)
$B map $E/excavator.mech S1 S2 --report S1 --nx 46 --ny 39 --flow 120 --out $R/map_boom_arm.csv 2> $R/map_boom_arm.txt
# 6. Monte Carlo with bootstrap confidence intervals
$B mc $E/excavator.mech $E/excavator_mc.txt --out $R/monte_carlo.csv --sensitivity $R/sensitivity.csv 2> $R/monte_carlo.txt
# 7. drawing with the solid model projected on the machine plane
$B draw $E/excavator.mech --step $E/excavator.step --parts $E/excavator_parts.csv --level 1 --color yes --scale 12 --ybottom -1150 --height 72 --font_pt 9 --items no --view -1400,-1150,4950,4300 --out $R/excavator.svg
$B draw $E/excavator.mech S1=2100 S2=1300 S3=1450 --step $E/excavator.step --parts $E/excavator_parts.csv --level 1 --color yes --scale 12 --ybottom -1150 --height 72 --font_pt 9 --items no --view -1400,-1150,4950,4300 --out $R/excavator_raised.svg
echo "results written to $R"
