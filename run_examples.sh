#!/bin/sh
# Runs the Horty example and writes all results to examples/horty/results
set -e
B=${CYLLOAD:-build/cylload}
E=examples/horty
R=$E/results
mkdir -p $R
for V in "" _outline; do
  T=$E/tower$V.mech; P=$E/post_lift$V.mech
  $B pose  $T S8=973 S10=1183                       > $R/pose_reference$V.txt
  $B pose  $T S8=973 S10=758                        > $R/pose_mast_horizontal$V.txt
  $B sweep $T S8  S10=1183 --n 141 --out $R/lift$V.csv      2> $R/lift$V.txt
  $B sweep $T S10 S8=973   --n 177 --out $R/tilt$V.csv      2> $R/tilt$V.txt
  $B sweep $T S19 S8=973 S10=1183 --n 101 --out $R/slide$V.csv 2> $R/slide$V.txt
  $B sweep $P Sp --n 121 --out $R/post_lift$V.csv           2> $R/post_lift$V.txt
  $B ramp  $T S10 S8=973 --n 177 --flow 56.1 --plim 140 --out $R/ramp_tilt$V.csv 2> $R/ramp_tilt$V.txt
  $B ramp  $P Sp --n 121 --flow 56.1 --plim 140 --out $R/ramp_post_lift$V.csv   2> $R/ramp_post_lift$V.txt
done
# hammer position along the mast (winch length 4942.5..10721.5 mm = 0..5.78 m), mast horizontal
$B sweep $E/tower.mech winch S8=973 S10=758 --from 4942.5 --to 10721.5 --n 59 --report S10 --out $R/hammer_position.csv 2> $R/hammer_position.txt
$B mc $E/tower.mech $E/tower_mc.txt --out $R/monte_carlo.csv 2> $R/monte_carlo.txt
$B mc $E/tower.mech $E/tower_mc_110kg.txt --out $R/monte_carlo_110kg.csv 2> $R/monte_carlo_110kg.txt
$B sweep $E/tower_outline.mech winch S8=973 S10=758 --from 4942.5 --to 10721.5 --n 59 --report S10 --out $R/hammer_position_outline.csv 2> $R/hammer_position_outline.txt
$B mc $E/tower.mech $E/tower_mc.txt --sensitivity $R/sensitivity.csv > /dev/null 2>&1
$B simulate $E/tower.mech S10 S8=973 S10=978 --flow 11.25 --ramp 0.1 --out $R/simulation_tr0.1.csv 2> $R/simulation_tr0.1.txt
$B simulate $E/tower.mech S10 S8=973 S10=978 --flow 11.25 --ramps 0.0001,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0 --out $R/ramp_scan.csv 2> $R/ramp_scan.txt
$B simulate $E/tower.mech S10 S8=973 S10=978 --flow 11.25 --ramps 0.0001,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08,0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12,0.125,0.13,0.135,0.14,0.145,0.15,0.155,0.16,0.165,0.17,0.175,0.18,0.185,0.19,0.195,0.2,0.205,0.21,0.215,0.22,0.225,0.23,0.235,0.24,0.245,0.25,0.255,0.26,0.265,0.27,0.275,0.28,0.285,0.29,0.295,0.3,0.305,0.31,0.315,0.32,0.325,0.33,0.335,0.34,0.345,0.35,0.355,0.36,0.365,0.37,0.375,0.38,0.385,0.39,0.395,0.4,0.42,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1 --out $R/ramp_scan_dense.csv 2> $R/ramp_scan_dense.txt
$B map $E/tower.mech S10 winch S8=973 --yfrom 4942.5 --yto 10721.5 --nx 177 --ny 59 --flow 56.1 --out $R/map_tilt_hammer.csv 2> $R/map_tilt_hammer.txt
$B draw $E/tower.mech S8=973 S10=1183 --width 80 --height 110 --out $R/mechanism_working.svg
$B draw $E/tower.mech S8=698 S10=758 --width 120 --height 60 --out $R/mechanism_transport.svg
$B draw $E/tower.mech S8=973 S10=1183 --scale 20 --ybottom -800 --height 120 --color yes --font_pt 8 --view -2350,-800,1050,4300 --title "Working position" --out $R/mechanism_working_colour.svg
$B draw $E/tower.mech S8=698 S10=758 --scale 20 --ybottom -800 --height 120 --color yes --font_pt 8 --items no --view -2600,-800,450,2000 --title "Transport position" --out $R/mechanism_transport_colour.svg
$B plot $E/figures/fig_static.txt --out $R/fig_static.svg
$B plot $E/figures/fig_verification.txt --out $R/fig_verification.svg
$B mc $E/tower.mech $E/tower_mc_hammer.txt --sensitivity $R/sensitivity_hammer.csv > /dev/null 2>&1
python3 - "$R" <<'PY' || echo "(python3 not found: monte_carlo_hammer.csv not regenerated)"
import sys
R=sys.argv[1]
rows=[l.strip().split(',') for l in open(R+'/sensitivity_hammer.csv') if l.startswith('#stats,sb_')]
with open(R+'/monte_carlo_hammer.csv','w') as f:
    f.write('s_b_m,p_min_bar,p5_bar,p50_bar,p95_bar,p_max_bar,mean_bar\n')
    for r in rows: f.write(','.join([r[1][3:]]+r[2:8])+'\n')
PY
U=examples/usarya; RU=$U/results
mkdir -p $RU
for f in $U/cases/lift_m*.mech; do $B sweep $f lift --from 1150.388 --to 1890.388 --n 75 --out $RU/$(basename $f .mech).csv 2> $RU/$(basename $f .mech).txt; done
for f in $U/cases/tip_m*.mech; do $B sweep $f tip --from 988.435 --to 1588.435 --n 75 --out $RU/$(basename $f .mech).csv 2> $RU/$(basename $f .mech).txt; done
$B sweep $U/rotor.mech S --from 1013 --to 1583 --n 59 --out $RU/rotor.csv 2> $RU/rotor.txt
$B sweep $U/drawbar.mech ED --from 737 --to 935 --n 45 --out $RU/drawbar.csv 2> $RU/drawbar.txt
$B sweep $U/drawbar_load.mech ED --from 737 --to 935 --n 45 --out $RU/drawbar_load.csv 2> $RU/drawbar_load.txt
$B pose $U/drawbar_load.mech ED=737 > $RU/drawbar_load_737.txt
$B plot $U/figures/fig_lift.txt --out $RU/fig_lift.svg
$B plot $U/figures/fig_tip.txt --out $RU/fig_tip.svg
$B plot $U/figures/fig_fold.txt --out $RU/fig_fold.svg
$B draw $U/hopper.mech lift=1150.388 tip=988.435 --sketch $U/hopper_sketch.txt --scale 10 --ybottom -1760 --color yes --font_pt 7 --items no --view -2200,-1760,1750,4420 --out $RU/hopper_lowered.svg
$B draw $U/hopper.mech lift=1890.388 tip=988.435 --sketch $U/hopper_sketch.txt --scale 10 --ybottom -1760 --color yes --font_pt 7 --items no --view -2250,-1760,1750,4420 --out $RU/hopper_raised.svg
$B draw $U/hopper.mech lift=1890.388 tip=1588.435 --sketch $U/hopper_sketch.txt --scale 10 --ybottom -1760 --color yes --font_pt 7 --items no --view -3150,-1760,1750,4420 --out $RU/hopper_tipped.svg
# version 1.2: reduced mass without the cylinder pin masses (as in the published Python model)
for V in "" _outline; do
  $B ramp  $E/tower${V}_mr_nopins.mech S10 S8=973 --n 177 --flow 56.1 --plim 140 --out $R/ramp_tilt${V}_mr_nopins.csv 2> $R/ramp_tilt${V}_mr_nopins.txt
  $B sweep $E/tower${V}_mr_nopins.mech S10 S8=973 --n 177 --out $R/tilt${V}_mr_nopins.csv 2> $R/tilt${V}_mr_nopins.txt
done
# damping: excess force over the quasi-steady force for zeta = 0, 0.05, 0.1, 0.2
for Z in 0 0.05 0.1 0.2; do
  $B simulate $E/tower.mech S10 S8=973 S10=978 --flow 11.25 --damping $Z --ramps 0.0001,0.05,0.1,0.2,0.3,0.5,0.75,1.0 --out $R/ramp_scan_zeta$Z.csv 2> /dev/null
done
# convergence of the Monte Carlo estimators with the number of samples
for N in 500 1000 2000 4000 8000; do
  $B mc $E/tower.mech $E/tower_mc.txt --samples $N --bootstrap 200 --sensitivity $R/sensitivity_N$N.csv > /dev/null 2>&1
done
# drawbar loaded by the frame weight at its centre of gravity (five heights of the centre of gravity)
for f in $U/cases/drawbar_frame_h*.mech; do
  $B sweep $f ED --from 737 --to 935 --n 45 --out $RU/$(basename $f .mech).csv 2> $RU/$(basename $f .mech).txt
  $B pose $f ED=737 > $RU/$(basename $f .mech)_737.txt
done
echo "results written to $R and $RU"
sh examples/excavator/run.sh
sh examples/cyklop/run.sh
python3 examples/excavator/check_independent.py || echo "(python3 with numpy not found: independent check not run)"
