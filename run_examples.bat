@echo off
rem Runs the Horty example and writes all results to examples\horty\results
setlocal
if "%CYLLOAD%"=="" set CYLLOAD=build\Release\cylload.exe
set E=examples\horty
set R=%E%\results
if not exist %R% mkdir %R%
for %%V in ("" "_outline") do call :one %%~V
"%CYLLOAD%" sweep %E%\tower.mech winch S8=973 S10=758 --from 4942.5 --to 10721.5 --n 59 --report S10 --out %R%\hammer_position.csv 2> %R%\hammer_position.txt
"%CYLLOAD%" mc %E%\tower.mech %E%\tower_mc.txt --out %R%\monte_carlo.csv 2> %R%\monte_carlo.txt
"%CYLLOAD%" mc %E%\tower.mech %E%\tower_mc_110kg.txt --out %R%\monte_carlo_110kg.csv 2> %R%\monte_carlo_110kg.txt
"%CYLLOAD%" sweep %E%\tower_outline.mech winch S8=973 S10=758 --from 4942.5 --to 10721.5 --n 59 --report S10 --out %R%\hammer_position_outline.csv 2> %R%\hammer_position_outline.txt
"%CYLLOAD%" mc %E%\tower.mech %E%\tower_mc.txt --sensitivity %R%\sensitivity.csv > NUL 2>&1
"%CYLLOAD%" simulate %E%\tower.mech S10 S8=973 S10=978 --flow 11.25 --ramp 0.1 --out %R%\simulation_tr0.1.csv 2> %R%\simulation_tr0.1.txt
"%CYLLOAD%" simulate %E%\tower.mech S10 S8=973 S10=978 --flow 11.25 --ramps 0.0001,0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.22,0.24,0.26,0.28,0.3,0.32,0.34,0.36,0.38,0.4,0.45,0.5,0.55,0.6,0.7,0.8,0.9,1.0 --out %R%\ramp_scan.csv 2> %R%\ramp_scan.txt
"%CYLLOAD%" simulate %E%\tower.mech S10 S8=973 S10=978 --flow 11.25 --ramps 0.0001,0.005,0.01,0.015,0.02,0.025,0.03,0.035,0.04,0.045,0.05,0.055,0.06,0.065,0.07,0.075,0.08,0.085,0.09,0.095,0.1,0.105,0.11,0.115,0.12,0.125,0.13,0.135,0.14,0.145,0.15,0.155,0.16,0.165,0.17,0.175,0.18,0.185,0.19,0.195,0.2,0.205,0.21,0.215,0.22,0.225,0.23,0.235,0.24,0.245,0.25,0.255,0.26,0.265,0.27,0.275,0.28,0.285,0.29,0.295,0.3,0.305,0.31,0.315,0.32,0.325,0.33,0.335,0.34,0.345,0.35,0.355,0.36,0.365,0.37,0.375,0.38,0.385,0.39,0.395,0.4,0.42,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1 --out %R%\ramp_scan_dense.csv 2> %R%\ramp_scan_dense.txt
"%CYLLOAD%" map %E%\tower.mech S10 winch S8=973 --yfrom 4942.5 --yto 10721.5 --nx 177 --ny 59 --flow 56.1 --out %R%\map_tilt_hammer.csv 2> %R%\map_tilt_hammer.txt
"%CYLLOAD%" draw %E%\tower.mech S8=973 S10=1183 --width 80 --height 110 --out %R%\mechanism_working.svg
"%CYLLOAD%" draw %E%\tower.mech S8=698 S10=758 --width 120 --height 60 --out %R%\mechanism_transport.svg
"%CYLLOAD%" plot %E%\figures\fig_static.txt --out %R%\fig_static.svg
"%CYLLOAD%" plot %E%\figures\fig_verification.txt --out %R%\fig_verification.svg
"%CYLLOAD%" draw %E%\tower.mech S8=973 S10=1183 --scale 20 --ybottom -800 --height 120 --color yes --font_pt 8 --view -2350,-800,1050,4300 --title "Working position" --out %R%\mechanism_working_colour.svg
"%CYLLOAD%" draw %E%\tower.mech S8=698 S10=758 --scale 20 --ybottom -800 --height 120 --color yes --font_pt 8 --items no --view -2600,-800,450,2000 --title "Transport position" --out %R%\mechanism_transport_colour.svg
"%CYLLOAD%" mc %E%\tower.mech %E%\tower_mc_hammer.txt --sensitivity %R%\sensitivity_hammer.csv > nul 2>&1
python -c "import sys;R=sys.argv[1];rows=[l.strip().split(',') for l in open(R+'/sensitivity_hammer.csv') if l.startswith('#stats,sb_')];f=open(R+'/monte_carlo_hammer.csv','w');f.write('s_b_m,p_min_bar,p5_bar,p50_bar,p95_bar,p_max_bar,mean_bar\n');[f.write(','.join([r[1][3:]]+r[2:8])+'\n') for r in rows]" %R% || echo (python not found: monte_carlo_hammer.csv not regenerated)
set U=examples\usarya
set RU=%U%\results
if not exist %RU% mkdir %RU%
for %%M in (0 6000) do for %%T in (0 3.1325 10.1929) do (
  "%CYLLOAD%" sweep %U%\cases\lift_m%%M_b%%T.mech lift --from 1150.388 --to 1890.388 --n 75 --out %RU%\lift_m%%M_b%%T.csv 2> %RU%\lift_m%%M_b%%T.txt
  "%CYLLOAD%" sweep %U%\cases\tip_m%%M_b%%T.mech tip --from 988.435 --to 1588.435 --n 75 --out %RU%\tip_m%%M_b%%T.csv 2> %RU%\tip_m%%M_b%%T.txt
)
"%CYLLOAD%" sweep %U%\rotor.mech S --from 1013 --to 1583 --n 59 --out %RU%\rotor.csv 2> %RU%\rotor.txt
"%CYLLOAD%" sweep %U%\drawbar.mech ED --from 737 --to 935 --n 45 --out %RU%\drawbar.csv 2> %RU%\drawbar.txt
"%CYLLOAD%" sweep %U%\drawbar_load.mech ED --from 737 --to 935 --n 45 --out %RU%\drawbar_load.csv 2> %RU%\drawbar_load.txt
"%CYLLOAD%" pose %U%\drawbar_load.mech ED=737 > %RU%\drawbar_load_737.txt
"%CYLLOAD%" plot %U%\figures\fig_lift.txt --out %RU%\fig_lift.svg
"%CYLLOAD%" plot %U%\figures\fig_tip.txt --out %RU%\fig_tip.svg
"%CYLLOAD%" plot %U%\figures\fig_fold.txt --out %RU%\fig_fold.svg
"%CYLLOAD%" draw %U%\hopper.mech lift=1150.388 tip=988.435 --sketch %U%\hopper_sketch.txt --scale 10 --ybottom -1760 --color yes --font_pt 7 --items no --view -2200,-1760,1750,4420 --out %RU%\hopper_lowered.svg
"%CYLLOAD%" draw %U%\hopper.mech lift=1890.388 tip=988.435 --sketch %U%\hopper_sketch.txt --scale 10 --ybottom -1760 --color yes --font_pt 7 --items no --view -2250,-1760,1750,4420 --out %RU%\hopper_raised.svg
"%CYLLOAD%" draw %U%\hopper.mech lift=1890.388 tip=1588.435 --sketch %U%\hopper_sketch.txt --scale 10 --ybottom -1760 --color yes --font_pt 7 --items no --view -3150,-1760,1750,4420 --out %RU%\hopper_tipped.svg
rem version 1.2: reduced mass without the cylinder pin masses, damping, Monte Carlo convergence
for %%V in ("" "_outline") do (
  "%CYLLOAD%" ramp  %E%\tower%%~V_mr_nopins.mech S10 S8=973 --n 177 --flow 56.1 --plim 140 --out %R%\ramp_tilt%%~V_mr_nopins.csv 2> %R%\ramp_tilt%%~V_mr_nopins.txt
  "%CYLLOAD%" sweep %E%\tower%%~V_mr_nopins.mech S10 S8=973 --n 177 --out %R%\tilt%%~V_mr_nopins.csv 2> %R%\tilt%%~V_mr_nopins.txt
)
for %%Z in (0 0.05 0.1 0.2) do "%CYLLOAD%" simulate %E%\tower.mech S10 S8=973 S10=978 --flow 11.25 --damping %%Z --ramps 0.0001,0.05,0.1,0.2,0.3,0.5,0.75,1.0 --out %R%\ramp_scan_zeta%%Z.csv 2> nul
for %%N in (500 1000 2000 4000 8000) do "%CYLLOAD%" mc %E%\tower.mech %E%\tower_mc.txt --samples %%N --bootstrap 200 --sensitivity %R%\sensitivity_N%%N.csv > nul 2>&1
for %%H in (0000 0500 1000 1500 2000) do (
  "%CYLLOAD%" sweep %U%\cases\drawbar_frame_h%%H.mech ED --from 737 --to 935 --n 45 --out %RU%\drawbar_frame_h%%H.csv 2> %RU%\drawbar_frame_h%%H.txt
  "%CYLLOAD%" pose %U%\cases\drawbar_frame_h%%H.mech ED=737 > %RU%\drawbar_frame_h%%H_737.txt
)
call :excavator
echo results written to %R%
goto :eof
:excavator
set X=examples\excavator
set RX=%X%\results
if not exist %RX% mkdir %RX%
"%CYLLOAD%" step %X%\excavator.step --level 1 --out %RX%\step_level1.csv 2> %RX%\step_level1.txt
"%CYLLOAD%" step %X%\excavator.step --level 1 --bspline fan --out %RX%\step_level1_fan.csv 2> %RX%\step_level1_fan.txt
"%CYLLOAD%" masstable %X%\excavator.step %X%\excavator_assign.csv --level 1 --remainder 60 --out %RX%\excavator_masses.csv 2> %RX%\masstable.txt
"%CYLLOAD%" info  %X%\excavator.mech > %RX%\info.txt
"%CYLLOAD%" pose  %X%\excavator.mech > %RX%\pose_reference.txt
"%CYLLOAD%" sweep %X%\excavator.mech S1 --n 91 --out %RX%\sweep_S1.csv 2> %RX%\sweep_S1.txt
"%CYLLOAD%" sweep %X%\excavator.mech S2 --n 76 --out %RX%\sweep_S2.csv 2> %RX%\sweep_S2.txt
"%CYLLOAD%" sweep %X%\excavator.mech S3 --n 91 --out %RX%\sweep_S3.csv 2> %RX%\sweep_S3.txt
"%CYLLOAD%" ramp %X%\excavator.mech S1 --n 91 --flow 120 --plim 250 --out %RX%\ramp_S1.csv 2> %RX%\ramp_S1.txt
"%CYLLOAD%" simulate %X%\excavator.mech S1 --flow 120 --ramp 0.05 --dir ext --out %RX%\simulate_S1.csv 2> %RX%\simulate_S1.txt
"%CYLLOAD%" map %X%\excavator.mech S1 S2 --report S1 --nx 46 --ny 39 --flow 120 --out %RX%\map_boom_arm.csv 2> %RX%\map_boom_arm.txt
"%CYLLOAD%" mc %X%\excavator.mech %X%\excavator_mc.txt --out %RX%\monte_carlo.csv --sensitivity %RX%\sensitivity.csv 2> %RX%\monte_carlo.txt
"%CYLLOAD%" draw %X%\excavator.mech --step %X%\excavator.step --parts %X%\excavator_parts.csv --level 1 --color yes --scale 12 --ybottom -1150 --height 72 --font_pt 9 --items no --view -1400,-1150,4950,4300 --out %RX%\excavator.svg
"%CYLLOAD%" draw %X%\excavator.mech S1=2100 S2=1300 S3=1450 --step %X%\excavator.step --parts %X%\excavator_parts.csv --level 1 --color yes --scale 12 --ybottom -1150 --height 72 --font_pt 9 --items no --view -1400,-1150,4950,4300 --out %RX%\excavator_raised.svg
set K=examples\cyklop
set RK=%K%\results
if not exist %RK% mkdir %RK%
"%CYLLOAD%" map   %K%\cyklop.mech S1 S2 --report S1 --nx 21 --ny 21 --plim 100 --out %RK%\map_S1.csv 2> %RK%\map_S1.txt
"%CYLLOAD%" map   %K%\cyklop.mech S1 S2 --report S2 --nx 21 --ny 21 --plim 100 --out %RK%\map_S2.csv 2> %RK%\map_S2.txt
"%CYLLOAD%" sweep %K%\cyklop.mech S1 S2=902 --n 21 --report S1 --out %RK%\sweep_S1.csv 2> %RK%\sweep_S1.txt
"%CYLLOAD%" draw  %K%\cyklop.mech --overlay S1=1002,S2=1532 --font_pt 9 --out %RK%\cyklop.svg
python examples\excavator\check_independent.py || echo (python with numpy not found: independent check not run)
goto :eof
:one
set V=%1
"%CYLLOAD%" pose  %E%\tower%V%.mech S8=973 S10=1183 > %R%\pose_reference%V%.txt
"%CYLLOAD%" pose  %E%\tower%V%.mech S8=973 S10=758  > %R%\pose_mast_horizontal%V%.txt
"%CYLLOAD%" sweep %E%\tower%V%.mech S8  S10=1183 --n 141 --out %R%\lift%V%.csv 2> %R%\lift%V%.txt
"%CYLLOAD%" sweep %E%\tower%V%.mech S10 S8=973   --n 177 --out %R%\tilt%V%.csv 2> %R%\tilt%V%.txt
"%CYLLOAD%" sweep %E%\tower%V%.mech S19 S8=973 S10=1183 --n 101 --out %R%\slide%V%.csv 2> %R%\slide%V%.txt
"%CYLLOAD%" sweep %E%\post_lift%V%.mech Sp --n 121 --out %R%\post_lift%V%.csv 2> %R%\post_lift%V%.txt
"%CYLLOAD%" ramp  %E%\tower%V%.mech S10 S8=973 --n 177 --flow 56.1 --plim 140 --out %R%\ramp_tilt%V%.csv 2> %R%\ramp_tilt%V%.txt
"%CYLLOAD%" ramp  %E%\post_lift%V%.mech Sp --n 121 --flow 56.1 --plim 140 --out %R%\ramp_post_lift%V%.csv 2> %R%\ramp_post_lift%V%.txt
goto :eof
