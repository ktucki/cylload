#!/bin/sh
# T-214 grapple loader boom: pressure maps of both cylinders, lift sweep with the arm straight, drawing
cd "$(dirname "$0")"
B=${CYLLOAD:-../../build/cylload}
R=results
mkdir -p $R
$B info cyklop.mech > $R/info.txt
$B map cyklop.mech S1 S2 --report S1 --nx 21 --ny 21 --plim 100 --out $R/map_S1.csv 2> $R/map_S1.txt
$B map cyklop.mech S1 S2 --report S2 --nx 21 --ny 21 --plim 100 --out $R/map_S2.csv 2> $R/map_S2.txt
$B sweep cyklop.mech S1 S2=902 --n 21 --report S1 --out $R/sweep_S1.csv 2> $R/sweep_S1.txt
$B draw cyklop.mech --overlay S1=1002,S2=1532 --font_pt 9 --out $R/cyklop.svg
