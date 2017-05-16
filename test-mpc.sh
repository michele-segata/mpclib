#!/usr/bin/env bash

MPC=build/src/test-mpc

if [ ! -f "$MPC" ]; then
	echo "Please build the project before running the samples"
	exit 1
fi

r_path=`which Rscript`
if [ "$r_path" == "" ]; then
	echo "R not found. Samples will be run but not plots will be produced"
	r_found=false
else
	r_found=true
fi

echo "Running test 01"
$MPC > test01.csv
if [ $r_found ]; then
	Rscript plot-mpc.R test01.csv &> /dev/null
fi

echo "Running test 02"
$MPC -u 0.01 -d 0.01 > test02.csv
if [ $r_found ]; then
	Rscript plot-mpc.R test02.csv &> /dev/null
fi

echo "Running test 03"
$MPC -u 0.01 -d 0.01 --u-min -1.0 --u-max 1.0 > test03.csv
if [ $r_found ]; then
	Rscript plot-mpc.R test03.csv &> /dev/null
fi

echo "Running test 04"
$MPC -u 0.01 -d 0.01 --u-min -1.0 --u-max 1.0 --du-min -0.5 --du-max 0.5 > test04.csv
if [ $r_found ]; then
	Rscript plot-mpc.R test04.csv &> /dev/null
fi

echo "Running test 05"
$MPC -u 0.01 -d 0.01 -b 0.6 > test05.csv
if [ $r_found ]; then
	Rscript plot-mpc.R test05.csv &> /dev/null
fi

echo "Running test 06"
$MPC -u 0.01 -d 0.01 -b 0.6 --out-slack 10 > test06.csv
if [ $r_found ]; then
	Rscript plot-mpc.R test06.csv 0.6 &> /dev/null
fi
