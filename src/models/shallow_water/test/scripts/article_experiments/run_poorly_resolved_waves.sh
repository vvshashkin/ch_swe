#!/bin/bash

source $(dirname "$0")"/gen_config_linear.sh"

EXE=$1/SHALLOW_WATER_MAIN
Nprocs=$2

OMEGA='1e-4'

CONFIG_TEMPLATE="
domain = {\n
    N = %%%N,\n
    Nz = 1,\n
    topology_type         = 'cube',\n
    horizontal_staggering = '%%%STAGGERING',\n
    metric                = {\n
            metric_type             = 'shallow_atmosphere_metric',\n
            metric_2d_type          = 'ecs',\n
            planet_radius           = Earth_radii,\n
            omega                   = $OMEGA,
            coriolis_constant       = .true.
}},\n
\n
dt = 1.0*%%%DT*sec,\n
simulation_time = 25.0*day,\n
tau_write = 1.0*hour,\n
tau_diagnostics = 1.0*hour,\n
\n
timescheme_name = 'rk4'\n
\n
testcase_section = 'Gauss_linear_config'\n
Gauss_linear_config = {\n
    testcase_name ='Oscillating_gauss_linear',
    nu = 32
}\n
"
HMEAN="(2*pi*Earth_radii\/(5*day))^2\/Earth_grav"

run_gaussian_linear(){
    CONFIG_TEMPLATE=$1
    N=$2
    DT=$3
    STAGGERING=$4
    SBP_ORDER=$5

    SUFFIX="N"$N"_dt"$DT"_"$STAGGERING"_"$SBP_ORDER$corner
    OUTFILE="swm_"$SUFFIX".out"
    HDAT="h_"$SUFFIX".dat"

    echo $SUFFIX
    echo $OUTFILE

	gen_config "$CONFIG_TEMPLATE" $N $DT $STAGGERING $SBP_ORDER $HMEAN > gaussian.conf
	mpirun -n $Nprocs $EXE gaussian.conf &> $OUTFILE
    # mv h.dat $HDAT
}

run_gaussian_linear "$CONFIG_TEMPLATE" 064 0450 Ch 21 $HMEAN

run_gaussian_linear "$CONFIG_TEMPLATE" 064 0450 Ch 42 $HMEAN

run_gaussian_linear "$CONFIG_TEMPLATE" 064 0450 Ch 63 $HMEAN

run_gaussian_linear "$CONFIG_TEMPLATE" 064 0450 Ah 63 $HMEAN