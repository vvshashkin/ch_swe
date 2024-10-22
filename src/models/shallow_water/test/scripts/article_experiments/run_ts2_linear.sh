#!/bin/bash

source $(dirname "$0")"/gen_config_linear.sh"

EXE=$1/SHALLOW_WATER_MAIN
Nprocs=$2

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
            omega                   = Earth_omega,\n
            rotation_axis          = [1.0/sqrt(2.0),0.0,1.0/sqrt(2.0)]\n
}},\n
\n
dt = 1.0*%%%DT*sec,\n
simulation_time = 10.0*day,\n
tau_write = 1.0*hour,\n
tau_diagnostics = 1.0*hour,\n
\n
timescheme_name = 'rk4'\n
\n
testcase_section = 'Williamson_test2_config'\n

Williamson_test2_config = {\n
    testcase_name ='Williamson_test2',\n
    linear_mode   = .true.\n
}\n"

HMEAN="29400\/Earth_grav"

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

	gen_config "$CONFIG_TEMPLATE" $N $DT $STAGGERING $SBP_ORDER $HMEAN > ts2_lin.conf
	mpirun -n $Nprocs $EXE ts2_lin.conf &> $OUTFILE
    mv h.dat $HDAT
}

# run_gaussian_linear "$CONFIG_TEMPLATE" 016 1800 Ch 21 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 048 0600 Ch 21 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 064 0450 Ch 21 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 096 0300 Ch 21 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 128 0225 Ch 21 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 192 0150 Ch 21 $HMEAN
# run_gaussian_linear "$CONFIG_TEMPLATE" 192 0150 Ch 21 $HMEAN

# run_gaussian_linear "$CONFIG_TEMPLATE" 016 1800 Ch 42 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 048 0600 Ch 42 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 064 0450 Ch 42 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 096 0300 Ch 42 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 128 0225 Ch 42 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 192 0150 Ch 42 $HMEAN
# run_gaussian_linear "$CONFIG_TEMPLATE" 192 0150 Ch 42 $HMEAN

# run_gaussian_linear "$CONFIG_TEMPLATE" 016 1800 Ch 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 048 0600 Ch 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 064 0450 Ch 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 096 0300 Ch 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 128 0225 Ch 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 192 0150 Ch 63 $HMEAN
# run_gaussian_linear "$CONFIG_TEMPLATE" 192 0150 Ch 63 $HMEAN

# run_gaussian_linear "$CONFIG_TEMPLATE" 016 1800 Ah 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 048 0600 Ah 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 064 0450 Ah 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 096 0300 Ah 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 128 0225 Ah 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 192 0150 Ah 63 $HMEAN
# run_gaussian_linear "$CONFIG_TEMPLATE" 192 0150 Ah 63 $HMEAN
