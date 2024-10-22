#!/bin/bash

source $(dirname "$0")"/gen_config.sh"

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
            rotation_axis          = [1.0/sqrt(2.0),0.0,1.0/sqrt(2.0)]}\n
},\n
\n
dt = 1.0*%%%DT*sec,\n
simulation_time = 10.0*day,\n
tau_write = 1.0*day,\n
tau_diagnostics = 1.0*hour,\n
\n
timescheme_name = 'rk4'\n
\n
testcase_section = 'Williamson_test2_config'\n
Williamson_test2_config = {\n
    testcase_name ='Williamson_test2'\n
}\n
"

run_ts2(){
    CONFIG_TEMPLATE=$1
    N=$2
    DT=$3
    STAGGERING=$4
    SBP_ORDER=$5
    UPSTREAM_ORDER=$6

    SUFFIX="N"$N"_dt"$DT"_"$STAGGERING"_"$SBP_ORDER
    OUTFILE="swm_"$SUFFIX".out"
    ERRFILE="errors_"$SUFFIX".txt"
    HDAT="h_"$SUFFIX".dat"

    echo $SUFFIX
    echo $OUTFILE

	gen_config "$CONFIG_TEMPLATE" $N $DT $STAGGERING $SBP_ORDER $UPSTREAM_ORDER > swm_WT2.conf
	mpirun -n $Nprocs $EXE swm_WT2.conf &> $OUTFILE
    grep "l2_h" $OUTFILE > $ERRFILE
    mv h.dat $HDAT
}
run_ts2 "$CONFIG_TEMPLATE" 020 800 Ch 21 up4
run_ts2 "$CONFIG_TEMPLATE" 040 400 Ch 21 up4
run_ts2 "$CONFIG_TEMPLATE" 080 200 Ch 21 up4
run_ts2 "$CONFIG_TEMPLATE" 160 100 Ch 21 up4
run_ts2 "$CONFIG_TEMPLATE" 020 800 Ch 42 up4
run_ts2 "$CONFIG_TEMPLATE" 040 400 Ch 42 up4
run_ts2 "$CONFIG_TEMPLATE" 080 200 Ch 42 up4
run_ts2 "$CONFIG_TEMPLATE" 160 100 Ch 42 up4
run_ts2 "$CONFIG_TEMPLATE" 020 800 Ch 63 up5
run_ts2 "$CONFIG_TEMPLATE" 040 400 Ch 63 up5
run_ts2 "$CONFIG_TEMPLATE" 080 200 Ch 63 up5
run_ts2 "$CONFIG_TEMPLATE" 160 100 Ch 63 up5