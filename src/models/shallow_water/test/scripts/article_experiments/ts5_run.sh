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
            omega                   = Earth_omega,}\n
},\n
\n
dt = 1.0*%%%DT*sec,\n
simulation_time = 15.0*day,\n
tau_write = 1.0*day,\n
tau_diagnostics = 1.0*hour,\n
\n
timescheme_name = 'rk4'\n
\n
testcase_section = 'Williamson_test5_config'\n
Williamson_test5_config = {\n
    testcase_name ='Williamson_test5'\n
}\n
"

run_ts5(){
    CONFIG_TEMPLATE=$1
    N=$2
    DT=$3
    STAGGERING=$4
    SBP_ORDER=$5
    UPSTREAM_ORDER=$6

    SUFFIX="N"$N"_dt"$DT"_"$STAGGERING"_"$SBP_ORDER
    OUTFILE="swm_"$SUFFIX".out"
    HDAT="h_"$SUFFIX".dat"
    ZDAT="curl_"$SUFFIX".dat"

	gen_config "$CONFIG_TEMPLATE" $N $DT $STAGGERING $SBP_ORDER $UPSTREAM_ORDER > swm_WT2.conf
	mpirun -n $Nprocs $EXE swm_WT2.conf &> $OUTFILE
    mv h.dat $HDAT
    mv curl.dat $ZDAT
}
run_ts5 "$CONFIG_TEMPLATE" 020 800 Ch 21 up4
run_ts5 "$CONFIG_TEMPLATE" 040 400 Ch 21 up4
run_ts5 "$CONFIG_TEMPLATE" 080 200 Ch 21 up4
run_ts5 "$CONFIG_TEMPLATE" 160 100 Ch 21 up4
run_ts5 "$CONFIG_TEMPLATE" 020 800 Ch 42 up4
run_ts5 "$CONFIG_TEMPLATE" 040 400 Ch 42 up4
run_ts5 "$CONFIG_TEMPLATE" 080 200 Ch 42 up4
run_ts5 "$CONFIG_TEMPLATE" 160 100 Ch 42 up4
run_ts5 "$CONFIG_TEMPLATE" 020 800 Ch 63 up5
run_ts5 "$CONFIG_TEMPLATE" 040 400 Ch 63 up5
run_ts5 "$CONFIG_TEMPLATE" 080 200 Ch 63 up5
run_ts5 "$CONFIG_TEMPLATE" 160 100 Ch 63 up5