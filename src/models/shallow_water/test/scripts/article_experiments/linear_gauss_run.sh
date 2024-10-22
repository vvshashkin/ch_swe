#!/bin/bash

source $(dirname "$0")"/gen_config_linear.sh"

EXE=$1/SHALLOW_WATER_MAIN
Nprocs=$2
corner=$3

OMEGA='0.0'
[[ $corner == "fcorner" || $corner == "fcorner1" ]] && OMEGA='1e-4'

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
testcase_section = 'Gauss_linear_config'\n"


if [[ $corner == "corner" || $corner == "fcorner" ]]; then
    corner="_"$corner
    CONFIG_TEMPLATE=$CONFIG_TEMPLATE"
    Gauss_linear_config = {\n
        testcase_name ='Gaussian_hill_linear',\n
        lam0 = 0.25*pi,
        phi0 = asin(1.0/sqrt(3.0)),
        refsol_path='$1/ref_solutions/gauss_linear/gauss_h${corner}_cub192_xy.dat'
    }\n
    "
elif [[ $corner == "corner1" || $corner == "fcorner1" ]]; then
    corner="_"$corner
    CONFIG_TEMPLATE=$CONFIG_TEMPLATE"
    Gauss_linear_config = {\n
        testcase_name ='Gaussian_hill_linear',\n
        lam0 = 0.27*pi,
        phi0 = asin(1.0/sqrt(3.0))-0.02*pi,
        refsol_path='$1/ref_solutions/gauss_linear/gauss_h${corner}_cub192_xy.dat'
    }\n
    "
else
    corner=""
    CONFIG_TEMPLATE=$CONFIG_TEMPLATE"
    Gauss_linear_config = {\n
        testcase_name ='Gaussian_hill_linear',\n
        refsol_path='$1/ref_solutions/gauss_linear/gauss_h_cub192_xy.dat'
    }\n
    "
fi
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
    mv h.dat $HDAT
}

# run_gaussian_linear "$CONFIG_TEMPLATE" 016 1800 Ch 21 $HMEAN
# run_gaussian_linear "$CONFIG_TEMPLATE" 032 0900 Ch 21 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 048 0600 Ch 21 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 064 0450 Ch 21 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 096 0300 Ch 21 $HMEAN
# run_gaussian_linear "$CONFIG_TEMPLATE" 128 0225 Ch 21 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 192 0150 Ch 21 $HMEAN

# run_gaussian_linear "$CONFIG_TEMPLATE" 016 1800 Ch 42 $HMEAN
# run_gaussian_linear "$CONFIG_TEMPLATE" 032 0900 Ch 42 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 048 0600 Ch 42 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 064 0450 Ch 42 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 096 0300 Ch 42 $HMEAN
# run_gaussian_linear "$CONFIG_TEMPLATE" 128 0225 Ch 42 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 192 0150 Ch 42 $HMEAN

# run_gaussian_linear "$CONFIG_TEMPLATE" 016 1800 Ch 63 $HMEAN
# run_gaussian_linear "$CONFIG_TEMPLATE" 032 0900 Ch 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 048 0600 Ch 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 064 0450 Ch 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 096 0300 Ch 63 $HMEAN
# run_gaussian_linear "$CONFIG_TEMPLATE" 128 0225 Ch 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 192 0150 Ch 63 $HMEAN

# run_gaussian_linear "$CONFIG_TEMPLATE" 016 1800 Ah 63 $HMEAN
# run_gaussian_linear "$CONFIG_TEMPLATE" 032 0900 Ah 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 048 0600 Ah 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 064 0450 Ah 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 096 0300 Ah 63 $HMEAN
# run_gaussian_linear "$CONFIG_TEMPLATE" 128 0225 Ah 63 $HMEAN
run_gaussian_linear "$CONFIG_TEMPLATE" 192 0150 Ah 63 $HMEAN