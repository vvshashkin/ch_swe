#!/bin/bash

TESTNAME="SWM Williamson test #2"
EXE=./$1/TS2_MAIN

echo $TESTNAME $EXE

cp ./$1"/src/models/shallow_water/test/ts2/TS2_auto.cfg" swm_model.cfg

mpirun -n 6 $EXE

[ $? -eq 0 ] && echo $TESTNAME passed || echo $TESTNAME failed
