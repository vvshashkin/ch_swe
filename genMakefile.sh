#!/bin/bash
LMKL="mkl_lapack95_lp64" #" -lmkl_intel_lp64 -lmkl_core -lmkl_gf_lp64 -lmkl_sequential -lmkl_lapack95_lp64"

CFLAGS=" -DIFORT -c -traceback -init=snan -init=arrays -check all -ftrapuv -fpp -fp-model strict -qmkl"
LFLAGS=" -traceback -init=snan -init=arrays -check all -ftrapuv -fpp -fp-model strict -qmkl "

CFLAGS_OPT=" -DIFORT -fpp -c -traceback -O3 -qmkl"
LFLAGS_OPT=" -fpp -traceback -O3 -qmkl "

FoBiS.py build -s ./src -m Makefile -compiler intel -fc "mpiifort" -lflags "$LFLAGS" -cflags "$CFLAGS" -ext_libs "$LMKL"

FoBiS.py build -s ./src -m Makefile.opt -compiler intel -fc "mpiifort" -lflags "$LFLAGS_OPT" -cflags "$CFLAGS_OPT"  -ext_libs "$LMKL" --obj_dir ./obj_opt --mod_dir ./mod_opt

#build all fortran programs from src
#please note that this is a dirty trick, we need a better solution
#maybe to dig in FoBiS.py options/code
echo "all: \$(addprefix \$(DEXE),\$(EXES))" >> Makefile
