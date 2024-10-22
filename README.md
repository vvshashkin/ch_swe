# Cubed sphere shallow water model

This code is used in *Summation-by-Parts Finite-Difference Method for Linear Shallow Water Equations on Staggered Curvilinear Grids in Closed Domains* article.

## Build
Intel fortran compiler (ifort) is required
```bash
make -f Makefile.opt SHALLOW_WATER_MAIN
```

## Running experiments

First, download reference solutions from https://disk.yandex.ru/d/_PjySTdCnQAIqA and place it ref_solutions directory in the model repository/

### Linear gaussian test

Run the following commands
```bash
mkdir gauss_linear; cd gauss_linear
path_to_model_repo/src/models/shallow_water/test/scripts/article_experiments/linear_gauss_run.sh path_to_exe n_mpi exp_variant
#visualization
python3 path_to_model_repo/src/models/shallow_water/test/scripts/article_experiments/plot_gauss_linear.py path_to_model_repo/ref_solutions exp_variant
```
where
- ``path_to_model_repo`` is the path to the model repository folder (e.q. ``~/repos/ch_swe/``)
- ``path_to_exe`` is the path to directory containing compiled ``SHALLOW_WATER_MAIN``
- ``n_mpi`` is the number of mpi processes to run  ``SHALLOW_WATER_MAIN``
- ``exp_variant`` is one of (empty), "corner", "fcorner"

The NCAR graphics python package https://pyngl.ucar.edu is required to draw plots.

### Solid rotation

```bash
mkdir gauss_linear; cd gauss_linear
path_to_model_repo/src/models/shallow_water/test/scripts/article_experiments/run_ts2_linear.sh path_to_exe n_mpi
#visualization
python3 path_to_model_repo/src/models/shallow_water/test/scripts/article_experiments/plot_gauss_linear.py .
```
