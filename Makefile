#!/usr/bin/make

#main building variables
DSRC    = src
DOBJ    = obj/
DMOD    = mod/
DEXE    = ./
LIBS    =  -lmkl_lapack95_lp64 
FC      = mpiifort
OPTSC   =  -DIFORT -c -traceback -init=snan -init=arrays -check all -ftrapuv -fpp -fp-model strict -qmkl -module mod
OPTSL   =  -traceback -init=snan -init=arrays -check all -ftrapuv -fpp -fp-model strict -qmkl  -module mod
VPATH   = $(DSRC) $(DOBJ) $(DMOD)
MKDIRS  = $(DOBJ) $(DMOD) $(DEXE)
LCEXES  = $(shell echo $(EXES) | tr '[:upper:]' '[:lower:]')
EXESPO  = $(addsuffix .o,$(LCEXES))
EXESOBJ = $(addprefix $(DOBJ),$(EXESPO))

#auxiliary variables
COTEXT  = "Compiling $(<F)"
LITEXT  = "Assembling $@"

#building rules
$(DEXE)TEST_STANDARD_ATMOSPHERE_MAIN: $(MKDIRS) $(DOBJ)test_standard_atmosphere_main.o
	@rm -f $(filter-out $(DOBJ)test_standard_atmosphere_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TEST_STANDARD_ATMOSPHERE_MAIN
$(DEXE)SHALLOW_WATER_MAIN: $(MKDIRS) $(DOBJ)shallow_water_main.o
	@rm -f $(filter-out $(DOBJ)shallow_water_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) SHALLOW_WATER_MAIN
$(DEXE)BAROTROPIC_INST_MAIN: $(MKDIRS) $(DOBJ)barotropic_inst_main.o
	@rm -f $(filter-out $(DOBJ)barotropic_inst_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) BAROTROPIC_INST_MAIN
$(DEXE)ELDRED_TEST_MAIN: $(MKDIRS) $(DOBJ)eldred_test_main.o
	@rm -f $(filter-out $(DOBJ)eldred_test_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) ELDRED_TEST_MAIN
$(DEXE)FORCING_TEST_MAIN: $(MKDIRS) $(DOBJ)forcing_test_main.o
	@rm -f $(filter-out $(DOBJ)forcing_test_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) FORCING_TEST_MAIN
$(DEXE)TS5_MAIN: $(MKDIRS) $(DOBJ)ts5_main.o
	@rm -f $(filter-out $(DOBJ)ts5_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TS5_MAIN
$(DEXE)GAUSSIAN_HILL_MAIN: $(MKDIRS) $(DOBJ)gaussian_hill_main.o
	@rm -f $(filter-out $(DOBJ)gaussian_hill_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) GAUSSIAN_HILL_MAIN
$(DEXE)TS2_MAIN: $(MKDIRS) $(DOBJ)ts2_main.o
	@rm -f $(filter-out $(DOBJ)ts2_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) TS2_MAIN
$(DEXE)RH4_WAVE_MAIN: $(MKDIRS) $(DOBJ)rh4_wave_main.o
	@rm -f $(filter-out $(DOBJ)rh4_wave_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) RH4_WAVE_MAIN
$(DEXE)GET_LATLON_MAIN: $(MKDIRS) $(DOBJ)get_latlon_main.o
	@rm -f $(filter-out $(DOBJ)get_latlon_main.o,$(EXESOBJ))
	@echo $(LITEXT)
	@$(FC) $(OPTSL) $(DOBJ)*.o $(LIBS) -o $@
EXES := $(EXES) GET_LATLON_MAIN

#compiling rules
$(DOBJ)operator_mod.o: src/operator_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)key_value_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)tile_mod.o: src/tile_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ideal_gas_law_mod.o: src/ideal_gas_law_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)namelist_read_mod.o: src/namelist_read_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)vec_math_mod.o: src/vec_math_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)grid_field_based_vector_mod.o: src/grid_field_based_vector_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)vector_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_quadrature_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)const_mod.o: src/const_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)grid_field_factory_mod.o: src/grid_field_factory_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tiles_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)grid_field_collection_mod.o: src/grid_field_collection_mod.f90 \
	$(DOBJ)linked_list_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)string_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)cmd_args_mod.o: src/cmd_args_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sph_coords_mod.o: src/sph_coords_mod.f90 \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)config_mod.o: src/config_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)global_diag_mod.o: src/global_diag_mod.f90 \
	$(DOBJ)container_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)stvec_abstract_mod.o: src/stvec_abstract_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)tiles_mod.o: src/tiles_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)grid_field_mod.o: src/grid_field_mod.f90 \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)outputer_abstract_mod.o: src/outputer/outputer_abstract_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)outputer_factory_mod.o: src/outputer/outputer_factory_mod.f90 \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)master_paneled_outputer_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)latlon_outputer_mod.o \
	$(DOBJ)regrid_factory_mod.o \
	$(DOBJ)latlon_grid_generator_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)mpi_paneled_outputer_mod.o: src/outputer/mpi_paneled_outputer_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)master_paneled_outputer_mod.o: src/outputer/master_paneled_outputer_mod.f90 \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)tiles_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)latlon_outputer_mod.o: src/outputer/latlon_outputer_mod.f90 \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)tiles_mod.o \
	$(DOBJ)abstract_regridder_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_quadrature_mod.o: src/quadrature/sbp_quadrature_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_quadrature_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)quadrature_factory_mod.o: src/quadrature/quadrature_factory_mod.f90 \
	$(DOBJ)abstract_quadrature_mod.o \
	$(DOBJ)default_quadrature_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sbp_quadrature_mod.o \
	$(DOBJ)sbp_operators_collection_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_quadrature_mod.o: src/quadrature/abstract_quadrature_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)default_quadrature_mod.o: src/quadrature/default_quadrature_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_quadrature_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)diagonal_preconditioner_mod.o: src/solver/diagonal_preconditioner_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)preconditioner_mod.o \
	$(DOBJ)vector_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)iterative_solver_factory_mod.o: src/solver/iterative_solver_factory_mod.f90 \
	$(DOBJ)iterative_solver_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)cg_solver_mod.o \
	$(DOBJ)bicgstab_solver_mod.o \
	$(DOBJ)chebyshev_solver_mod.o \
	$(DOBJ)richardson_solver_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_based_vector_mod.o \
	$(DOBJ)abstract_quadrature_mod.o \
	$(DOBJ)quadrature_factory_mod.o \
	$(DOBJ)prolongation_factory_mod.o \
	$(DOBJ)restriction_factory_mod.o \
	$(DOBJ)linear_operator_factory_mod.o \
	$(DOBJ)preconditioner_factory_mod.o \
	$(DOBJ)multigrid_solver_mod.o \
	$(DOBJ)geosci_config_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)linear_operator_factory_mod.o: src/solver/linear_operator_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)linear_operator_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)swm_helm_operator_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)grad_factory_mod.o \
	$(DOBJ)co2contra_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)iterative_solver_mod.o: src/solver/iterative_solver_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)linear_operator_mod.o \
	$(DOBJ)vector_mod.o \
	$(DOBJ)preconditioner_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)bicgstab_solver_mod.o: src/solver/bicgstab_solver_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)iterative_solver_mod.o \
	$(DOBJ)linear_operator_mod.o \
	$(DOBJ)vector_mod.o \
	$(DOBJ)preconditioner_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)preconditioner_mod.o: src/solver/preconditioner_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)vector_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)chebyshev_solver_mod.o: src/solver/chebyshev_solver_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)iterative_solver_mod.o \
	$(DOBJ)linear_operator_mod.o \
	$(DOBJ)vector_mod.o \
	$(DOBJ)preconditioner_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)richardson_solver_mod.o: src/solver/richardson_solver_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)iterative_solver_mod.o \
	$(DOBJ)linear_operator_mod.o \
	$(DOBJ)vector_mod.o \
	$(DOBJ)preconditioner_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)linear_operator_mod.o: src/solver/linear_operator_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)vector_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)preconditioner_factory_mod.o: src/solver/preconditioner_factory_mod.f90 \
	$(DOBJ)preconditioner_mod.o \
	$(DOBJ)grid_field_based_vector_mod.o \
	$(DOBJ)linear_operator_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_quadrature_mod.o \
	$(DOBJ)diagonal_preconditioner_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)vector_mod.o: src/solver/vector_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)cg_solver_mod.o: src/solver/cg_solver_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)iterative_solver_mod.o \
	$(DOBJ)linear_operator_mod.o \
	$(DOBJ)vector_mod.o \
	$(DOBJ)preconditioner_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_prolongation_mod.o: src/solver/multigrid/abstract_prolongation_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)vector_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)prolongation_factory_mod.o: src/solver/multigrid/prolongation_factory_mod.f90 \
	$(DOBJ)abstract_prolongation_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)prolongation_bilinear_const_mod.o \
	$(DOBJ)exchange_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)restriction_factory_mod.o: src/solver/multigrid/restriction_factory_mod.f90 \
	$(DOBJ)abstract_restriction_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)restriction_4point_cell_average_mod.o \
	$(DOBJ)exchange_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_restriction_mod.o: src/solver/multigrid/abstract_restriction_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)vector_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)restriction_4point_cell_average_mod.o: src/solver/multigrid/restriction_4point_cell_average_mod.f90 \
	$(DOBJ)abstract_restriction_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)vector_mod.o \
	$(DOBJ)grid_field_based_vector_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)multigrid_solver_mod.o: src/solver/multigrid/multigrid_solver_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)iterative_solver_mod.o \
	$(DOBJ)linear_operator_mod.o \
	$(DOBJ)vector_mod.o \
	$(DOBJ)preconditioner_mod.o \
	$(DOBJ)abstract_restriction_mod.o \
	$(DOBJ)abstract_prolongation_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)prolongation_bilinear_const_mod.o: src/solver/multigrid/prolongation_bilinear_const_mod.f90 \
	$(DOBJ)abstract_prolongation_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)vector_mod.o \
	$(DOBJ)grid_field_based_vector_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)mesh_mod.o: src/mesh/mesh_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)mesh_collection_mod.o: src/mesh/mesh_collection_mod.f90 \
	$(DOBJ)linked_list_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)mesh_factory_mod.o: src/mesh/mesh_factory_mod.f90 \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tiles_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)orography_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)metric_mod.o \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)topology_factory_mod.o: src/topology/topology_factory_mod.f90 \
	$(DOBJ)topology_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)cubed_sphere_topology_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)cubed_sphere_topology_mod.o: src/topology/cubed_sphere_topology_mod.f90 \
	$(DOBJ)topology_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)topology_mod.o: src/topology/topology_mod.f90 \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ecs_halo_vec_a_factory_mod.o: src/equiang_cs/ecs_halo_vec_a_factory_mod.f90 \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)ecs_halo_vec_a_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)ecs_halo_factory_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)metric_2d_ecs_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ecs_halo_vec_xy_mod.o: src/equiang_cs/ecs_halo_vec_xy_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)ecs_halo_xy_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ecs_halo_mod.o: src/equiang_cs/ecs_halo_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ecs_halo_vec_c_mod.o: src/equiang_cs/ecs_halo_vec_c_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ecs_halo_vec_c_factory_mod.o: src/equiang_cs/ecs_halo_vec_c_factory_mod.f90 \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)ecs_halo_vec_c_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)metric_2d_ecs_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ecs_halo_vec_xy_factory_mod.o: src/equiang_cs/ecs_halo_vec_xy_factory_mod.f90 \
	$(DOBJ)ecs_halo_vec_xy_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)ecs_halo_factory_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)metric_2d_ecs_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ecs_halo_xy_mod.o: src/equiang_cs/ecs_halo_xy_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ecs_halo_vec_ch_factory_mod.o: src/equiang_cs/ecs_halo_vec_Ch_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)ecs_halo_vec_ch_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)particles_mod.o \
	$(DOBJ)distribute_particles_mod.o \
	$(DOBJ)geosci_config_mod.o \
	$(DOBJ)string_mod.o \
	$(DOBJ)particle_interp_factory_mod.o \
	$(DOBJ)metric_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ecs_ah_vec_sync_factory_mod.o: src/equiang_cs/ecs_Ah_vec_sync_factory_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)ecs_ah_vec_sync_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)metric_2d_ecs_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ecs_ah_vec_sync_mod.o: src/equiang_cs/ecs_Ah_vec_sync_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ecs_halo_vec_ch_mod.o: src/equiang_cs/ecs_halo_vec_Ch_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)abstract_particle_interpolation_mod.o \
	$(DOBJ)particle_values_mod.o \
	$(DOBJ)array_tools_mod.o \
	$(DOBJ)assembly_particles_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ecs_halo_vec_a_mod.o: src/equiang_cs/ecs_halo_vec_a_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ecs_halo_factory_mod.o: src/equiang_cs/ecs_halo_factory_mod.f90 \
	$(DOBJ)ecs_halo_mod.o \
	$(DOBJ)ecs_halo_xy_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)halo_factory_mod.o: src/halo/halo_factory_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)ecs_halo_factory_mod.o \
	$(DOBJ)ecs_halo_vec_a_factory_mod.o \
	$(DOBJ)ecs_halo_vec_xy_factory_mod.o \
	$(DOBJ)ecs_halo_vec_c_factory_mod.o \
	$(DOBJ)ecs_ah_vec_sync_factory_mod.o \
	$(DOBJ)ecs_halo_vec_ch_factory_mod.o \
	$(DOBJ)halo_default_exchange_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)halo_ah_scalar_sync_mod.o \
	$(DOBJ)halo_c_default_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)halo_ah_scalar_sync_mod.o: src/halo/halo_Ah_scalar_sync_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)halo_mod.o: src/halo/halo_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)halo_default_exchange_mod.o: src/halo/halo_default_exchange_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)halo_c_default_mod.o: src/halo/halo_C_default_mod.f90 \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)crank_nikolson_mod.o: src/time_schemes/crank_nikolson_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)imex_m2_mod.o: src/time_schemes/imex_M2_mod.f90 \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)timescheme_factory_mod.o: src/time_schemes/timescheme_factory_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)explicit_eul1_mod.o \
	$(DOBJ)implicit_eul_mod.o \
	$(DOBJ)rk4_mod.o \
	$(DOBJ)ars343.o \
	$(DOBJ)imex_m2_mod.o \
	$(DOBJ)crank_nikolson_mod.o \
	$(DOBJ)lsrk_mod.o \
	$(DOBJ)sisl_settls_mod.o \
	$(DOBJ)sl_ukmo_mod.o \
	$(DOBJ)abstract_generic_config_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)lsrk_mod.o: src/time_schemes/lsrk_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)explicit_eul1_mod.o: src/time_schemes/explicit_Eul1_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sisl_settls_mod.o: src/time_schemes/SISL_SETTLS_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)rk4_mod.o: src/time_schemes/rk4_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)implicit_eul_mod.o: src/time_schemes/implicit_Eul_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)timescheme_mod.o: src/time_schemes/timescheme_mod.f90 \
	$(DOBJ)operator_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sl_ukmo_mod.o: src/time_schemes/SL_UKMO_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)timer_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ars343.o: src/time_schemes/ars343.f90 \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)cubsphere_scalar_field_mod.o: src/idealized_fields/cubsphere_scalar_field_mod.f90 \
	$(DOBJ)scalar_field_native_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)simple_latlon_scalar_field_mod.o: src/idealized_fields/simple_latlon_scalar_field_mod.f90 \
	$(DOBJ)scalar_field_latlon_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)basic_idealized_fields_mod.o: src/idealized_fields/basic_idealized_fields_mod.f90 \
	$(DOBJ)idealized_field_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)baroclinic_instability_jw2006_mod.o: src/idealized_fields/baroclinic_instability_JW2006_mod.f90 \
	$(DOBJ)vector_field_latlon_mod.o \
	$(DOBJ)scalar_field_latlon_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)xyz_scalar_field_mod.o: src/idealized_fields/xyz_scalar_field_mod.f90 \
	$(DOBJ)scalar_field_cartesian_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)solid_rotation_fields_mod.o: src/idealized_fields/solid_rotation_fields_mod.f90 \
	$(DOBJ)vector_field_latlon_mod.o \
	$(DOBJ)vector_field_cartesian_mod.o \
	$(DOBJ)scalar_field_cartesian_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)idealized_fields_factory_mod.o: src/idealized_fields/idealized_fields_factory_mod.f90 \
	$(DOBJ)idealized_field_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)string_mod.o \
	$(DOBJ)solid_rotation_fields_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)baroclinic_instability_jw2006_mod.o \
	$(DOBJ)gw_testcase_fields_mod.o \
	$(DOBJ)basic_idealized_fields_mod.o \
	$(DOBJ)vertical_profiles_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)vertical_profiles_mod.o: src/idealized_fields/vertical_profiles_mod.f90 \
	$(DOBJ)scalar_field_cartesian_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)gw_testcase_fields_mod.o: src/idealized_fields/GW_testcase_fields_mod.f90 \
	$(DOBJ)scalar_field_cartesian_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)solid_rotation_fields_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)scalar_field_cartesian_mod.o: src/idealized_fields/abstract_classes/scalar_field_cartesian_mod.f90 \
	$(DOBJ)idealized_field_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)vector_idealized_field_mod.o: src/idealized_fields/abstract_classes/vector_idealized_field_mod.f90 \
	$(DOBJ)idealized_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)vector_field_latlon_mod.o: src/idealized_fields/abstract_classes/vector_field_latlon_mod.f90 \
	$(DOBJ)vector_idealized_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)vector_field_cartesian_mod.o: src/idealized_fields/abstract_classes/vector_field_cartesian_mod.f90 \
	$(DOBJ)vector_idealized_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)scalar_field_native_mod.o: src/idealized_fields/abstract_classes/scalar_field_native_mod.f90 \
	$(DOBJ)idealized_field_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)idealized_field_mod.o: src/idealized_fields/abstract_classes/idealized_field_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)scalar_field_latlon_mod.o: src/idealized_fields/abstract_classes/scalar_field_latlon_mod.f90 \
	$(DOBJ)idealized_field_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)latlon_grid_generator_mod.o: src/stuff/latlon_grid_generator_mod.f90 \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)key_value_mod.o: src/stuff/key_value_mod.f90 \
	$(DOBJ)string_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)array_tools_mod.o: src/stuff/array_tools_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)basic_collection_mod.o: src/stuff/basic_collection_mod.f90 \
	$(DOBJ)linked_list_mod.o \
	$(DOBJ)string_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)smart_array_mod.o: src/stuff/smart_array_mod.f90 \
	$(DOBJ)array_tools_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)error_mod.o: src/stuff/error_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)timer_mod.o: src/stuff/timer_mod.f90 \
	$(DOBJ)str_util_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)linked_list_mod.o: src/stuff/linked_list_mod.f90 \
	$(DOBJ)string_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)time_util_mod.o: src/stuff/time_util_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)conv_order_mod.o: src/stuff/conv_order_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)string_mod.o: src/stuff/string_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_scorer_mod.o: src/stuff/abstract_scorer_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)str_util_mod.o: src/stuff/str_util_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_generic_config_mod.o: src/stuff/config/abstract_generic_config_mod.f90 src/stuff/config/generic_config_get_interface.fpp src/stuff/config/generic_config_get_array_interface.fpp \
	$(DOBJ)string_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)config_tools_mod.o: src/stuff/config/config_tools_mod.f90 \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)geosci_config_mod.o \
	$(DOBJ)namelist_read_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)geosci_config_mod.o: src/stuff/config/geosci_config/geosci_config_mod.f90 src/stuff/config/geosci_config/geosci_config_get_var.fpp src/stuff/config/geosci_config/geosci_config_get_array.fpp \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)string_mod.o \
	$(DOBJ)str_util_mod.o \
	$(DOBJ)linked_list_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)geosci_config_parser_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)geosci_config_parser_mod.o: src/stuff/config/geosci_config/geosci_config_parser_mod.f90 \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)string_mod.o \
	$(DOBJ)str_util_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)arithmetic_parser_mod.o \
	$(DOBJ)arithmetic_parser_val_mod.o \
	$(DOBJ)error_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)arithmetic_parser_tools_mod.o: src/stuff/config/geosci_config/arithmetic_parser/arithmetic_parser_tools_mod.f90 \
	$(DOBJ)arithmetic_parser_val_mod.o \
	$(DOBJ)error_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)arithmetic_parser_mod.o: src/stuff/config/geosci_config/arithmetic_parser/arithmetic_parser_mod.f90 \
	$(DOBJ)arithmetic_parser_tools_mod.o \
	$(DOBJ)arithmetic_parser_val_mod.o \
	$(DOBJ)str_util_mod.o \
	$(DOBJ)error_mod.o \
	$(DOBJ)abstract_generic_config_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)arithmetic_parser_val_mod.o: src/stuff/config/geosci_config/arithmetic_parser/arithmetic_parser_val_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)standard_atmosphere_mod.o: src/stuff/standard_atmosphere/standard_atmosphere_mod.f90 \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)test_standard_atmosphere_main.o: src/stuff/standard_atmosphere/test_standard_atmosphere_main.f90 \
	$(DOBJ)standard_atmosphere_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)solve_tridiag_mod.o: src/stuff/linalg/solve_tridiag_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_trajectory_solver_mod.o: src/SL/trajectory_solver/abstract_trajectory_solver_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_flexible_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)shallow_atm_trajectory_tools_mod.o: src/SL/trajectory_solver/shallow_atm_trajectory_tools_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)metric_mod.o \
	$(DOBJ)shallow_atm_metric_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)shallow_spherical_trajectory_solver_mod.o: src/SL/trajectory_solver/shallow_spherical_trajectory_solver_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_flexible_mod.o \
	$(DOBJ)abstract_trajectory_solver_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)shallow_atm_trajectory_tools_mod.o \
	$(DOBJ)abstract_dep_points_interp_driver_mod.o \
	$(DOBJ)timer_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)trajectory_solver_factory_mod.o: src/SL/trajectory_solver/trajectory_solver_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_trajectory_solver_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)geosci_config_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)dep_points_interp_driver_factory_mod.o \
	$(DOBJ)stvec_flexible_factory_mod.o \
	$(DOBJ)string_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)shallow_spherical_trajectory_solver_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)dep_points_interp_driver_factory_mod.o: src/SL/dep_points_interp_driver/dep_points_interp_driver_factory_mod.f90 \
	$(DOBJ)abstract_dep_points_interp_driver_mod.o \
	$(DOBJ)dp_interpolator_factory_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)wd_interp_driver_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)tiles_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)string_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)abstract_generic_config_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)wd_interp_driver_mod.o: src/SL/dep_points_interp_driver/WD_interp_driver_mod.f90 \
	$(DOBJ)abstract_dep_points_interp_driver_mod.o \
	$(DOBJ)abstract_dp_interpolator_mod.o \
	$(DOBJ)stvec_flexible_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)string_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_dep_points_interp_driver_mod.o: src/SL/dep_points_interp_driver/abstract_dep_points_interp_driver_mod.f90 \
	$(DOBJ)stvec_flexible_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)dp_interpolators_mod.o: src/SL/dp_interpolators/dp_interpolators_mod.f90 \
	$(DOBJ)abstract_dp_interpolator_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_dp_interpolator_mod.o: src/SL/dp_interpolators/abstract_dp_interpolator_mod.f90 \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)dp_interpolator_factory_mod.o: src/SL/dp_interpolators/dp_interpolator_factory_mod.f90 \
	$(DOBJ)abstract_dp_interpolator_mod.o \
	$(DOBJ)dp_interpolators_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)shallow_water_main.o: src/models/shallow_water/shallow_water_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)cmd_args_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)config_tools_mod.o \
	$(DOBJ)shallow_water_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)shallow_water_mod.o: src/models/shallow_water/shallow_water_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)stvec_swm_factory_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)operator_swm_factory_mod.o \
	$(DOBJ)swm_output_diag_mod.o \
	$(DOBJ)swm_output_diag_factory_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)operator_swm_diff_mod.o \
	$(DOBJ)operator_swm_diff_factory_mod.o \
	$(DOBJ)swm_forcing_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)abstract_scorer_mod.o \
	$(DOBJ)timer_mod.o \
	$(DOBJ)williamson_test2_mod.o \
	$(DOBJ)williamson_test5_mod.o \
	$(DOBJ)rh4_testcase_mod.o \
	$(DOBJ)barotropic_instability_testcase_mod.o \
	$(DOBJ)eldred_testcase_mod.o \
	$(DOBJ)gaussian_hill_linear_test_mod.o \
	$(DOBJ)oscillating_gaussian_hill_linear_test.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)swm_forcing_mod.o: src/models/shallow_water/swm_forcing_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)swm_helm_operator_mod.o: src/models/shallow_water/operator/swm_helm_operator_mod.f90 \
	$(DOBJ)linear_operator_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)vector_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_based_vector_mod.o \
	$(DOBJ)abstract_massflux_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)operator_swm_lin_mod.o: src/models/shallow_water/operator/operator_swm_lin_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)iterative_solver_mod.o \
	$(DOBJ)linear_operator_mod.o \
	$(DOBJ)vector_mod.o \
	$(DOBJ)abstract_quadrature_mod.o \
	$(DOBJ)swm_helm_operator_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)grid_field_based_vector_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)coriolis_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)operator_swm_mod.o: src/models/shallow_water/operator/operator_swm_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)abstract_ke_mod.o \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)abstract_quadrature_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)coriolis_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)operator_swm_factory_mod.o: src/models/shallow_water/operator/operator_swm_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)operator_swm_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)grad_factory_mod.o \
	$(DOBJ)curl_factory_mod.o \
	$(DOBJ)coriolis_factory_mod.o \
	$(DOBJ)ke_factory_mod.o \
	$(DOBJ)massflux_factory_mod.o \
	$(DOBJ)co2contra_factory_mod.o \
	$(DOBJ)quadrature_factory_mod.o \
	$(DOBJ)hordiff_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)operator_swm_vec_inv_imex_mod.o \
	$(DOBJ)stvec_swm_factory_mod.o \
	$(DOBJ)cg_solver_mod.o \
	$(DOBJ)bicgstab_solver_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)swm_helm_operator_mod.o \
	$(DOBJ)grid_field_based_vector_mod.o \
	$(DOBJ)iterative_solver_factory_mod.o \
	$(DOBJ)operator_adv_swm_mod.o \
	$(DOBJ)vector_advection_factory_mod.o \
	$(DOBJ)operator_swm_sisl_mod.o \
	$(DOBJ)dep_points_interp_driver_factory_mod.o \
	$(DOBJ)trajectory_solver_factory_mod.o \
	$(DOBJ)string_mod.o \
	$(DOBJ)stvec_flexible_factory_mod.o \
	$(DOBJ)skew_swm_operator_mod.o \
	$(DOBJ)flux_div_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)operator_swm_lin_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)operator_swm_diff_mod.o: src/models/shallow_water/operator/operator_swm_diff_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)abstract_hordiff_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)operator_swm_vec_inv_imex_mod.o: src/models/shallow_water/operator/operator_swm_vec_inv_imex_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)abstract_ke_mod.o \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)abstract_quadrature_mod.o \
	$(DOBJ)iterative_solver_mod.o \
	$(DOBJ)linear_operator_mod.o \
	$(DOBJ)vector_mod.o \
	$(DOBJ)swm_helm_operator_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)grid_field_based_vector_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)coriolis_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)skew_swm_operator_mod.o: src/models/shallow_water/operator/skew_swm_operator_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_flux_div_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)abstract_quadrature_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)shallow_atm_trajectory_tools_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)key_value_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)operator_adv_swm_mod.o: src/models/shallow_water/operator/operator_adv_swm_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)abstract_quadrature_mod.o \
	$(DOBJ)abstract_vector_advection_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)key_value_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)operator_swm_sisl_mod.o: src/models/shallow_water/operator/operator_swm_SISL_mod.f90 \
	$(DOBJ)stvec_flexible_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)abstract_ke_mod.o \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)abstract_quadrature_mod.o \
	$(DOBJ)iterative_solver_mod.o \
	$(DOBJ)linear_operator_mod.o \
	$(DOBJ)vector_mod.o \
	$(DOBJ)swm_helm_operator_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)grid_field_based_vector_mod.o \
	$(DOBJ)abstract_dep_points_interp_driver_mod.o \
	$(DOBJ)abstract_trajectory_solver_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)shallow_atm_trajectory_tools_mod.o \
	$(DOBJ)coriolis_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)operator_swm_diff_factory_mod.o: src/models/shallow_water/operator/operator_swm_diff_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)operator_swm_diff_mod.o \
	$(DOBJ)hordiff_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)swm_output_diag_mod.o: src/models/shallow_water/output_and_diag/swm_output_diag_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_flexible_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)abstract_quadrature_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)swm_output_diag_factory_mod.o: src/models/shallow_water/output_and_diag/swm_output_diag_factory_mod.f90 \
	$(DOBJ)swm_output_diag_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)curl_factory_mod.o \
	$(DOBJ)co2contra_factory_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)quadrature_factory_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)coriolis_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)config_swm_mod.o: src/models/shallow_water/config/config_swm_mod.f90 \
	$(DOBJ)config_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)barotropic_inst_main.o: src/models/shallow_water/test/barotropic_instability/barotropic_inst_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)barotropic_inst_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)barotropic_instability_u_mod.o: src/models/shallow_water/test/barotropic_instability/barotropic_instability_u_mod.f90 \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)barotropic_inst_mod.o: src/models/shallow_water/test/barotropic_instability/barotropic_inst_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)stvec_swm_factory_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)operator_swm_factory_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)config_tools_mod.o \
	$(DOBJ)operator_swm_mod.o \
	$(DOBJ)operator_swm_vec_inv_imex_mod.o \
	$(DOBJ)operator_swm_diff_mod.o \
	$(DOBJ)operator_swm_diff_factory_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)curl_factory_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)co2contra_factory_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)namelist_read_mod.o \
	$(DOBJ)barotropic_instability_u_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)eldred_test_mod.o: src/models/shallow_water/test/Eldred_test/Eldred_test_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)stvec_swm_factory_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)operator_swm_factory_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)config_tools_mod.o \
	$(DOBJ)operator_swm_diff_mod.o \
	$(DOBJ)operator_swm_diff_factory_mod.o \
	$(DOBJ)operator_swm_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)random_friction_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)namelist_read_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)eldred_test_main.o: src/models/shallow_water/test/Eldred_test/Eldred_test_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)eldred_test_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)random_friction_mod.o: src/models/shallow_water/test/forcing_test/random_friction_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)forcing_test_main.o: src/models/shallow_water/test/forcing_test/forcing_test_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)forcing_test_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)forcing_test_mod.o: src/models/shallow_water/test/forcing_test/forcing_test_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)stvec_swm_factory_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)operator_swm_factory_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)config_swm_mod.o \
	$(DOBJ)operator_swm_mod.o \
	$(DOBJ)operator_swm_diff_mod.o \
	$(DOBJ)operator_swm_diff_factory_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)random_friction_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)namelist_read_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ts5_main.o: src/models/shallow_water/test/ts5/ts5_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)ts5_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ts5_mod.o: src/models/shallow_water/test/ts5/ts5_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)stvec_swm_factory_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)operator_swm_mod.o \
	$(DOBJ)operator_swm_factory_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)config_tools_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)operator_swm_diff_mod.o \
	$(DOBJ)operator_swm_diff_factory_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)namelist_read_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)gaussian_hill_main.o: src/models/shallow_water/test/gaussian_hill_lin/gaussian_hill_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)gaussian_hill_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)gaussian_hill_mod.o: src/models/shallow_water/test/gaussian_hill_lin/gaussian_hill_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)stvec_swm_factory_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)operator_swm_factory_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)config_tools_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)namelist_read_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)config_gaussian_hill_mod.o: src/models/shallow_water/test/gaussian_hill_lin/config_gaussian_hill_mod.f90 \
	$(DOBJ)config_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ts2_main.o: src/models/shallow_water/test/ts2/ts2_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)ts2_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ts2_mod.o: src/models/shallow_water/test/ts2/ts2_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)stvec_swm_factory_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)operator_swm_factory_mod.o \
	$(DOBJ)operator_swm_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)config_tools_mod.o \
	$(DOBJ)timer_mod.o \
	$(DOBJ)operator_swm_diff_mod.o \
	$(DOBJ)operator_swm_diff_factory_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)namelist_read_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)rh4_wave_mod.o: src/models/shallow_water/test/RH4_wave/RH4_wave_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)stvec_swm_factory_mod.o \
	$(DOBJ)operator_mod.o \
	$(DOBJ)operator_swm_factory_mod.o \
	$(DOBJ)timescheme_mod.o \
	$(DOBJ)timescheme_factory_mod.o \
	$(DOBJ)outputer_abstract_mod.o \
	$(DOBJ)outputer_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)config_tools_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)operator_swm_diff_mod.o \
	$(DOBJ)operator_swm_diff_factory_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)key_value_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)namelist_read_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)rh4_wave_main.o: src/models/shallow_water/test/RH4_wave/RH4_wave_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)rh4_wave_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)williamson_test2_mod.o: src/models/shallow_water/testcases/Williamson_test2_mod.f90 \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)stvec_flexible_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_scorer_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)oscillating_gaussian_hill_linear_test.o: src/models/shallow_water/testcases/oscillating_gaussian_hill_linear_test.f90 \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)stvec_flexible_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_scorer_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)williamson_test5_mod.o: src/models/shallow_water/testcases/Williamson_test5_mod.f90 \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)stvec_flexible_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)rh4_testcase_mod.o: src/models/shallow_water/testcases/RH4_testcase_mod.f90 \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)stvec_flexible_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)gaussian_hill_linear_test_mod.o: src/models/shallow_water/testcases/gaussian_hill_linear_test_mod.f90 \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)stvec_flexible_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_scorer_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)stvec_swm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)eldred_testcase_mod.o: src/models/shallow_water/testcases/Eldred_testcase_mod.f90 \
	$(DOBJ)swm_forcing_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)stvec_flexible_mod.o \
	$(DOBJ)stvec_swm_factory_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)random_friction_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)barotropic_instability_testcase_mod.o: src/models/shallow_water/testcases/barotropic_instability_testcase_mod.f90 \
	$(DOBJ)test_fields_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)stvec_flexible_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)stvec_swm_mod.o: src/models/shallow_water/stvec/stvec_swm_mod.f90 \
	$(DOBJ)stvec_flexible_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)stvec_swm_factory_mod.o: src/models/shallow_water/stvec/stvec_swm_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_swm_mod.o \
	$(DOBJ)stvec_flexible_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)string_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)operator_iomega_mod.o: src/models/iomega_model/operator_iomega_mod.f90 \
	$(DOBJ)operator_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)stvec_iomega_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)stvec_iomega_mod.o: src/models/iomega_model/stvec_iomega_mod.f90 \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)domain_factory_mod.o: src/domain/domain_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)topology_factory_mod.o \
	$(DOBJ)cubed_sphere_topology_mod.o \
	$(DOBJ)metric_mod.o \
	$(DOBJ)metric_factory_mod.o \
	$(DOBJ)orography_mod.o \
	$(DOBJ)orography_factory_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)geosci_config_mod.o \
	$(DOBJ)mesh_factory_mod.o \
	$(DOBJ)parcomm_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)string_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)partition_factory_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)mesh_collection_mod.o \
	$(DOBJ)tiles_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)domain_mod.o: src/domain/domain_mod.f90 \
	$(DOBJ)topology_mod.o \
	$(DOBJ)metric_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)orography_mod.o \
	$(DOBJ)mesh_collection_mod.o \
	$(DOBJ)abstract_generic_config_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)config_orography_mod.o: src/domain/orography/config_orography_mod.f90 \
	$(DOBJ)config_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)orography_factory_mod.o: src/domain/orography/orography_factory_mod.f90 \
	$(DOBJ)orography_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)orography_test_field_mod.o \
	$(DOBJ)schar_orography_field_mod.o \
	$(DOBJ)mirw3d_orography_field_mod.o \
	$(DOBJ)baroclinic_instability_test_orography_mod.o \
	$(DOBJ)test_fieds_3d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)orography_mod.o: src/domain/orography/orography_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)linked_list_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)partition_mod.o: src/parallel/partition_mod.f90 \
	$(DOBJ)tiles_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)smart_array_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)partition_factory_mod.o: src/parallel/partition_factory_mod.f90 \
	$(DOBJ)partition_mod.o \
	$(DOBJ)tiles_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)exchange_halo_ch_mod.o: src/parallel/exchange_halo_Ch_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)buffer_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)parcomm_factory_mod.o: src/parallel/parcomm_factory_mod.f90 \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)exchange_gather_mod.o: src/parallel/exchange_gather_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)buffer_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)exchange_halo_c_mod.o: src/parallel/exchange_halo_C_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)buffer_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)parcomm_mod.o: src/parallel/parcomm_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)buffer_mod.o: src/parallel/buffer_mod.f90 \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)exchange_factory_mod.o: src/parallel/exchange_factory_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)partition_mod.o \
	$(DOBJ)topology_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)tiles_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)exchange_halo_ch_mod.o \
	$(DOBJ)exchange_halo_c_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)exchange_gather_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)exchange_abstract_mod.o: src/parallel/exchange_abstract_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)exchange_halo_mod.o: src/parallel/exchange_halo_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)buffer_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)get_latlon_main.o: src/mini_apps/get_latlon/get_latlon_main.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)domain_factory_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_norm_mod.o: src/differential_operators/sbp_operators/sbp_norm_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_diff_c2i_21_mod.o: src/differential_operators/sbp_operators/sbp_diff_c2i_21_mod.f90 \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_diff_c2i_63_mod.o: src/differential_operators/sbp_operators/sbp_diff_c2i_63_mod.f90 \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sbp_diff_i2c_63_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_factory_mod.o: src/differential_operators/sbp_operators/sbp_factory_mod.f90 \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sbp_operators_collection_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_interp_i2c_21_mod.o: src/differential_operators/sbp_operators/sbp_interp_i2c_21_mod.f90 \
	$(DOBJ)sbp_interp_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_interp_c2i_63_mod.o: src/differential_operators/sbp_operators/sbp_interp_c2i_63_mod.f90 \
	$(DOBJ)sbp_interp_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sbp_interp_i2c_63_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_diff_mod.o: src/differential_operators/sbp_operators/sbp_diff_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_diff_21_mod.o: src/differential_operators/sbp_operators/sbp_diff_21_mod.f90 \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_diff_43_mod.o: src/differential_operators/sbp_operators/sbp_diff_43_mod.f90 \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_interp_mod.o: src/differential_operators/sbp_operators/sbp_interp_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_interp_i2c_42_mod.o: src/differential_operators/sbp_operators/sbp_interp_i2c_42_mod.f90 \
	$(DOBJ)sbp_interp_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_interp_c2i_42_mod.o: src/differential_operators/sbp_operators/sbp_interp_c2i_42_mod.f90 \
	$(DOBJ)sbp_interp_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sbp_interp_i2c_42_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_diff_c2i_42_mod.o: src/differential_operators/sbp_operators/sbp_diff_c2i_42_mod.f90 \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_interp_i2c_63_mod.o: src/differential_operators/sbp_operators/sbp_interp_i2c_63_mod.f90 \
	$(DOBJ)sbp_interp_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_diff_i2c_63_mod.o: src/differential_operators/sbp_operators/sbp_diff_i2c_63_mod.f90 \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_operators_collection_mod.o: src/differential_operators/sbp_operators/sbp_operators_collection_mod.f90 \
	$(DOBJ)sbp_diff_i2c_63_mod.o \
	$(DOBJ)sbp_diff_c2i_63_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_operator_mod.o: src/differential_operators/sbp_operators/sbp_operator_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_diff_63_mod.o: src/differential_operators/sbp_operators/sbp_diff_63_mod.f90 \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_diff_i2c_21_mod.o: src/differential_operators/sbp_operators/sbp_diff_i2c_21_mod.f90 \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_diff_i2c_42_mod.o: src/differential_operators/sbp_operators/sbp_diff_i2c_42_mod.f90 \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_diff_42_mod.o: src/differential_operators/sbp_operators/sbp_diff_42_mod.f90 \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_interp_c2i_21_mod.o: src/differential_operators/sbp_operators/sbp_interp_c2i_21_mod.f90 \
	$(DOBJ)sbp_interp_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)laplace_true_hor_mod.o: src/differential_operators/3d/laplace/laplace_true_hor_mod.f90 \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)abstract_co2contra_3d_mod.o \
	$(DOBJ)abstract_div_3d_mod.o \
	$(DOBJ)abstract_grad_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)laplace_true_hor_factory_mod.o: src/differential_operators/3d/laplace/laplace_true_hor_factory_mod.f90 \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)laplace_true_hor_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)grad_3d_factory_mod.o \
	$(DOBJ)div_3d_factory_mod.o \
	$(DOBJ)co2contra_3d_factory_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)interpolators_w2uv_mod.o: src/differential_operators/3d/interpolators/interpolators_w2uv_mod.f90 \
	$(DOBJ)abstract_interpolators3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)interpolators_uv2w_mod.o: src/differential_operators/3d/interpolators/interpolators_uv2w_mod.f90 \
	$(DOBJ)abstract_interpolators3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)interpolator_w2uv_factory_mod.o: src/differential_operators/3d/interpolators/interpolator_w2uv_factory_mod.f90 \
	$(DOBJ)abstract_interpolators3d_mod.o \
	$(DOBJ)interpolators_w2uv_mod.o \
	$(DOBJ)vertical_operator_factory_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_interpolators3d_mod.o: src/differential_operators/3d/interpolators/abstract_interpolators3d_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)interpolator_uv2w_factory_mod.o: src/differential_operators/3d/interpolators/interpolator_uv2w_factory_mod.f90 \
	$(DOBJ)abstract_interpolators3d_mod.o \
	$(DOBJ)interpolators_uv2w_mod.o \
	$(DOBJ)vertical_operator_factory_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)div_3d_hor_vert_mod.o: src/differential_operators/3d/divergence/div_3d_hor_vert_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_div_3d_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)vec_math_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)div_3d_factory_mod.o: src/differential_operators/3d/divergence/div_3d_factory_mod.f90 \
	$(DOBJ)abstract_div_3d_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)div_3d_hor_vert_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)vertical_operator_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_div_3d_mod.o: src/differential_operators/3d/divergence/abstract_div_3d_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)mixvec_transform_staggered_mod.o: src/differential_operators/3d/mixvec_transform/mixvec_transform_staggered_mod.f90 \
	$(DOBJ)abstract_mixvec_transform_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_mixvec_transform_mod.o: src/differential_operators/3d/mixvec_transform/abstract_mixvec_transform_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)mixvec_transform_hor_colocated_mod.o: src/differential_operators/3d/mixvec_transform/mixvec_transform_hor_colocated_mod.f90 \
	$(DOBJ)abstract_mixvec_transform_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)config_mixvec_transform_mod.o: src/differential_operators/3d/mixvec_transform/config_mixvec_transform_mod.f90 \
	$(DOBJ)config_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)mixvec_transform_factory_mod.o: src/differential_operators/3d/mixvec_transform/mixvec_transform_factory_mod.f90 \
	$(DOBJ)abstract_mixvec_transform_mod.o \
	$(DOBJ)mixvec_transform_colocated_mod.o \
	$(DOBJ)mixvec_transform_hor_colocated_mod.o \
	$(DOBJ)mixvec_transform_staggered_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)config_mod.o \
	$(DOBJ)vertical_operator_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)config_mixvec_transform_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)mixvec_transform_colocated_mod.o: src/differential_operators/3d/mixvec_transform/mixvec_transform_colocated_mod.f90 \
	$(DOBJ)abstract_mixvec_transform_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_massflux_3d_mod.o: src/differential_operators/3d/massflux/abstract_massflux_3d_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)massflux_3d_mod.o: src/differential_operators/3d/massflux/massflux_3d_mod.f90 \
	$(DOBJ)abstract_massflux_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)abstract_massflux_vertical_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)massflux_3d_factory_mod.o: src/differential_operators/3d/massflux/massflux_3d_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_massflux_3d_mod.o \
	$(DOBJ)massflux_factory_mod.o \
	$(DOBJ)massflux_vertical_factory_mod.o \
	$(DOBJ)massflux_3d_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)string_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_grad_3d_mod.o: src/differential_operators/3d/gradient/abstract_grad_3d_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)grad_3d_hor_vert_mod.o: src/differential_operators/3d/gradient/grad_3d_hor_vert_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_grad_3d_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)grad_3d_factory_mod.o: src/differential_operators/3d/gradient/grad_3d_factory_mod.f90 \
	$(DOBJ)abstract_grad_3d_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grad_3d_hor_vert_mod.o \
	$(DOBJ)grad_factory_mod.o \
	$(DOBJ)vertical_operator_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)coriolis_3d_factory_mod.o: src/differential_operators/3d/coriolis/coriolis_3d_factory_mod.f90 \
	$(DOBJ)abstract_coriolis_3d_mod.o \
	$(DOBJ)shallow_atm_coriolis_colocated_mod.o \
	$(DOBJ)shallow_atm_coriolis_c_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)shallow_atm_coriolis_c_mod.o: src/differential_operators/3d/coriolis/shallow_atm_coriolis_C_mod.f90 \
	$(DOBJ)abstract_coriolis_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_coriolis_3d_mod.o: src/differential_operators/3d/coriolis/abstract_coriolis_3d_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)shallow_atm_coriolis_colocated_mod.o: src/differential_operators/3d/coriolis/shallow_atm_coriolis_colocated_mod.f90 \
	$(DOBJ)abstract_coriolis_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)co2contra_3d_factory_mod.o: src/differential_operators/3d/co2contra/co2contra_3d_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_co2contra_3d_mod.o \
	$(DOBJ)co2contra_3d_colocated_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)co2contra_3d_cgrid_mod.o \
	$(DOBJ)vertical_operator_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)co2contra_3d_h_colocated_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)co2contra_3d_colocated_mod.o: src/differential_operators/3d/co2contra/co2contra_3d_colocated_mod.f90 \
	$(DOBJ)abstract_co2contra_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)co2contra_3d_h_colocated_mod.o: src/differential_operators/3d/co2contra/co2contra_3d_h_colocated_mod.f90 \
	$(DOBJ)abstract_co2contra_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_co2contra_3d_mod.o: src/differential_operators/3d/co2contra/abstract_co2contra_3d_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)co2contra_3d_cgrid_mod.o: src/differential_operators/3d/co2contra/co2contra_3d_Cgrid_mod.f90 \
	$(DOBJ)abstract_co2contra_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)advection_p3d_mod.o: src/differential_operators/3d/advection/advection_p3d_mod.f90 \
	$(DOBJ)abstract_scalar_advection3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)abstract_v_nabla_mod.o \
	$(DOBJ)abstract_adv_z_mod.o \
	$(DOBJ)exchange_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)scalar_advection_factory_mod.o: src/differential_operators/3d/advection/scalar_advection_factory_mod.f90 \
	$(DOBJ)abstract_scalar_advection3d_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)v_nabla_factory_mod.o \
	$(DOBJ)adv_z_factory_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)interpolator_uv2w_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)vertical_operator_factory_mod.o \
	$(DOBJ)advection_p3d_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)advection_w3d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)config_advection_3d_mod.o: src/differential_operators/3d/advection/config_advection_3d_mod.f90 \
	$(DOBJ)config_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_vector_advection3d_mod.o: src/differential_operators/3d/advection/abstract_vector_advection3d_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)vector_advection3d_factory_mod.o: src/differential_operators/3d/advection/vector_advection3d_factory_mod.f90 \
	$(DOBJ)abstract_vector_advection3d_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)shallow_atm_vecadv_mod.o \
	$(DOBJ)vector_advection_factory_mod.o \
	$(DOBJ)adv_z_factory_mod.o \
	$(DOBJ)interpolator_w2uv_factory_mod.o \
	$(DOBJ)scalar_advection_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)config_advection_3d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_scalar_advection3d_mod.o: src/differential_operators/3d/advection/abstract_scalar_advection3d_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)advection_w3d_mod.o: src/differential_operators/3d/advection/advection_w3d_mod.f90 \
	$(DOBJ)abstract_scalar_advection3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)abstract_interpolators3d_mod.o \
	$(DOBJ)abstract_v_nabla_mod.o \
	$(DOBJ)abstract_adv_z_mod.o \
	$(DOBJ)exchange_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)shallow_atm_vecadv_mod.o: src/differential_operators/3d/advection/shallow_atm_vecadv_mod.f90 \
	$(DOBJ)abstract_vector_advection3d_mod.o \
	$(DOBJ)abstract_vector_advection_mod.o \
	$(DOBJ)abstract_scalar_advection3d_mod.o \
	$(DOBJ)abstract_interpolators3d_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_adv_z_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)hor_christofel_factory_mod.o: src/differential_operators/horizontal/Christofel/hor_Christofel_factory_mod.f90 \
	$(DOBJ)abstract_hor_christofel_mod.o \
	$(DOBJ)hor_christofel_shallow_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)metric_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_hor_christofel_mod.o: src/differential_operators/horizontal/Christofel/abstract_hor_Christofel_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)hor_christofel_shallow_mod.o: src/differential_operators/horizontal/Christofel/hor_Christofel_shallow_mod.f90 \
	$(DOBJ)abstract_hor_christofel_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)laplace_ch_halo_mod.o: src/differential_operators/horizontal/laplace/laplace_ch_halo_mod.f90 \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_laplace_mod.o: src/differential_operators/horizontal/laplace/abstract_laplace_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)divgrad_laplace_mod.o: src/differential_operators/horizontal/laplace/divgrad_laplace_mod.f90 \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_co2contra_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)laplace_factory_mod.o: src/differential_operators/horizontal/laplace/laplace_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)laplace_true_hor_factory_mod.o \
	$(DOBJ)divgrad_laplace_mod.o \
	$(DOBJ)grad_factory_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)co2contra_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)laplace_ch_halo_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)laplace_ah_sbp21_narrow_mod.o \
	$(DOBJ)laplace_ah_sbp42_narrow_mod.o \
	$(DOBJ)laplace_ah_sbp63_narrow_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)sbp_diff_21_mod.o \
	$(DOBJ)sbp_diff_42_mod.o \
	$(DOBJ)sbp_diff_63_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)laplace_ah_sbp21_narrow_mod.o: src/differential_operators/horizontal/laplace/laplace_Ah_sbp21_narrow_mod.f90 \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)laplace_ah_sbp63_narrow_mod.o: src/differential_operators/horizontal/laplace/laplace_Ah_sbp63_narrow_mod.f90 \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)laplace_ah_sbp21_narrow_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)laplace_ah_sbp42_narrow_mod.o: src/differential_operators/horizontal/laplace/laplace_Ah_sbp42_narrow_mod.f90 \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)laplace_ah_sbp21_narrow_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_v_nabla_mod.o: src/differential_operators/horizontal/v_dot_nabla/abstract_v_nabla_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)v_nabla_mod.o: src/differential_operators/horizontal/v_dot_nabla/v_nabla_mod.f90 \
	$(DOBJ)abstract_v_nabla_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)v_nabla_sbp_factory_mod.o: src/differential_operators/horizontal/v_dot_nabla/v_nabla_sbp_factory_mod.f90 \
	$(DOBJ)v_nabla_ah_sbp_mod.o \
	$(DOBJ)sbp_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)v_nabla_ah_sbp_mod.o: src/differential_operators/horizontal/v_dot_nabla/v_nabla_Ah_sbp_mod.f90 \
	$(DOBJ)abstract_v_nabla_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)v_nabla_factory_mod.o: src/differential_operators/horizontal/v_dot_nabla/v_nabla_factory_mod.f90 \
	$(DOBJ)abstract_v_nabla_mod.o \
	$(DOBJ)v_nabla_mod.o \
	$(DOBJ)v_nabla_sbp_factory_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_div_mod.o: src/differential_operators/horizontal/divergence/abstract_div_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)div_a2_mod.o: src/differential_operators/horizontal/divergence/div_a2_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)div_c2_mod.o: src/differential_operators/horizontal/divergence/div_c2_mod.f90 \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)div_sbp_proj_mod.o: src/differential_operators/horizontal/divergence/div_sbp_proj_mod.f90 \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)div_sbp_sat_mod.o: src/differential_operators/horizontal/divergence/div_sbp_sat_mod.f90 \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_halo_mod.o \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)div_factory_mod.o: src/differential_operators/horizontal/divergence/div_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)div_c2_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)div_sbp_proj_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)sbp_diff_21_mod.o \
	$(DOBJ)sbp_diff_42_mod.o \
	$(DOBJ)sbp_diff_43_mod.o \
	$(DOBJ)sbp_diff_63_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)div_sbp_sat_mod.o \
	$(DOBJ)sbp_diff_i2c_21_mod.o \
	$(DOBJ)sbp_diff_i2c_42_mod.o \
	$(DOBJ)sbp_diff_i2c_63_mod.o \
	$(DOBJ)sbp_diff_c2i_21_mod.o \
	$(DOBJ)sbp_diff_c2i_42_mod.o \
	$(DOBJ)sbp_diff_c2i_63_mod.o \
	$(DOBJ)div_a2_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)massflux_ch_mod.o: src/differential_operators/horizontal/massflux/massflux_ch_mod.f90 \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)massflux_factory_mod.o: src/differential_operators/horizontal/massflux/massflux_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)massflux_colocated_mod.o \
	$(DOBJ)massflux_cgrid_mod.o \
	$(DOBJ)massflux_ch_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)sbp_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_massflux_mod.o: src/differential_operators/horizontal/massflux/abstract_massflux_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)massflux_colocated_mod.o: src/differential_operators/horizontal/massflux/massflux_colocated_mod.f90 \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)massflux_cgrid_mod.o: src/differential_operators/horizontal/massflux/massflux_Cgrid_mod.f90 \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)grad_ch_halo_mod.o: src/differential_operators/horizontal/gradient/grad_ch_halo_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)sbp_operator_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)grad_c2_ecs_mod.o: src/differential_operators/horizontal/gradient/grad_c2_ecs_mod.f90 \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)grad_sbp_proj_mod.o: src/differential_operators/horizontal/gradient/grad_sbp_proj_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)grad_sbp_sat_mod.o: src/differential_operators/horizontal/gradient/grad_sbp_sat_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_grad_mod.o: src/differential_operators/horizontal/gradient/abstract_grad_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)grad_factory_mod.o: src/differential_operators/horizontal/gradient/grad_factory_mod.f90 \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grad_c2_ecs_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)grad_sbp_proj_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)sbp_diff_21_mod.o \
	$(DOBJ)sbp_diff_42_mod.o \
	$(DOBJ)sbp_diff_43_mod.o \
	$(DOBJ)sbp_diff_63_mod.o \
	$(DOBJ)grad_sbp_sat_mod.o \
	$(DOBJ)sbp_diff_c2i_21_mod.o \
	$(DOBJ)sbp_diff_c2i_42_mod.o \
	$(DOBJ)sbp_diff_c2i_63_mod.o \
	$(DOBJ)sbp_diff_i2c_21_mod.o \
	$(DOBJ)sbp_diff_i2c_42_mod.o \
	$(DOBJ)sbp_diff_i2c_63_mod.o \
	$(DOBJ)grad_a2_mod.o \
	$(DOBJ)grad_ch_halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)grad_a2_mod.o: src/differential_operators/horizontal/gradient/grad_a2_mod.f90 \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)hordiff_cgrid_mod.o: src/differential_operators/horizontal/hordiff/hordiff_Cgrid_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_hordiff_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_grad_perp_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)vec_math_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_hordiff_mod.o: src/differential_operators/horizontal/hordiff/abstract_hordiff_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)hordiff_factory_mod.o: src/differential_operators/horizontal/hordiff/hordiff_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)abstract_hordiff_mod.o \
	$(DOBJ)hordiff_cgrid_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)grad_factory_mod.o \
	$(DOBJ)co2contra_factory_mod.o \
	$(DOBJ)curl_factory_mod.o \
	$(DOBJ)grad_perp_factory_mod.o \
	$(DOBJ)hordiff_scalar_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)laplace_factory_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)hordiff_colocated_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)hordiff_no_metric_mod.o: src/differential_operators/horizontal/hordiff/hordiff_no_metric_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)abstract_hordiff_mod.o \
	$(DOBJ)abstract_laplace_mod.o \
	$(DOBJ)vec_math_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)hordiff_colocated_mod.o: src/differential_operators/horizontal/hordiff/hordiff_colocated_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_hordiff_mod.o \
	$(DOBJ)vec_math_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)hordiff_scalar_mod.o: src/differential_operators/horizontal/hordiff/hordiff_scalar_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)abstract_hordiff_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)abstract_laplace_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)grad_perp_factory_mod.o: src/differential_operators/horizontal/grad_perp/grad_perp_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_grad_perp_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grad_perp_sbp_sat_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)sbp_diff_i2c_21_mod.o \
	$(DOBJ)sbp_diff_i2c_42_mod.o \
	$(DOBJ)sbp_diff_i2c_63_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)grad_perp_sbp_sat_mod.o: src/differential_operators/horizontal/grad_perp/grad_perp_sbp_sat_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_grad_perp_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_grad_perp_mod.o: src/differential_operators/horizontal/grad_perp/abstract_grad_perp_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ke_factory_mod.o: src/differential_operators/horizontal/kinetic_energy/KE_factory_mod.f90 \
	$(DOBJ)abstract_ke_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)ke_colocated_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)ke_cgrid_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)halo_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_ke_mod.o: src/differential_operators/horizontal/kinetic_energy/abstract_KE_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ke_colocated_mod.o: src/differential_operators/horizontal/kinetic_energy/KE_colocated_mod.f90 \
	$(DOBJ)abstract_ke_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)ke_cgrid_mod.o: src/differential_operators/horizontal/kinetic_energy/KE_Cgrid_mod.f90 \
	$(DOBJ)abstract_ke_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)coriolis_factory_mod.o: src/differential_operators/horizontal/coriolis/coriolis_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)coriolis_cgrid_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)interpolator_w2h_factory_mod.o \
	$(DOBJ)coriolis_cgrid_noncons_mod.o \
	$(DOBJ)sbp_factory_mod.o \
	$(DOBJ)co2contra_factory_mod.o \
	$(DOBJ)coriolis_ch_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)coriolis_colocated_mod.o \
	$(DOBJ)metric_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)coriolis_cgrid_noncons_mod.o: src/differential_operators/horizontal/coriolis/coriolis_Cgrid_noncons_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)abstract_co2contra_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)coriolis_cgrid_mod.o: src/differential_operators/horizontal/coriolis/coriolis_Cgrid_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)interpolator_w2h_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)coriolis_colocated_mod.o: src/differential_operators/horizontal/coriolis/coriolis_colocated_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_coriolis_mod.o: src/differential_operators/horizontal/coriolis/abstract_coriolis_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)coriolis_c_contra_lincons_mod.o: src/differential_operators/horizontal/coriolis/coriolis_C_contra_lincons_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)vec_math_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)coriolis_ch_mod.o: src/differential_operators/horizontal/coriolis/coriolis_Ch_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_coriolis_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)interpolator_w2h_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_abstract_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)p2uv_colocated_mod.o: src/differential_operators/horizontal/interpolator_2d/p2uv_colocated_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)interpolator_p2uv_mod.o: src/differential_operators/horizontal/interpolator_2d/interpolator_p2uv_mod.f90 \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sbp_interp_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)interpolator_uv2p_mod.o: src/differential_operators/horizontal/interpolator_2d/interpolator_uv2p_mod.f90 \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sbp_interp_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)interpolator_uv2q_mod.o: src/differential_operators/horizontal/interpolator_2d/interpolator_uv2q_mod.f90 \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sbp_interp_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)interpolator_q2uv_mod.o: src/differential_operators/horizontal/interpolator_2d/interpolator_q2uv_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sbp_interp_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)interpolator_w2h_mod.o: src/differential_operators/horizontal/interpolator_2d/interpolator_w2h_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)sbp_interp_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_interpolators2d_mod.o: src/differential_operators/horizontal/interpolator_2d/abstract_interpolators2d_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)interpolator2d_factory_mod.o: src/differential_operators/horizontal/interpolator_2d/interpolator2d_factory_mod.f90 \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)interpolator_p2uv_mod.o \
	$(DOBJ)p2uv_colocated_mod.o \
	$(DOBJ)interpolator_uv2p_mod.o \
	$(DOBJ)interpolator_q2uv_mod.o \
	$(DOBJ)interpolator_uv2q_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)sbp_factory_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sbp_interp_mod.o \
	$(DOBJ)sbp_interp_c2i_21_mod.o \
	$(DOBJ)sbp_interp_i2c_21_mod.o \
	$(DOBJ)sbp_interp_i2c_42_mod.o \
	$(DOBJ)sbp_interp_c2i_42_mod.o \
	$(DOBJ)sbp_interp_i2c_63_mod.o \
	$(DOBJ)sbp_interp_c2i_63_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)interpolator_w2h_factory_mod.o: src/differential_operators/horizontal/interpolator_2d/interpolator_w2h_factory_mod.f90 \
	$(DOBJ)interpolator_w2h_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)sbp_factory_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)sbp_interp_i2c_21_mod.o \
	$(DOBJ)sbp_interp_i2c_42_mod.o \
	$(DOBJ)sbp_interp_i2c_63_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)vector_advection_factory_mod.o: src/differential_operators/horizontal/vec_advection/vector_advection_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_vector_advection_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)v_nabla_mod.o \
	$(DOBJ)hor_christofel_factory_mod.o \
	$(DOBJ)vector_advection_ah_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)v_nabla_sbp_factory_mod.o \
	$(DOBJ)vector_advection_c_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)abstract_v_nabla_mod.o \
	$(DOBJ)vector_advection_ch_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_vector_advection_mod.o: src/differential_operators/horizontal/vec_advection/abstract_vector_advection_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)vector_advection_c_mod.o: src/differential_operators/horizontal/vec_advection/vector_advection_C_mod.f90 \
	$(DOBJ)abstract_v_nabla_mod.o \
	$(DOBJ)abstract_vector_advection_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)abstract_hor_christofel_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)vector_advection_ch_mod.o: src/differential_operators/horizontal/vec_advection/vector_advection_Ch_mod.f90 \
	$(DOBJ)abstract_v_nabla_mod.o \
	$(DOBJ)abstract_vector_advection_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)abstract_hor_christofel_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)vector_advection_ah_mod.o: src/differential_operators/horizontal/vec_advection/vector_advection_Ah_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)abstract_vector_advection_mod.o \
	$(DOBJ)abstract_v_nabla_mod.o \
	$(DOBJ)abstract_hor_christofel_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)co2contra_ch_mod.o: src/differential_operators/horizontal/co2contra/co2contra_Ch_mod.f90 \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)co2contra_colocated_mod.o: src/differential_operators/horizontal/co2contra/co2contra_colocated_mod.f90 \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_co2contra_mod.o: src/differential_operators/horizontal/co2contra/abstract_co2contra_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)co2contra_cgrid_mod.o: src/differential_operators/horizontal/co2contra/co2contra_Cgrid_mod.f90 \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o \
	$(DOBJ)tile_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)co2contra_factory_mod.o: src/differential_operators/horizontal/co2contra/co2contra_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_co2contra_mod.o \
	$(DOBJ)co2contra_colocated_mod.o \
	$(DOBJ)co2contra_cgrid_mod.o \
	$(DOBJ)co2contra_ch_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)sbp_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)curl_div_based_mod.o: src/differential_operators/horizontal/curl/curl_div_based_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)curl_factory_mod.o: src/differential_operators/horizontal/curl/curl_factory_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)curl_div_based_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)curl_sbp_sat_mod.o \
	$(DOBJ)sbp_diff_c2i_21_mod.o \
	$(DOBJ)sbp_diff_c2i_42_mod.o \
	$(DOBJ)sbp_diff_c2i_63_mod.o \
	$(DOBJ)sbp_diff_i2c_21_mod.o \
	$(DOBJ)sbp_diff_i2c_42_mod.o \
	$(DOBJ)sbp_diff_i2c_63_mod.o \
	$(DOBJ)exchange_factory_mod.o \
	$(DOBJ)halo_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)curl_sbp_sat_mod.o: src/differential_operators/horizontal/curl/curl_sbp_sat_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_curl_mod.o \
	$(DOBJ)sbp_diff_mod.o \
	$(DOBJ)exchange_abstract_mod.o \
	$(DOBJ)halo_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_curl_mod.o: src/differential_operators/horizontal/curl/abstract_curl_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)massflux_div_mod.o: src/differential_operators/horizontal/flux_div/massflux_div_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)abstract_flux_div_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_flux_div_mod.o: src/differential_operators/horizontal/flux_div/abstract_flux_div_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)skew_symmetric_flux_div_mod.o: src/differential_operators/horizontal/flux_div/skew_symmetric_flux_div_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_div_mod.o \
	$(DOBJ)abstract_grad_mod.o \
	$(DOBJ)abstract_massflux_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)vec_math_mod.o \
	$(DOBJ)abstract_flux_div_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)flux_div_factory_mod.o: src/differential_operators/horizontal/flux_div/flux_div_factory_mod.f90 \
	$(DOBJ)abstract_flux_div_mod.o \
	$(DOBJ)skew_symmetric_flux_div_mod.o \
	$(DOBJ)massflux_div_mod.o \
	$(DOBJ)grad_factory_mod.o \
	$(DOBJ)div_factory_mod.o \
	$(DOBJ)massflux_factory_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)massflux_vertical_factory_mod.o: src/differential_operators/vertical/massflux_vertical_factory_mod.f90 \
	$(DOBJ)abstract_massflux_vertical_mod.o \
	$(DOBJ)massflux_vertical_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)massflux_vertical_mod.o: src/differential_operators/vertical/massflux_vertical_mod.f90 \
	$(DOBJ)abstract_massflux_vertical_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)adv_z_mod.o: src/differential_operators/vertical/adv_z_mod.f90 \
	$(DOBJ)abstract_adv_z_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)identity_vertical_operator_mod.o: src/differential_operators/vertical/identity_vertical_operator_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)sbp_vertical_operator_mod.o: src/differential_operators/vertical/sbp_vertical_operator_mod.f90 \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)sbp_operator_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_adv_z_mod.o: src/differential_operators/vertical/abstract_adv_z_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)adv_z_factory_mod.o: src/differential_operators/vertical/adv_z_factory_mod.f90 \
	$(DOBJ)abstract_adv_z_mod.o \
	$(DOBJ)adv_z_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)vertical_operator_factory_mod.o: src/differential_operators/vertical/vertical_operator_factory_mod.f90 \
	$(DOBJ)abstract_vertical_operator_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)identity_vertical_operator_mod.o \
	$(DOBJ)sbp_vertical_operator_mod.o \
	$(DOBJ)sbp_factory_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_vertical_operator_mod.o: src/differential_operators/vertical/abstract_vertical_operator_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_massflux_vertical_mod.o: src/differential_operators/vertical/abstract_massflux_vertical_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)vertical_transform_factory_mod.o: src/metric/vertical_transform_factory_mod.f90 \
	$(DOBJ)abstract_vertical_transform_mod.o \
	$(DOBJ)vertical_transform_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)metric_mod.o: src/metric/metric_mod.f90 \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)orography_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)config_metric_mod.o: src/metric/config_metric_mod.f90 \
	$(DOBJ)config_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)metric_2d_mod.o: src/metric/metric_2d_mod.f90 \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)metric_factory_mod.o: src/metric/metric_factory_mod.f90 \
	$(DOBJ)metric_mod.o \
	$(DOBJ)topology_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)config_metric_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)string_mod.o \
	$(DOBJ)shallow_atm_metric_mod.o \
	$(DOBJ)vertical_transform_factory_mod.o \
	$(DOBJ)metric_2d_mod.o \
	$(DOBJ)cubed_sphere_topology_mod.o \
	$(DOBJ)metric_2d_ecs_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_vertical_transform_mod.o: src/metric/abstract_vertical_transform_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)shallow_atm_metric_mod.o: src/metric/shallow_atm_metric_mod.f90 \
	$(DOBJ)metric_mod.o \
	$(DOBJ)metric_2d_mod.o \
	$(DOBJ)abstract_vertical_transform_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)orography_mod.o \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)metric_2d_ecs_mod.o: src/metric/metric_2d_ecs_mod.f90 \
	$(DOBJ)metric_2d_mod.o \
	$(DOBJ)cubed_sphere_topology_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)vertical_transform_mod.o: src/metric/vertical_transform_mod.f90 \
	$(DOBJ)abstract_vertical_transform_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)latlon_regrid_mod.o: src/regridders/latlon_regrid_mod.f90 \
	$(DOBJ)abstract_regridder_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)abstract_interpolators2d_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)regrid_factory_mod.o: src/regridders/regrid_factory_mod.f90 \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)abstract_regridder_mod.o \
	$(DOBJ)latlon_regrid_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)halo_mod.o \
	$(DOBJ)halo_factory_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)tiles_mod.o \
	$(DOBJ)interpolator2d_factory_mod.o \
	$(DOBJ)latlon_grid_generator_mod.o \
	$(DOBJ)array_tools_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_regridder_mod.o: src/regridders/abstract_regridder_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)stvec_flexible_mod.o: src/stvec/stvec_flexible_mod.f90 \
	$(DOBJ)grid_field_collection_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)string_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)basic_collection_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)stvec_mod.o: src/stvec/stvec_mod.f90 \
	$(DOBJ)container_abstract_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)stvec_flexible_factory_mod.o: src/stvec/stvec_flexible_factory_mod.f90 \
	$(DOBJ)stvec_flexible_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)stvec_mod.o \
	$(DOBJ)string_mod.o \
	$(DOBJ)grid_field_factory_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)container_abstract_mod.o: src/stvec/container_abstract_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)assembly_particles_mod.o: src/particles/assembly_particles_mod.f90 \
	$(DOBJ)smart_array_mod.o \
	$(DOBJ)particle_values_mod.o \
	$(DOBJ)distribute_particles_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)particle_values_mod.o: src/particles/particle_values_mod.f90 \
	$(DOBJ)string_mod.o \
	$(DOBJ)particles_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)particles_mod.o: src/particles/particles_mod.f90 \
	$(DOBJ)smart_array_mod.o \
	$(DOBJ)string_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)distribute_particles_mod.o: src/particles/distribute_particles_mod.f90 \
	$(DOBJ)particles_mod.o \
	$(DOBJ)smart_array_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)abstract_generic_config_mod.o \
	$(DOBJ)tile_mod.o \
	$(DOBJ)tiles_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)bicubic_interp_mod.o: src/particles/interpolations/bicubic_interp_mod.f90 \
	$(DOBJ)abstract_particle_interpolation_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)particles_mod.o \
	$(DOBJ)particle_values_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)array_tools_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)particle_interp_factory_mod.o: src/particles/interpolations/particle_interp_factory_mod.f90 \
	$(DOBJ)abstract_particle_interpolation_mod.o \
	$(DOBJ)bilinear_interp_mod.o \
	$(DOBJ)bicubic_interp_mod.o \
	$(DOBJ)domain_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_particle_interpolation_mod.o: src/particles/interpolations/abstract_particle_interpolation_mod.f90 \
	$(DOBJ)domain_mod.o \
	$(DOBJ)particles_mod.o \
	$(DOBJ)particle_values_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)string_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)bilinear_interp_mod.o: src/particles/interpolations/bilinear_interp_mod.f90 \
	$(DOBJ)abstract_particle_interpolation_mod.o \
	$(DOBJ)domain_mod.o \
	$(DOBJ)particles_mod.o \
	$(DOBJ)particle_values_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)array_tools_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)latlon_functions_mod.o: src/test_fields/latlon_functions_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)test_fields_mod.o: src/test_fields/test_fields_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)sph_coords_mod.o \
	$(DOBJ)latlon_functions_mod.o \
	$(DOBJ)barotropic_instability_u_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)isotermal_profile_mod.o: src/test_fields/vertical_thermodynamic_profiles/isotermal_profile_mod.f90 \
	$(DOBJ)abstract_vertical_profile_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)const_n_profile_mod.o: src/test_fields/vertical_thermodynamic_profiles/const_N_profile_mod.f90 \
	$(DOBJ)abstract_vertical_profile_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)abstract_vertical_profile_mod.o: src/test_fields/vertical_thermodynamic_profiles/abstract_vertical_profile_mod.f90
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)orography_test_field_mod.o: src/test_fields/test_fields_3d/orography_test_field_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)schar_orography_field_mod.o: src/test_fields/test_fields_3d/Schar_orography_field_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)div3d_test_field_mod.o: src/test_fields/test_fields_3d/div3d_test_field_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)vertical_div_test_field_mod.o: src/test_fields/test_fields_3d/vertical_div_test_field_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)grad3d_test_field_mod.o: src/test_fields/test_fields_3d/grad3d_test_field_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)mirw3d_orography_field_mod.o: src/test_fields/test_fields_3d/MIRW3d_orography_field_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)vertical_test_field_mod.o: src/test_fields/test_fields_3d/vertical_test_field_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)abstract_vertical_profile_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)test_fieds_3d_mod.o: src/test_fields/test_fields_3d/test_fieds_3d_mod.f90 \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)baroclinic_instability_solve_h_sigma_mod.o: src/test_fields/test_fields_3d/baroclinic_instability/baroclinic_instability_solve_h_sigma_mod.f90 \
	$(DOBJ)baroclinic_instability_test_parameters_mod.o \
	$(DOBJ)baroclinic_instability_test_therm_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)baroclinic_instability_testcase_mod.o: src/test_fields/test_fields_3d/baroclinic_instability/baroclinic_instability_testcase_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)baroclinic_instability_solve_h_sigma_mod.o \
	$(DOBJ)baroclinic_instability_test_parameters_mod.o \
	$(DOBJ)baroclinic_instability_test_therm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)baroclinic_instability_test_parameters_mod.o: src/test_fields/test_fields_3d/baroclinic_instability/baroclinic_instability_test_parameters_mod.f90 \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)baroclinic_instability_test_therm_mod.o: src/test_fields/test_fields_3d/baroclinic_instability/baroclinic_instability_test_therm_mod.f90 \
	$(DOBJ)baroclinic_instability_test_parameters_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)baroclinic_instability_test_orography_mod.o: src/test_fields/test_fields_3d/baroclinic_instability/baroclinic_instability_test_orography_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)const_mod.o \
	$(DOBJ)sph_coords_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)solid_rotation_fields_factory_mod.o: src/test_fields/test_fields_3d/solid_rotation/solid_rotation_fields_factory_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)solid_rotation_therm_mod.o \
	$(DOBJ)solid_rotation_wind_field_mod.o \
	$(DOBJ)parcomm_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)solid_rotation_wind_field_mod.o: src/test_fields/test_fields_3d/solid_rotation/solid_rotation_wind_field_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

$(DOBJ)solid_rotation_therm_mod.o: src/test_fields/test_fields_3d/solid_rotation/solid_rotation_therm_mod.f90 \
	$(DOBJ)test_fieds_3d_mod.o \
	$(DOBJ)mesh_mod.o \
	$(DOBJ)grid_field_mod.o \
	$(DOBJ)const_mod.o
	@echo $(COTEXT)
	@$(FC) $(OPTSC) -Isrc/stuff/config -Isrc/stuff/config/geosci_config  $< -o $@

#phony auxiliary rules
.PHONY : $(MKDIRS)
$(MKDIRS):
	@mkdir -p $@
.PHONY : cleanobj
cleanobj:
	@echo deleting objects
	@rm -fr $(DOBJ)
.PHONY : cleanmod
cleanmod:
	@echo deleting mods
	@rm -fr $(DMOD)
.PHONY : cleanexe
cleanexe:
	@echo deleting exes
	@rm -f $(addprefix $(DEXE),$(EXES))
.PHONY : clean
clean: cleanobj cleanmod
.PHONY : cleanall
cleanall: clean cleanexe
all: $(addprefix $(DEXE),$(EXES))
