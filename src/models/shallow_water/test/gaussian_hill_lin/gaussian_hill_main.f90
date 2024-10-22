program gaussian_hill_test

use parcomm_mod, only : init_global_parallel_enviroment, &
                        deinit_global_parallel_enviroment, parcomm_global

use gaussian_hill_mod, only : run_gaussian_hill

call init_global_parallel_enviroment()

call run_gaussian_hill()

call deinit_global_parallel_enviroment()

end program gaussian_hill_test
