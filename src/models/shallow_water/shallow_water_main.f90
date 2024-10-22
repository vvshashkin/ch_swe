program shallow_water_main

    use parcomm_mod,          only : init_global_parallel_enviroment,   &
                                     deinit_global_parallel_enviroment, &
                                     parcomm_global
    use cmd_args_mod,         only : cmd_arg_t, get_cmd_args
    use generic_config_mod,   only : generic_config_t
    use config_tools_mod,     only : parse_config_file

    use shallow_water_mod,    only : run_shallow_water_model

    implicit none

    type(cmd_arg_t),         allocatable :: cmd_args(:)
    integer(kind=4)                      ::nargs
    character(len=:),        allocatable :: config_file_name
    class(generic_config_t), allocatable :: config

    call init_global_parallel_enviroment()

    call get_cmd_args(cmd_args, nargs)

    config_file_name = "swm_model.cfg"
    if(nargs >= 2) config_file_name = cmd_args(2)%str

    config = parse_config_file(config_file_name, parcomm_global)

    call run_shallow_water_model(config)

    call deinit_global_parallel_enviroment()

end program