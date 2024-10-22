recursive subroutine SUBROUTINE_NAME(this,q,var_name,default,found,error_code)
    class(geosci_config_t),  intent(inout)         :: this
    VAR_TYPE, allocatable,   intent(out)           :: q(:)
    character(len=*),        intent(in)            :: var_name
    VAR_TYPE,                intent(in),  optional :: default(:)
    logical,                 intent(out), optional :: found
    integer(kind=4),         intent(out), optional :: error_code

    logical         :: found_loc
    integer(kind=4) :: length, i
    class(generic_config_t), allocatable :: array_conf, element_conf

    call this%get(array_conf,var_name,error_code=error_code,found=found_loc)

    if(present(error_code)) then
        if(error_code /= NO_ERROR) return
        error_code = NO_ERROR
    end if

    if(present(found)) found = found_loc

    if(.not. found_loc) then

        if(present(default)) then
            q = default
            return
        end if

        if(present(found)) return

        if(present(error_code)) then
            error_code = NOT_FOUND_ERROR
        end if

        return
    end if

    call this%get(length,var_name//"%length",error_code=error_code)
    if(present(error_code)) then
        if(error_code /= NO_ERROR) return
    end if

#ifndef SUBCONFIG
    allocate(q(1:length))
    do i = 1, length
        call this%get(q(i),var_name//"%"//integer_to_str(i),error_code=error_code)
        if(present(error_code)) then
            if(error_code /= NO_ERROR) return
        end if
    end do
#else
    allocate(geosci_config_t :: q(1:length))
    do i = 1, length

        call this%get(element_conf,var_name//"%"//integer_to_str(i),error_code=error_code)
        if(present(error_code)) then
            if(error_code /= NO_ERROR) return
        end if

        select type (element_conf)
        type is (geosci_config_t)
            select type (array_element=>q(i))
            type is (geosci_config_t)
                array_element = element_conf
            end select
        class default
            if(present(error_code)) then
                error_code = TYPE_ERROR
                return
            end if
            call parcomm_global%abort("geosci_config_t, get_subconfig_array type error")
        end select
    end do
#endif
end subroutine
