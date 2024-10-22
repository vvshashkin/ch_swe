recursive subroutine SUBROUTINE_NAME(this,q,var_name,default,found,error_code)
    class(geosci_config_t),  intent(inout)         :: this
    VAR_TYPE,                intent(out)           :: q
    character(len=*),        intent(in)            :: var_name
#if defined(SUBCONFIG)
    class(generic_config_t), intent(in),  optional :: default
#elif defined(CHAR)
    character(len=*),        intent(in),  optional :: default
#else
    VAR_TYPE,                intent(in),  optional :: default
#endif
    logical,                 intent(out), optional :: found
    integer(kind=4),         intent(out), optional :: error_code

    class(*), pointer :: value
    logical :: crash
    integer(kind=4) :: ind
    class(generic_config_t), allocatable :: subconfig

    ind = index(var_name,"%")
    if(ind /= 0) then
        call this%get(subconfig,var_name(1:ind-1))
        call subconfig%get(q,var_name(ind+1:),default=default,found=found,error_code=error_code)
        return
    end if

    call this%data%get_item(value,var_name)

    if(present(found)) found = .false.
    if(present(error_code)) error_code = NO_ERROR

    if(.not. associated(value)) then

        if(present(default)) then
            q = default
            return
        end if

        if(present(error_code)) then
            error_code = NOT_FOUND_ERROR
            return
        end if

        if(.not. present(found)) &
            call parcomm_global%abort(__FILE__//": variable "//var_name//" not found in config")

        return
    end if

    if(present(found)) found = .true.

    select type(value)
#if defined(STRING_T)
    type is (string_t)
#elif defined(CHAR)
    type is (character(len=*))
#elif defined(SUBCONFIG)
    class is (generic_config_t)
#else
    type is (VAR_TYPE)
#endif
        q = value
#if defined(STRING_T)
    type is(character(len=*))
        q%str = value
#elif defined(CHAR)
    type is(string_t)
        q = value%str
#endif
    class default
        if(present(error_code)) then
          error_code = TYPE_ERROR
        else
         call parcomm_global%abort(__FILE__//": variable "//var_name//" type error")
        end if
    end select

end subroutine
