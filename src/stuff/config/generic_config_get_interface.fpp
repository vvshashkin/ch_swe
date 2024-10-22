#ifndef VAR_TYPE
#define VAR_TYPE integer(kind=4)
#endif

#ifndef INTERFACE_NAME
#define INTERFACE_NAME get_i4
#endif

subroutine INTERFACE_NAME(this,q,var_name,default,found,error_code)
    import string_t, generic_config_t
    class(generic_config_t), intent(inout)         :: this
    VAR_TYPE,                intent(out)           :: q
    character(len=*),        intent(in)            :: var_name
#ifndef CHAR
    VAR_TYPE,                intent(in),  optional :: default
#else
    character(len=*),        intent(in) , optional :: default
#endif
    logical,                 intent(out), optional :: found
    integer(kind=4),         intent(out), optional :: error_code
end subroutine

#undef VAR_TYPE1
