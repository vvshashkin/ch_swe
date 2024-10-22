module arithmetic_parser_val_mod

implicit none

integer(kind=4), parameter :: INT_TYPE_ID = 1, REAL_TYPE_ID = 2, UNKNOWN_TYPE_ID = -1

type, abstract :: value_t
end type

type, extends(value_t) :: real_value_t
    real(kind=8) :: val
end type

type, extends(value_t) :: int_value_t
    integer(kind=4) :: val
end type

end module