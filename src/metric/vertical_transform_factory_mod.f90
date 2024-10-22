module vertical_transform_factory_mod

use abstract_vertical_transform_mod, only : vertical_transform_t
use vertical_transform_mod,          only : vertical_transform_default_t, &
                                            vertical_transform_quadratic_t
use parcomm_mod,                     only : parcomm_global

contains

function create_vertical_transform(vert_transform_name) result(vert_transform)

    character(len=*), intent(in) :: vert_transform_name
    class(vertical_transform_t), allocatable :: vert_transform
    integer(kind=4) :: name_len

    name_len = len(trim(vert_transform_name))

    if(vert_transform_name == "vertical_transform_default") then
        vert_transform = vertical_transform_default_t()
    else if(vert_transform_name(1:min(28,name_len)) == "vertical_transform_quadratic") then
        vert_transform = create_quadratic_vertical_transform(vert_transform_name)
    else
        call parcomm_global%abort("create_vertical_transform, unknown vertical transform:"// &
                                                                        vert_transform_name)
    end if
end function create_vertical_transform

function create_quadratic_vertical_transform(vert_transform_name) result(vert_transform)
    character(len=*), intent(in) :: vert_transform_name
    class(vertical_transform_t), allocatable :: vert_transform

    real(kind=8)    :: a
    integer(kind=4) :: name_len

    name_len = len(vert_transform_name)

    if(vert_transform_name(min(29,name_len):min(34,name_len)) == "_wlin=") then
        read(vert_transform_name(min(35,name_len):),*) a
    else
        call parcomm_global%abort("create_vertical_transform, incorrect quadratic vertical transform name:"// &
                                                                        vert_transform_name)
    end if

    vert_transform = vertical_transform_quadratic_t(a,1.0_8-a)
end function

end module vertical_transform_factory_mod
