module conv_order_mod

contains

real(kind=8) function calculate_convergence_rate(Ns,e) result(conv)
    integer(kind=4), intent(in) :: Ns(:)
    real(kind=8),    intent(in) :: e(:)

    integer(kind=4) :: n, i
    real(kind=8)    :: x(size(Ns)), y(size(Ns))
    real(kind=8)    :: sx, sx2, sy, syx, det

    n = size(Ns)
    do i = 1, n
        y(i) = log(max(e(i),1e-14))
    end do

    x = log(1.0*Ns)

    !solve min{(y-ax-b)**2} for a:
    !conv rate = -a
    sx = sum(x)
    sy = sum(y)
    sx2 = sum(x**2)
    syx = sum(y*x)
    det = n*sx2-sx**2
    conv =-(n*syx-sx*sy) / det

end function calculate_convergence_rate

end module
