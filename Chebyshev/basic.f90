module basic
    use libHd
    use libwfn
    implicit none

!Input variable
    real*8 ::mass
    integer::order
    real*8 ::left, right, dq, &
             Hmin, Hmax, &
             kmin

end module basic