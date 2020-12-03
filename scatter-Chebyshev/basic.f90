module basic
    use libHd
    use libwfn
    implicit none

!Input variable
    real*8 ::mass, &
             total_time, output_interval, &
             left, right, dq, &
             order, Hmin, Hmax

end module basic