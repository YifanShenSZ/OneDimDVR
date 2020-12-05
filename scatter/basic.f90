module basic
    use libHd
    use libwfn
    implicit none

!Input variable
    real*8 ::mass, &
             total_time, dt, output_interval, &
             left, right, dq, &
             kmin, kmax

end module basic