!3-state dibatization based on Euler angle rotation
!Reference: chaoqun undergraduate thesis
program main
    use Mathematics
    implicit none
    integer::i
    real*8::x, dx
    real*8, dimension(3)::angle
    real*8, dimension(3, 3)::H

    x = -30d0
    dx = 0.01d0
    angle = 1d-12 ! converged compared to 1d-8

    open(unit=99 , file="angle.in", form="unformatted", status="replace")
    open(unit=100, file="energy.txt", status="replace")
    open(unit=101, file="coupling.txt", status="replace")
        write(99)x, angle
        call compute_Hd(x, angle, H)
        write(100,'(4F20.15)')x, H(1,1), H(2,2), H(3,3)
        write(101,'(4F20.15)')x, H(2,1), H(3,1), H(3,2)
        do i = 1, ceiling(60d0 / dx)
            call RK4(x, angle, angle, compute_derivative, dx, 3)
            x = x + dx
            write(99)x, angle
            call compute_Hd(x, angle, H)
            write(100,'(4F20.15)')x, H(1,1), H(2,2), H(3,3)
            write(101,'(4F20.15)')x, H(2,1), H(3,1), H(3,2)
        end do
    close(99); close(100); close(101)

contains
real*8 function d21(x)
    real*8, intent(in)::x
    d21 = 0d0
end function d21

real*8 function d31(x)
    real*8, intent(in)::x
    d31 = 1d0 * x * exp(-x * x / 2d0)
end function d31

real*8 function d32(x)
    real*8, intent(in)::x
    d32 = 1d1 * x * exp(-x * x / 2d0)
end function d32

subroutine compute_derivative(d, x, angle, dim)
    real*8, dimension(3), intent(out)::d
    real*8, intent(in)::x
    real*8, dimension(3), intent(in)::angle
    integer, intent(in)::dim
    d(1) = -d21(x) - d31(x) * sin(angle(1)) * cotan(angle(2)) + d32(x) * cos(angle(1)) * cotan(angle(2))
    d(2) = d31(x) * cos(angle(1)) + d32(x) * sin(angle(1))
    d(3) = d31(x) * sin(angle(1)) / sin(angle(2)) - d32(x) * cos(angle(1)) / sin(angle(2))
end subroutine compute_derivative

subroutine rotation_matrix(angle, U)
    real*8, dimension(3), intent(in)::angle
    real*8, dimension(3, 3), intent(out)::U
    real*8::alpha, beta, gamma
    alpha = angle(1)
    beta  = angle(2)
    gamma = angle(3)
    U(1, 1) =  cos(alpha) * cos(beta ) * cos(gamma) - sin(alpha) * sin(gamma)
    U(2, 1) =  sin(alpha) * cos(beta ) * cos(gamma) + cos(alpha) * sin(gamma)
    U(3, 1) = -sin(beta ) * cos(gamma)
    U(1, 2) = -cos(alpha) * cos(beta ) * sin(gamma) - sin(alpha) * cos(gamma)
    U(2, 2) = -sin(alpha) * cos(beta ) * sin(gamma) + cos(alpha) * cos(gamma)
    U(3, 2) =  sin(beta ) * sin(gamma)
    U(1, 3) =  cos(alpha) * sin(beta )
    U(2, 3) =  sin(alpha) * sin(beta )
    U(3, 3) =  cos(beta )
end subroutine rotation_matrix

subroutine compute_Hd(x, angle, H)
    real*8, intent(in)::x
    real*8, dimension(3), intent(in)::angle
    real*8, dimension(3, 3), intent(out)::H
    real*8, dimension(3, 3)::U
    H = 0d0
    H(1, 1) = -0.005d0
    H(3, 3) =  0.005d0
    call rotation_matrix(angle, U)
    H = matmul(transpose(U), matmul(H, U))
end subroutine compute_Hd

end program main