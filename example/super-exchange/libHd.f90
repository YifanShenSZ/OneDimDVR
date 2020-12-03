!My super-exchange model modified from linjun's
!His diabatic PES becomes my adiabatic
!His diabatic coupling becomes my NAC
!My diabatization is a 3-state dibatization based on Euler angle rotation
!Reference: chaoqun undergraduate thesis
module libHd
    implicit none

integer, parameter::NStates = 3

!Global variable
    real*8, dimension(6001)::libHd_xs
    real*8, dimension(3, 6001)::libHd_angles

contains
subroutine initialize_libHd()
    integer::i
    open(unit=99, file="angle.in", form="unformatted")
    do i = 1, 6001
        read(99)libHd_xs(i), libHd_angles(:, i)
    end do
    close(99)
end subroutine initialize_libHd

!Use linear interpolation to calculate the angles at a specific x
!Calculate the rotation matrix from Euler angles
!Reference: chaoqun undergraduate thesis
subroutine rotation_matrix(x, U)
    real*8, intent(in)::x
    real*8, dimension(3, 3), intent(out)::U
    integer::location
    real*8, dimension(3)::angle
    real*8::alpha, beta, gamma
    !Use linear interpolation to calculate the angles at a specific x
    location = floor(((x + 30d0) / 0.01d0))
    if (1 <= location .and. location <= 6001) then
        angle = libHd_angles(:, location + 1) &
              + (libHd_angles(:, location + 2) - libHd_angles(:, location + 1)) &
              / (libHd_xs(location + 2) - libHd_xs(location + 1)) &
              * (x + 30d0 - floor(((x + 30d0) / 0.01d0)) * 0.01d0)
    else
        angle = 0d0
    end if
    alpha = angle(1)
    beta  = angle(2)
    gamma = angle(3)
    !Calculate the rotation matrix from Euler angles
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

subroutine compute_Hd(x, H)
    real*8, intent(in)::x
    real*8, dimension(3, 3), intent(out)::H
    real*8, dimension(3, 3)::U
    H = 0d0
    H(1, 1) = -0.005d0
    H(3, 3) =  0.005d0
    call rotation_matrix(x, U)
    H = matmul(transpose(U), matmul(H, U))
end subroutine compute_Hd

end module libHd