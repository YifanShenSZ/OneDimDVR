!My super-exchange model modified from linjun's
!His diabatic PES becomes my adiabatic
!His diabatic coupling becomes my NAC
!My diabatization is a 3-state dibatization based on Euler angle rotation
!Reference: chaoqun undergraduate thesis
module libHd
    implicit none

integer, parameter::NStates = 3

!Global variable
    real*8, dimension(10001)::libHd_xs
    real*8, dimension(3, 10001)::libHd_angles

contains
subroutine initialize_libHd()
    integer::i
    open(unit=99, file="angle.in", form="unformatted", status="old")
    do i = 1, size(libHd_xs)
        read(99)libHd_xs(i), libHd_angles(:, i)
    end do
    close(99)
end subroutine initialize_libHd

!Use linear interpolation to calculate the angles at a specific q
!Calculate the rotation matrix from Euler angles
!Reference: chaoqun undergraduate thesis
subroutine rotation_matrix(q, R)
    real*8, intent(in)::q
    real*8, dimension(3, 3), intent(out)::R
    integer::location
    real*8, dimension(3)::angle
    real*8::alpha, beta, gamma
    !Use linear interpolation to calculate the angles at a specific q
    location = floor(((q + 50d0) / 0.01d0))
    if (1 <= location .and. location <= size(libHd_xs)) then
        angle = libHd_angles(:, location + 1) &
              + (libHd_angles(:, location + 2) - libHd_angles(:, location + 1)) &
              / (libHd_xs(location + 2) - libHd_xs(location + 1)) &
              * (q + 50d0 - floor(((q + 50d0) / 0.01d0)) * 0.01d0)
    else
        angle = 0d0
    end if
    alpha = angle(1)
    beta  = angle(2)
    gamma = angle(3)
    !Calculate the rotation matrix from Euler angles
    R(1, 1) =  cos(alpha) * cos(beta ) * cos(gamma) - sin(alpha) * sin(gamma)
    R(2, 1) =  sin(alpha) * cos(beta ) * cos(gamma) + cos(alpha) * sin(gamma)
    R(3, 1) = -sin(beta ) * cos(gamma)
    R(1, 2) = -cos(alpha) * cos(beta ) * sin(gamma) - sin(alpha) * cos(gamma)
    R(2, 2) = -sin(alpha) * cos(beta ) * sin(gamma) + cos(alpha) * cos(gamma)
    R(3, 2) =  sin(beta ) * sin(gamma)
    R(1, 3) =  cos(alpha) * sin(beta )
    R(2, 3) =  sin(alpha) * sin(beta )
    R(3, 3) =  cos(beta )
end subroutine rotation_matrix

subroutine compute_Hd(q, H)
    real*8, intent(in)::q
    real*8, dimension(3, 3), intent(out)::H
    real*8, dimension(3, 3)::R
    !Compute adiabatic H
    H = 0d0
    H(1, 1) = -0.005d0
    H(3, 3) =  0.005d0
    !Rotate from adiabatz to diabatz
    call rotation_matrix(q, R)
    H = matmul(R, matmul(H, transpose(R)))
end subroutine compute_Hd

end module libHd