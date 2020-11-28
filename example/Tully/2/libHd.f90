! Tully 2nd model: dual avoided crossing
module libHd
    implicit none

contains
subroutine initialize_libHd()
end subroutine initialize_libHd

subroutine compute_Hd(x, H)
    real*8, dimension(3), intent(in):: x
    real*8, dimension(2, 2), intent(out):: H
    real*8, parameter:: A = 0.1, B = 0.28, C = 0.015, D = 0.06, E0 = 0.05
    H(1, 1) = 0d0
    H(2, 1) =  C * exp(-D * x(1)*x(1))
    H(2, 2) = -A * exp(-B * x(1)*x(1)) + E0
end subroutine compute_Hd

subroutine compute_Hd_dHd(x, H, dH)
    real*8, dimension(3), intent(in):: x
    real*8, dimension(2, 2), intent(out):: H
    real*8, dimension(3, 2, 2), intent(out):: dH
    real*8, parameter:: A = 0.1, B = 0.28, C = 0.015, D = 0.06, E0 = 0.05
    ! Hd
    H(1, 1) = 0d0
    H(2, 1) =  C * exp(-D * x(1)*x(1))
    H(2, 2) = -A * exp(-B * x(1)*x(1)) + E0
    ! dHd
    dH(1, 1, 1) = 0d0
    dH(1, 2, 1) = -2d0 * C * D * x(1) * exp(-D * x(1)*x(1))
    dH(1, 2, 2) =  2d0 * A * B * x(1) * exp(-B * x(1)*x(1))
    dH(2:3, :, :) = 0d0
end subroutine compute_Hd_dHd

end module libHd