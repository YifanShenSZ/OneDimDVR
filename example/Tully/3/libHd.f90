! Tully 3rd model: extended coupling with reflection
module libHd
    implicit none

contains
subroutine initialize_libHd()
end subroutine initialize_libHd

subroutine compute_Hd(x, H)
    real*8, dimension(3), intent(in):: x
    real*8, dimension(2, 2), intent(out):: H
    real*8, parameter:: A = 6d-4, B = 0.1, C = 0.9
    H(1, 1) = A
    if (x(1) < 0d0) then
        H(2, 1) =  B * exp( C * x(1))
    else
        H(2, 1) = -B * exp(-C * x(1)) + 2d0 * B
    end if
    H(2, 2) = -A
end subroutine compute_Hd

subroutine compute_Hd_dHd(x, H, dH)
    real*8, dimension(3), intent(in):: x
    real*8, dimension(2, 2), intent(out):: H
    real*8, dimension(3, 2, 2), intent(out):: dH
    real*8, parameter:: A = 6d-4, B = 0.1, C = 0.9
    ! Hd
    H(1, 1) = A
    if (x(1) < 0d0) then
        H(2, 1) =  B * exp( C * x(1))
    else
        H(2, 1) = -B * exp(-C * x(1)) + 2d0 * B
    end if
    H(2, 2) = -A
    ! dHd
    dH(1, 1, 1) = 0d0
    if (x(1) < 0d0) then
        dH(1, 2, 1) = B * C * exp( C * x(1))
    else
        dH(1, 2, 1) = B * C * exp(-C * x(1))
    end if
    dH(1, 2, 2) = 0d0
    dH(2:3, :, :) = 0d0
end subroutine compute_Hd_dHd

end module libHd