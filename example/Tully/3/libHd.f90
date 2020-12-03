! Tully 3rd model: extended coupling with reflection
module libHd
    implicit none

integer, parameter::NStates = 2

contains
subroutine initialize_libHd()
end subroutine initialize_libHd

subroutine compute_Hd(q, H)
    real*8, intent(in)::q
    real*8, dimension(2, 2), intent(out)::H
    real*8, parameter::A = 6d-4, B = 0.1, C = 0.9
    H(1, 1) = A
    if (q < 0d0) then
        H(2, 1) = B * exp(C * q)
    else
        H(2, 1) = B * (2d0 - exp(-C * q))
    end if
    H(2, 2) = -A
end subroutine compute_Hd

end module libHd