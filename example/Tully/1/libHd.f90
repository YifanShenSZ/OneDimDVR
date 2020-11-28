! Tully 1st model: simple avoided crossing
module libHd
    implicit none

integer, parameter::NStates = 2

contains
subroutine initialize_libHd()
end subroutine initialize_libHd

subroutine compute_Hd(q, H)
    real*8, intent(in):: q
    real*8, dimension(2, 2), intent(out):: H
    real*8, parameter:: A = 0.01, B = 1.6, C = 0.005, D = 1.0
    if (q > 0d0) then
        H(1, 1) =  A * (1d0 - exp(-B * q))
    else
        H(1, 1) = -A * (1d0 - exp( B * q))
    end if
    H(2, 2) = -H(1, 1)
    H(2, 1) = C * exp(-D * q * q)
end subroutine compute_Hd

end module libHd