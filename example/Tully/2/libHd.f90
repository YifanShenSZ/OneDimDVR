! Tully 2nd model: dual avoided crossing
module libHd
    implicit none

integer, parameter::NStates = 2

contains
subroutine initialize_libHd()
end subroutine initialize_libHd

subroutine compute_Hd(q, H)
    real*8, intent(in)::q
    real*8, dimension(2, 2), intent(out)::H
    real*8, parameter::A = 0.1, B = 0.28, C = 0.015, D = 0.06, E0 = 0.05
    H(1, 1) = 0d0
    H(2, 1) =  C * exp(-D * q *q)
    H(2, 2) = -A * exp(-B * q *q) + E0
end subroutine compute_Hd

end module libHd