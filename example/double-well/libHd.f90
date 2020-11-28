! A double well describing the proton tranfer between two H2O
! The water molecules are staggered and fixed, with O-O distance = 5 Bohr
! The proton is on the line between two oxygens
! The total point group is C2v (D2d when the proton is at centre)
! Electronic structure: B3LYP-GD3BJ / 6-311++G*
! Fit with 6th order polynomial
module libHd
    implicit none

integer, parameter::NStates = 1

contains
subroutine initialize_libHd()
end subroutine initialize_libHd

subroutine compute_Hd(q, H)
    real*8, intent(in)::q
    real*8, dimension(1, 1), intent(out)::H
    real*8::q2, q4
    q2 = q * q
    q4 = q2 * q2
    H(1, 1) = 3d-2 * q2 * (q4 + 2d0 * q2 - 1d0)
end subroutine compute_Hd

end module libHd