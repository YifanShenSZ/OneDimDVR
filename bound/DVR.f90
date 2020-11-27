module DVR
    use libHd
    implicit none

    !Input variable required in this module
    real*8::mass

contains
!Compute DVR kinetic energy matrix T
!W. H. Miller 1992 J. Chem. Phys.
subroutine compute_kinetic(grids, T, NGrids)
    integer, intent(in)::NGrids
    real*8, dimension(NGrids), intent(in)::grids
    real*8, dimension(NGrids, NGrids), intent(out)::T
    integer::i, j, k
    real*8::dq
    do i = 1, NGrids
        T(i, i) = 3.2898681336964524d0 ! pi^2 / 3
        do j = i + 1, NGrids
            k = j - i
            T(j, i) = 2d0 / (k * k)
            if (mod(k, 2) == 1) then
                T(j, i) = -T(j, i)
            end if
        end do
    end do
    dq = grids(2) - grids(1)
    T = T / 2d0 / mass / dq / dq
    forall (i = 1 : NGrids, j = 1 : NGrids, j > i)
        T(i, j) = T(j, i)
    end forall
end subroutine compute_kinetic

!Compute DVR Hamiltonian matrix
subroutine compute_Hamiltonian(grids, H, NGrids, NStates)
    integer, intent(in)::NGrids, NStates
    real*8, dimension(NGrids), intent(in)::grids
    complex*16, dimension(NGrids, NGrids, NStates, NStates), intent(out)::H
    integer::i, j
    real*8, dimension(NStates, NStates)::Hd
    real*8, dimension(NGrids, NGrids)::T
    !Fill kinetic energy into Hamiltonian matrix
    call compute_kinetic(grids, T, NGrids)
    forall (i = 1 : NStates, j = 1 : NStates, j >= i)
        H(:, :, j, i) = T
    end forall
    !Add potential energy to Hamiltonian matrix
    do i = 1, NGrids
        call compute_Hd(grids(i), Hd)
        H(i, i, :, :) = H(i, i, :, :) + Hd
    end do
end subroutine compute_Hamiltonian

end module DVR