module DVR
    use basic
    implicit none

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

!Compute DVR Hamiltonian matrix with absorbing potential
!D. E. Manolopoulos 2002 J. Chem. Phys.
subroutine compute_AbsorbedHamiltonian(grids, H, NGrids, NStates)
    integer, intent(in)::NGrids, NStates
    real*8, dimension(NGrids), intent(in)::grids
    complex*16, dimension(NGrids, NStates, NGrids, NStates), intent(out)::H
    integer::i, j
    real*8, dimension(NGrids, NGrids)::T
    real*8, dimension(NStates, NStates)::Hd
    !Absorbing potential
    real*8::c1, y, denominator
    complex*16::c2, absorbance
    !c is a const Manolopoulos obtained from elliptic integral
    real*8,parameter::csq = 6.87519864356d0, &
                      inv_csq = 0.14545034286924935d0 ! 1 / csq
    !Fill kinetic energy into Hamiltonian matrix
    call compute_kinetic(grids, T, NGrids)
    forall (i = 1 : NStates)
        H(:, i, :, i) = T
    end forall
    forall (i = 1 : NStates, j = 1 : NStates, j > i)
        H(:, j, :, i) = 0d0
    end forall
    !Add potential energy to Hamiltonian matrix
    c1 = csq * kmin * kmin / 39.47841760435743d0 ! 4 pi^2
    c2 = (0d0, -4d0) * kmin * kmin / mass
    do i = 1, NGrids
        call compute_Hd(grids(i), Hd)
        H(i, :, i, :) = H(i, :, i, :) + Hd
        !Add absorbance
        if(grids(i) > right) then
            y = grids(i) - right
            y = y * y
            y = c1 * y
            denominator = csq - y
            denominator = denominator * denominator
            absorbance = c2 * ((csq + y) / denominator - inv_csq)
            forall (j = 1 : NStates)
                H(i, j, i, j) = H(i, j, i, j) + absorbance
            end forall
        else if(grids(i) < left) then
            y = grids(i) - left
            y = y * y
            y = c1 * y
            denominator = csq - y
            denominator = denominator * denominator
            absorbance = c2 * ((csq + y) / denominator - inv_csq)
            forall (j = 1 : NStates)
                H(i, j, i, j) = H(i, j, i, j) + absorbance
            end forall
        end if
    end do
end subroutine compute_AbsorbedHamiltonian

!Compute DVR Hamiltonian matrix without/with absorbing potential
!D. E. Manolopoulos 2002 J. Chem. Phys.
subroutine compute_Hamiltonian_AbsorbedHamiltonian(grids, H, H_ab, NGrids, NStates)
    integer, intent(in)::NGrids, NStates
    real*8, dimension(NGrids), intent(in)::grids
    complex*16, dimension(NGrids, NStates, NGrids, NStates), intent(out)::H, H_ab
    integer::i, j
    real*8, dimension(NGrids, NGrids)::T
    real*8, dimension(NStates, NStates)::Hd
    !Absorbing potential
    real*8::c1, y, denominator
    complex*16::c2, absorbance
    !c is a const Manolopoulos obtained from elliptic integral
    real*8,parameter::csq = 6.87519864356d0, &
                      inv_csq = 0.14545034286924935d0 ! 1 / csq
    !Fill kinetic energy into Hamiltonian matrix
    call compute_kinetic(grids, T, NGrids)
    forall (i = 1 : NStates)
        H(:, i, :, i) = T
        H_ab(:, i, :, i) = T
    end forall
    forall (i = 1 : NStates, j = 1 : NStates, j > i)
        H(:, j, :, i) = 0d0
        H_ab(:, j, :, i) = 0d0
    end forall
    !Add potential energy to Hamiltonian matrix
    c1 = csq * kmin * kmin / 39.47841760435743d0 ! 4 pi^2
    c2 = (0d0, -4d0) * kmin * kmin / mass
    do i = 1, NGrids
        call compute_Hd(grids(i), Hd)
        H(i, :, i, :) = H(i, :, i, :) + Hd
        H_ab(i, :, i, :) = H_ab(i, :, i, :) + Hd
        !Add absorbance
        if(grids(i) > right) then
            y = grids(i) - right
            y = y * y
            y = c1 * y
            denominator = csq - y
            denominator = denominator * denominator
            absorbance = c2 * ((csq + y) / denominator - inv_csq)
            forall (j = 1 : NStates)
                H_ab(i, j, i, j) = H_ab(i, j, i, j) + absorbance
            end forall
        else if(grids(i) < left) then
            y = grids(i) - left
            y = y * y
            y = c1 * y
            denominator = csq - y
            denominator = denominator * denominator
            absorbance = c2 * ((csq + y) / denominator - inv_csq)
            forall (j = 1 : NStates)
                H_ab(i, j, i, j) = H_ab(i, j, i, j) + absorbance
            end forall
        end if
    end do
end subroutine compute_Hamiltonian_AbsorbedHamiltonian

end module DVR