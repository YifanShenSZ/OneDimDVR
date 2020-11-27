module DVR
    use libHd
    use Basic
    implicit none

!Parameter
    real*8::maxpopdev =1d-6, &!Stop wave function evolution when population no longer equals to 1
            minpop    =1d-2   !Stop transmission & reflection calculation when all population are absorbed
                              !Note: absorbing potential can only absorb 99%, 1% will be reflected by boundary

contains
!DVR Hamiltonian + absorbing potential
!D. E. Manolopoulos 2002 J. Chem. Phys.
subroutine scatter_Hamiltonian(grids, H, NGrids, NStates)
    integer, intent(in)::NGrids, NStates
    real*8, dimension(NGrids), intent(in)::grids
    complex*16, dimension(NGrids, NGrids, NStates, NStates), intent(out)::H
    integer::i, j, k
    real*8::dq, y, absorbance, denominator
    real*8, dimension(NStates, NStates)::Hd
    real*8, dimension(NGrids, NGrids)::T
    real*8, parameter::pipid3 = 3.2898681336964524d0 ! pi^2 / 3
    real*8, parameter::pipim4 = 39.47841760435743d0  ! 4 pi^2
    !A const Manolopoulos obtained from elliptic integral
    real*8,parameter::csq = 6.87519864356d0, &
                      inv_csq = 0.14545034286924935d0 ! 1 / csq
    !Compute kinetic energy matrix T
    do i = 1, NGrids
        T(i, i) = pipid3
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
    !Fill kinetic energy into Hamiltonian matrix
    forall (i = 1 : NStates, j = 1 : NStates, j > i)
        H(:, :, j, i) = T
    end forall
    !Add potential energy to Hamiltonian matrix
    do i = 1, NGrids
        call compute_Hd(grids(i), Hd)
        H(i, i, :, :) = H(i, i, :, :) + Hd
        !Add absorbance
        if(grids(i) > right) then
            y = csq * kmins / pipim4 * (grids(i)-right)**2
            denominator = csq - y
            denominator = denominator * denominator
            absorbance = -ci * kmins * 4d0 / mass * ((csq + y) / denominator - inv_csq)
            forall (j = 1 : NStates)
                H(i, i, j, j) = H(i, i, j, j) + absorbance
            end forall
        else if(grids(i) < left) then
            y = csq * kmins / pipim4 * (grids(i)-left)**2
            denominator = csq - y
            denominator = denominator * denominator
            absorbance = -ci * kmins * 4d0 / mass * ((csq + y) / denominator - inv_csq)
            forall (j = 1 : NStates)
                H(i, i, j, j) = H(i, i, j, j) + absorbance
            end forall
        end if
    end do
end subroutine scatter_Hamiltonian

end module DVR