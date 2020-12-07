module DVR
    use basic
    implicit none

contains
!Compute DVR momentum matrix p
!Derived by infinite order finite difference
subroutine compute_momentum(dq, p, NGrids)
    integer, intent(in)::NGrids
    real*8, intent(in)::dq
    complex*16, dimension(NGrids, NGrids), intent(out)::p
    integer::i, j, k
    do i = 1, NGrids
        p(i, i) = 0d0
        do j = i + 1, NGrids
            k = j - i
            p(j, i) = 1d0 / k
            if (mod(k, 2) == 1) then
                p(j, i) = -p(j, i)
            end if
        end do
    end do
    p = p * (0d0,-1d0) / dq
    forall (i = 1 : NGrids, j = 1 : NGrids, i < j)
        p(i, j) = conjg(p(j, i))
    end forall
end subroutine compute_momentum

!Compute DVR kinetic energy matrix T
!W. H. Miller 1992 J. Chem. Phys.
subroutine compute_kinetic(dq, T, NGrids)
    integer, intent(in)::NGrids
    real*8, intent(in)::dq
    complex*16, dimension(NGrids, NGrids), intent(out)::T
    integer::i, j, k
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
    T = T / 2d0 / mass / dq / dq
end subroutine compute_kinetic

!Compute potential energy matrix V with absorbing potential
!D. E. Manolopoulos 2004 J. Chem. Phys.
subroutine compute_potential(grids, V, NGrids, NStates)
    integer, intent(in)::NGrids, NStates
    real*8, dimension(NGrids), intent(in)::grids
    complex*16, dimension(NGrids, NStates, NStates), intent(out)::V
    integer::i, j
    real*8, dimension(NStates, NStates)::Hd
    !Absorbing potential
    real*8::c1, y, denominator
    complex*16::c2, absorbance
    !c is a const Manolopoulos obtained from elliptic integral
    real*8, parameter::csq = 6.87519864356d0, &
                       inv_csq = 0.14545034286924935d0 ! 1 / csq
    c1 = csq * kmin * kmin / 39.47841760435743d0 ! 4 pi^2
    c2 = (0d0, -4d0) * kmin * kmin / mass
    do i = 1, NGrids
        call compute_Hd(grids(i), Hd)
        V(i, :, :) = Hd
        !Add absorbance
        if(grids(i) > right) then
            y = grids(i) - right
            y = y * y
            y = c1 * y
            denominator = csq - y
            denominator = denominator * denominator
            absorbance = c2 * ((csq + y) / denominator - inv_csq)
            forall (j = 1 : NStates)
                V(i, j, j) =  V(i, j, j) + absorbance
            end forall
        else if(grids(i) < left) then
            y = grids(i) - left
            y = y * y
            y = c1 * y
            denominator = csq - y
            denominator = denominator * denominator
            absorbance = c2 * ((csq + y) / denominator - inv_csq)
            forall (j = 1 : NStates)
                V(i, j, j) =  V(i, j, j) + absorbance
            end forall
        end if
    end do
end subroutine compute_potential

end module DVR