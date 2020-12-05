module DVR
    use basic
    implicit none

contains
!Compute edge wave function damping
!H. Guo 2006 Phys. Rev. A
subroutine compute_damping(grids, D, NGrids)
    integer, intent(in)::NGrids
    real*8, dimension(NGrids), intent(in)::grids
    complex*16, dimension(NGrids), intent(out)::D
    integer::i
    real*8::sigma, y
    !This parameterization is borrowed from D. E. Manolopoulos 2002 J. Chem. Phys.
    !The length of absorption region = 2pi / kmin
    !And we want sufficiently small damping at the boundary (e^-16 < 1e-6)
    sigma = 6.283185307179586d0 / kmin / 4d0
    do i = 1, NGrids
        if(grids(i) > right) then
            y = (grids(i) - right) / sigma
            D(i) = exp(-y * y)
        else if(grids(i) < left) then
            y = (grids(i) - left) / sigma
            D(i) = exp(-y * y)
        else
            D(i) = (1d0,0d0)
        end if
    end do
end subroutine compute_damping

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

!Compute potential energy matrix V
subroutine compute_potential(grids, V, NGrids, NStates)
    integer, intent(in)::NGrids, NStates
    real*8, dimension(NGrids), intent(in)::grids
    complex*16, dimension(NGrids, NStates, NStates), intent(out)::V
    integer::i
    real*8, dimension(NStates, NStates)::Hd
    do i = 1, NGrids
        call compute_Hd(grids(i), Hd)
        V(i, :, :) = Hd
    end do
end subroutine compute_potential

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
end subroutine compute_momentum

end module DVR