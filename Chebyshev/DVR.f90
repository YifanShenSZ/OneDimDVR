module DVR
    use basic
    implicit none

contains
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

!Compute edge wave function damping
!H. Guo 2006 Phys. Rev. A
subroutine compute_damping(grids, D, NGrids)
    integer, intent(in)::NGrids
    real*8, dimension(NGrids), intent(in)::grids
    complex*16, dimension(NGrids), intent(out)::D
    integer::i
    real*8::sigma_left, sigma_right, y
    !The parameters in H. Guo 2006 Phys. Rev. A suggests
    !e^-0.09 or e^-0.08 at boundary
    sigma_left = (left - grids(1)) / 0.3d0
    sigma_right = (grids(NGrids) - right) / 0.3d0
    do i = 1, NGrids
        if(grids(i) < left) then
            y = (grids(i) - left) / sigma_left
            D(i) = exp(-y * y)
        else if(grids(i) > right) then
            y = (grids(i) - right) / sigma_right
            D(i) = exp(-y * y)
        else
            D(i) = (1d0,0d0)
        end if
    end do
end subroutine compute_damping

end module DVR