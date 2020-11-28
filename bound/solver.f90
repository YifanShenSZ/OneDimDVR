module solver
    use LinearAlgebra
    use basic
    use DVR
    implicit none

contains
!This is the most basic solver:
!    1. Build Hamiltonian matrix
!    2. Diagonalize Hamiltonian matrix
!    3. Propagate by exp(-i E) * H eigen vector
subroutine propagate_wavefunction()
    !Grid points
    integer::NSnapshots, NGrids
    real*8, allocatable, dimension(:)::grids
    !DVR Hamiltonian
    integer::NTotal
    real*8, allocatable, dimension(:)::energy
    complex*16, allocatable, dimension(:, :, :, :)::H
    !DVR wave function
    complex*16, allocatable, dimension(:)::phase, adiabatic_wfn
    complex*16, allocatable, dimension(:, :)::diabatic_wfn
    !Work variable
    integer::i, j, k
    !Discretize time and space
    NSnapshots = floor(total_time / output_interval) + 1
    NGrids = floor((right - left) / dq) + 1
    dq = (right - left) / (NGrids - 1)
    allocate(grids(NGrids))
    grids(1) = left
    do i = 2, NGrids
        grids(i) = grids(i - 1) + dq
    end do
    open(unit=99, file="grids.out", form="unformatted", status="replace")
        write(99)grids
    close(99)
    !Build and diagonalize Hamiltonian
    NTotal = NGrids * NStates
    allocate(energy(NTotal))
    allocate(H(NGrids, NStates, NGrids, NStates))
    call compute_Hamiltonian(grids, H, NGrids, NStates)
    call My_zheev('V', H, energy, NTotal)
    !Propagate wave function
    open(unit=99, file="wfn.out", form="unformatted", status="replace")
        !Initial condition
        allocate(diabatic_wfn(NGrids, NStates))
        do i = 1, NGrids
            call init_wfn(grids(i), diabatic_wfn(i, :))
        end do
        write(99)diabatic_wfn
        !Propagate by exp(-i E) * H eigen vector
        allocate(phase(NTotal))
        forall(i = 1 : NTotal)
            phase(i) = exp((0d0,-1d0) * energy(i) * output_interval)
        end forall
        allocate(adiabatic_wfn(NTotal))
        !adiabatic_wfn = H^dagger . diabatic_wfn
        call zgemv('C', NTotal, NTotal, (1d0,0d0), H, NTotal, diabatic_wfn, 1, (0d0,0d0), adiabatic_wfn, 1)
        do i = 2, NSnapshots
            adiabatic_wfn = adiabatic_wfn * phase
            !diabatic_wfn = H . adiabatic_wfn
            call zgemv('N', NTotal, NTotal, (1d0,0d0), H, NTotal, adiabatic_wfn, 1, (0d0,0d0), diabatic_wfn, 1)
            write(99)diabatic_wfn
        end do
    close(99)
    !Leave a checkpoint
    open(unit=99, file="checkpoint.out", status="replace")
        write(99,*)"Number of electronic states:"
        write(99,*)NStates
        write(99,*)"Number of time snapshots:"
        write(99,*)NSnapshots
        write(99,*)"Number of grid points:"
        write(99,*)NGrids
    close(99)
    !Clean up
    deallocate(grids)
    deallocate(energy)
    deallocate(H)
    deallocate(phase)
    deallocate(adiabatic_wfn)
    deallocate(diabatic_wfn)
end subroutine propagate_wavefunction

end module solver