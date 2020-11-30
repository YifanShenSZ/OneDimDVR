module solver
    use LinearAlgebra
    use basic
    use DVR
    implicit none

!Parameter
    logical::auto_spacing = .true.!Determine grid spacing by absorbance condition
    real*8::maxpopdev = 1d-6, &!Stop wave function propagation when population no longer equals to 1
            minpop    = 1d-4   !Stop transmission & reflection calculation when all population are absorbed
                               !Note that absorbing potential can only absorb 99%, 1% will be reflected by boundary,
                               !so a round trip through absorbing region absorbs 99.99%

contains
subroutine RK4(H, old, new, dt, dim)
    integer,intent(in)::dim
    complex*16, dimension(dim, dim), intent(in)::H
    complex*16, dimension(dim), intent(in)::old
    complex*16, dimension(dim), intent(out)::new
    real*8, intent(in)::dt
    real*8::dtd2
    complex*16, dimension(dim)::k1, k2, k3, k4
    dtd2 = dt / 2d0
    call TDSE(k1, old            , dim)
    call TDSE(k2, old + k1 * dtd2, dim)
    call TDSE(k3, old + k2 * dtd2, dim)
    call TDSE(k4, old + k3 * dt  , dim)
    new = old + dt / 6d0 * (k1 + 2d0 * k2 + 2d0 * k3 + k4)
    contains
    !Time-dependent Schrodinger equation
    subroutine TDSE(gradient, wfn, dim)
        integer, intent(in)::dim
        complex*16, dimension(dim), intent(inout)::gradient
        complex*16, dimension(dim), intent(in)::wfn
        call zsymv('L', dim, (0d0,-1d0), H, dim, wfn, 1, (0d0,0d0), gradient, 1)
    end subroutine TDSE
end subroutine RK4

!The basic solver cannot be used since the absorbing Hamiltonian is no longer Hermitian
!This is an RK4 solver:
!    1. Build absorbing Hamiltonian matrix
!    2. Propagate by RK4
subroutine propagate_wavefunction()
    !Grid points
    integer::OutputStep, NSnapshots, NUsualGrids, NAbsorbGrids, NGrids
    real*8, allocatable, dimension(:)::snapshots, grids
    !DVR Hamiltonian with absorbing potential
    integer::NTotal
    complex*16, allocatable, dimension(:, :, :, :)::H
    !DVR wave function
    complex*16, allocatable, dimension(:, :)::wfn
    !Work variable
    integer::i
    !Discretize time
    OutputStep = floor(output_interval / dt)
    dt = output_interval / OutputStep
    NSnapshots = floor(total_time / output_interval) + 1
    !Discretize space
    if (auto_spacing) then
        !D. E. Manolopoulos 2002 J. Chem. Phys. suggests at least 5 grid points every 2 pi / kmax
        dq = min(6.283185307179586d0 / kmax / 5d0, dq)
        write(*,*)"Consider absorbance condition and user input, grid spacing is set to ", dq
    end if
    NUsualGrids = floor((right - left) / dq) + 1
    dq = (right - left) / (NUsualGrids - 1)
    NAbsorbGrids = floor(6.283185307179586d0 / kmin / dq) - 1
    NGrids = NUsualGrids + 2 * NAbsorbGrids
    allocate(grids(NGrids))
    grids(NAbsorbGrids) = left - dq
    grids(NAbsorbGrids + NUsualGrids + 1) = right + dq
    do i = 1, NAbsorbGrids - 1
        grids(NAbsorbGrids - i) = grids(NAbsorbGrids - i + 1) - dq
        grids(NAbsorbGrids + NUsualGrids + i + 1) = grids(NAbsorbGrids + NUsualGrids + i) + dq
    end do
    grids(NAbsorbGrids + 1) = left
    do i = 2, NUsualGrids
        grids(NAbsorbGrids + i) = grids(NAbsorbGrids + i - 1) + dq
    end do
    !Build absorbing Hamiltonian
    NTotal = NGrids * NStates
    allocate(H(NGrids, NStates, NGrids, NStates))
    call compute_AbsorbedHamiltonian(grids, H, NGrids, NStates)
    !Propagate wave function
    open(unit=99, file="wfn.out", form="unformatted", status="replace")
        !Initial condition
        allocate(wfn(NGrids, NStates))
        do i = 1, NGrids
            call init_wfn(grids(i), wfn(i, :))
        end do
        write(99)wfn(NAbsorbGrids + 1 : NAbsorbGrids + NUsualGrids, :)
        !Propagate by RK4
        do i = 1, (NSnapshots - 1) * OutputStep
            call RK4(H, wfn, wfn, dt, NTotal)
            !Write trajectory
            if(mod(i, OutputStep) == 0) then
                write(99)wfn(NAbsorbGrids + 1 : NAbsorbGrids + NUsualGrids, :)
            end if
            !Check whether the population has been absorbed too much
            if (abs(population() - 1d0) > maxpopdev) then
                write(*,*)"Stop propagation at time ", i * dt
                write(*,*)"Because population deviates from 1 by more than ", maxpopdev
                write(*,*)"Current population = ", population()
                write(*,*)"Possible reasons are:"
                write(*,*)"1. Grid spacing is too large to integrate accurately"
                write(*,*)"2. Some population has escaped the interaction region"
                NSnapshots = 1 + i / OutputStep
                exit
            end if
        end do
    close(99)
    !Leave a checkpoint
    open(unit=99, file="checkpoint.out", status="replace")
        write(99,*)"Number of electronic states:"
        write(99,*)NStates
        write(99,*)"Number of time snapshots:"
        write(99,*)NSnapshots
        write(99,*)"Number of grid points:"
        write(99,*)NUsualGrids
    close(99)
    !Output the grids used in calculation
    allocate(snapshots(NSnapshots))
    snapshots(1) = 0d0
    do i = 2, NSnapshots
        snapshots(i) = snapshots(i - 1) + output_interval
    end do
    open(unit=99, file="snapshots.out", form="unformatted", status="replace")
        write(99)snapshots
    close(99)
    open(unit=99, file="grids.out", form="unformatted", status="replace")
        write(99)grids(NAbsorbGrids + 1 : NAbsorbGrids + NUsualGrids)
    close(99)
    !Clean up
    deallocate(snapshots)
    deallocate(grids)
    deallocate(H)
    deallocate(wfn)
    contains
    real*8 function population()
        integer::i
        population = 0d0
        do i = 1, NStates
            population = population + dq * dot_product( &
                         wfn(NAbsorbGrids + 1 : NAbsorbGrids + NUsualGrids, i), &
                         wfn(NAbsorbGrids + 1 : NAbsorbGrids + NUsualGrids, i))
        end do
    end function population
end subroutine propagate_wavefunction

!The basic solver cannot be used since the absorbing Hamiltonian is no longer Hermitian
!This is an RK4 solver:
!    1. Build absorbing Hamiltonian matrix
!    2. Propagate by RK4
!The difference to subroutine propagate_wavefunction
!is this subroutine cares only about the absorption at the boundary,
!i.e. transmission and reflection rather than saving the wavefunction
subroutine transmit_reflect()
    !grid points
    integer::NSnapshots, NUsualGrids, NAbsorbGrids, NGrids
    real*8, allocatable, dimension(:)::grids
    !DVR Hamiltonian without/with absorbing potential
    integer::NTotal
    complex*16, allocatable, dimension(:, :, :, :)::H, H_ab
    !DVR wave function
    complex*16, allocatable, dimension(:, :)::wfn, wfn_ab
    !transmission and reflection
    real*8, dimension(NStates)::transmission, reflection
    !work variable
    integer::i, j
    !Discretize time
    NSnapshots = floor(total_time / dt) + 1
    !Discretize space
    if (auto_spacing) then
        !D. E. Manolopoulos 2002 J. Chem. Phys. suggests at least 5 grid points every 2 pi / kmax
        dq = min(6.283185307179586d0 / kmax / 5d0, dq)
        write(*,*)"Consider absorbance condition and user input, grid spacing is set to ", dq
    end if
    NUsualGrids = floor((right - left) / dq) + 1
    dq = (right - left) / (NUsualGrids - 1)
    NAbsorbGrids = floor(6.283185307179586d0 / kmin / dq) - 1
    NGrids = NUsualGrids + 2 * NAbsorbGrids
    allocate(grids(NGrids))
    grids(NAbsorbGrids) = left - dq
    grids(NAbsorbGrids + NUsualGrids + 1) = right + dq
    do i = 1, NAbsorbGrids - 1
        grids(NAbsorbGrids - i) = grids(NAbsorbGrids - i + 1) - dq
        grids(NAbsorbGrids + NUsualGrids + i + 1) = grids(NAbsorbGrids + NUsualGrids + i) + dq
    end do
    grids(NAbsorbGrids + 1) = left
    do i = 2, NUsualGrids
        grids(NAbsorbGrids + i) = grids(NAbsorbGrids + i - 1) + dq
    end do
    !Build absorbing Hamiltonian
    NTotal = NGrids * NStates
    allocate(H(NGrids, NStates, NGrids, NStates))
    allocate(H_ab(NGrids, NStates, NGrids, NStates))
    call compute_Hamiltonian_AbsorbedHamiltonian(grids, H, H_ab, NGrids, NStates)
    !Set initial condition
    allocate(wfn(NGrids, NStates))
    allocate(wfn_ab(NGrids, NStates))
    do i = 1, NGrids
        call init_wfn(grids(i), wfn_ab(i, :))
    end do
    transmission = 0d0
    reflection   = 0d0
    minpop = minpop / dq
    !Propagate wave function
    do i = 1, NSnapshots
        call RK4(H   , wfn_ab, wfn   , dt, NTotal)
        call RK4(H_ab, wfn_ab, wfn_ab, dt, NTotal)
        !Calculate absorption
        forall (j = 1 : NStates)
            transmission(j) = transmission(j) &
                            + dq * dot_product( &
                              wfn(NAbsorbGrids + NUsualGrids + 1 : NGrids, j), &
                              wfn(NAbsorbGrids + NUsualGrids + 1 : NGrids, j)) &
                            - dq * dot_product( &
                              wfn_ab(NAbsorbGrids + NUsualGrids + 1 : NGrids, j), &
                              wfn_ab(NAbsorbGrids + NUsualGrids + 1 : NGrids, j))
            reflection(j)   = reflection(j) &
                            + dq * dot_product( &
                              wfn(1 : NAbsorbGrids, j), &
                              wfn(1 : NAbsorbGrids, j)) &
                            - dq * dot_product( &
                              wfn_ab(1 : NAbsorbGrids, j), &
                              wfn_ab(1 : NAbsorbGrids, j))
        end forall
        if (population(dq, wfn_ab, NTotal) < minpop) then
            write(*,*)"All population has been absorbed"
            write(*,*)"Stop propagation at time ", i * dt
            exit
        end if
    end do
    if (i > NSnapshots) then
        write(*,*)"Warning: not all population was absorbed within total propagation time"
        write(*,*)"Transmission + reflection may be less than 1"
    end if
    !Output
    open(unit=99, file="transmission.txt", status="replace")
        write(99,*)"state"//char(9)//"transmission"//char(9)//"reflection"
        do i = 1, NStates
            write(99,*)i, char(9), transmission(i), char(9), reflection(i)
        end do
    close(99)
    !Clean up
    deallocate(grids)
    deallocate(H)
    deallocate(wfn)
    deallocate(wfn_ab)
    contains
    real*8 function population(dq, wfn, N)
        real*8, intent(in)::dq
        integer, intent(in)::N
        complex*16, dimension(N), intent(in)::wfn
        population = dq * dot_product(wfn, wfn)
    end function population
end subroutine transmit_reflect

end module solver