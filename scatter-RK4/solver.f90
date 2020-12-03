module solver
    use LinearAlgebra
    use basic
    use DVR
    implicit none

!Parameter
    logical::auto_spacing = .true.!Determine grid spacing by absorbance condition
    real*8::maxpopdev = 1d-6, &!Stop wave function propagation when population no longer equals to 1
            minpop    = 1d-4   !Stop transmission & reflection calculation when all population are absorbed
                               !Absorbing potential can only absorb 99% of kmin (larger k has better absorption)

contains
subroutine RK4(T, V, old, new, dt, NGrids, NStates)
    integer,intent(in)::NGrids, NStates
    complex*16, dimension(NGrids, NGrids), intent(in)::T
    complex*16, dimension(NGrids, NStates, NStates), intent(in)::V
    complex*16, dimension(NGrids, NStates), intent(in)::old
    complex*16, dimension(NGrids, NStates), intent(out)::new
    real*8, intent(in)::dt
    real*8::dtd2
    complex*16, dimension(NGrids, NStates)::k1, k2, k3, k4
    dtd2 = dt / 2d0
    call TDSE(k1, old            )
    call TDSE(k2, old + k1 * dtd2)
    call TDSE(k3, old + k2 * dtd2)
    call TDSE(k4, old + k3 * dt  )
    new = old + dt / 6d0 * (k1 + 2d0 * k2 + 2d0 * k3 + k4)
    contains
    !Time-dependent Schrodinger equation
    subroutine TDSE(gradient, wfn)
        complex*16, dimension(NGrids, NStates), intent(out)::gradient
        complex*16, dimension(NGrids, NStates), intent(in)::wfn
        integer::i, j
        !$OMP PARALLEL DO PRIVATE(i, j)
        do i = 1, NStates
            !V . wfn
            gradient(:,i) = V(:,i,1) * wfn(:,1)
            do j = 2, NStates
                if (j > i) then
                    gradient(:,i) = gradient(:,i) + V(:,j,i) * wfn(:,j)
                else
                    gradient(:,i) = gradient(:,i) + V(:,i,j) * wfn(:,j)
                end if
            end do
            !-i * T . wfn + -i * V . wfn
            call zsymv('L', NGrids, (0d0,-1d0), T, NGrids, wfn(:,i), 1, (0d0,-1d0), gradient(:,i), 1)
        end do
        !$OMP END PARALLEL DO
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
    !DVR Hamiltonian: kinetic energy T and potential energy V
    complex*16, allocatable, dimension(:,:)::T
    complex*16, allocatable, dimension(:,:,:)::V
    !DVR wave function
    complex*16, allocatable, dimension(:,:)::wfn
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
    allocate(T(NGrids, NGrids))
    call compute_kinetic(dq, T, NGrids)
    allocate(V(NGrids, NStates, NStates))
    call compute_potential(grids, V, NGrids, NStates)
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
            call RK4(T, V, wfn, wfn, dt, NGrids, NStates)
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
    deallocate(T)
    deallocate(V)
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
    !DVR matrices: kinetic energy T, momentum p, potential energy V
    complex*16, allocatable, dimension(:,:)::T, p
    complex*16, allocatable, dimension(:,:,:)::V
    !DVR wave function
    complex*16, allocatable, dimension(:, :)::wfn
    !transmission and reflection
    real*8, dimension(NStates)::transmission, reflection
    complex*16, allocatable, dimension(:)::pwfn, pwfnstar
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
    !Build DVR matrices
    allocate(T(NGrids, NGrids))
    call compute_kinetic(dq, T, NGrids)
    allocate(p(NGrids, NGrids))
    call compute_momentum(dq, p, NGrids)
    allocate(V(NGrids, NStates, NStates))
    call compute_potential(grids, V, NGrids, NStates)
    !Set initial condition
    allocate(wfn(NGrids, NStates))
    do i = 1, NGrids
        call init_wfn(grids(i), wfn(i, :))
    end do
    transmission = 0d0
    reflection   = 0d0
    allocate(pwfn    (NGrids))
    allocate(pwfnstar(NGrids))
    !Propagate wave function
    do i = 1, NSnapshots
        call RK4(T, V, wfn, wfn, dt, NGrids, NStates)
        !Calculate absorption
        do j = 1, NStates
            call zhemv('L', NGrids, (1d0,0d0), p, NGrids,       wfn(:,j) , 1, (0d0,0d0), pwfn    , 1)
            call zhemv('L', NGrids, (1d0,0d0), p, NGrids, conjg(wfn(:,j)), 1, (0d0,0d0), pwfnstar, 1)
            transmission(j) = transmission(j) - dt / 2d0 / mass &
                            * (wfn(NAbsorbGrids + NUsualGrids, j) * pwfnstar(NAbsorbGrids + NUsualGrids) &
                            - conjg(wfn(NAbsorbGrids + NUsualGrids, j)) * pwfn(NAbsorbGrids + NUsualGrids))
            reflection(j) = reflection(j) + dt / 2d0 / mass &
                            * (wfn(NAbsorbGrids + 1, j) * pwfnstar(NAbsorbGrids + 1) &
                            - conjg(wfn(NAbsorbGrids + 1, j)) * pwfn(NAbsorbGrids + 1))
        end do
        if (population() < minpop) then
            write(*,*)"All population has been absorbed"
            write(*,*)"Stop propagation at time ", i * dt
            exit
        end if
    end do
    if (i > NSnapshots) then
        write(*,*)"Warning: not all population was absorbed within total propagation time"
        write(*,*)population(), " population remains"
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
    deallocate(T)
    deallocate(V)
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
end subroutine transmit_reflect

end module solver