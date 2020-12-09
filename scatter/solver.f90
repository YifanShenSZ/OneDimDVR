module solver
    use LinearAlgebra
    use basic
    use DVR
    implicit none

!Parameter
    logical::auto_spacing = .true.!Determine grid spacing by absorbance condition

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
    !DVR Hamiltonian: kinetic energy T, potential energy V
    complex*16, allocatable, dimension(:,:)::T
    complex*16, allocatable, dimension(:,:,:)::V
    !DVR wave function
    real*8::population
    complex*16, allocatable, dimension(:,:)::wfn
    !Work variable
    logical::early_stop
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
    !This -1 prevents the grids from touching the boundary
    !since absorbing potential diverges there
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
    !Build DVR Hamiltonian
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
        write(99)wfn
        !Propagate by RK4
        early_stop = .false.
        do i = 1, (NSnapshots - 1) * OutputStep
            call RK4(T, V, wfn, wfn, dt, NGrids, NStates)
            !Write trajectory
            if(mod(i, OutputStep) == 0) then
                write(99)wfn
            end if
            !Check population
            call compute_population(population) 
            if (population > 1.0001) then
                write(*,*)"Error: total population = ", population
                write(*,*)"Maybe grid spacing and/or time step is not sufficiently small"
                early_stop = .true.
            else if (population < 0.0001) then
                write(*,*)"All population has been absorbed"
                early_stop = .true.
            end if
            if (early_stop) then
                write(*,*)"Stop propagation at time ", i * dt
                NSnapshots = 1 + i / OutputStep
                exit
            end if
        end do
        if (.not. early_stop) then
            write(*,*)"Warning: not all population was absorbed within total propagation time"
            write(*,*)population, " population remains"
            write(*,*)"The scattering process has not reached its end yet"
        end if
    close(99)
    !Leave a checkpoint
    open(unit=99, file="checkpoint.txt", status="replace")
        write(99,*)"Number of electronic states:"
        write(99,*)NStates
        write(99,*)"Number of time snapshots:"
        write(99,*)NSnapshots
        write(99,*)"Number of grid points:"
        write(99,*)NGrids
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
        write(99)grids
    close(99)
    !Clean up
    deallocate(snapshots)
    deallocate(grids)
    deallocate(T)
    deallocate(V)
    deallocate(wfn)
    contains
    subroutine compute_population(population)
        real*8, intent(out)::population
        integer::i
        population = 0d0
        do i = 1, NStates
            population = population + dot_product( &
                         wfn(NAbsorbGrids + 1 : NAbsorbGrids + NUsualGrids, i), &
                         wfn(NAbsorbGrids + 1 : NAbsorbGrids + NUsualGrids, i))
        end do
        population = population * dq
    end subroutine compute_population
end subroutine propagate_wavefunction

!The difference to subroutine propagate_wavefunction is
!this subroutine cares only about transmission and reflection rather than saving the wavefunction
subroutine transmit_reflect()
    !grid points
    integer::NSnapshots, NUsualGrids, NAbsorbGrids, NGrids
    real*8, allocatable, dimension(:)::grids
    !DVR Hamiltonian: kinetic energy T, potential energy V
    complex*16, allocatable, dimension(:,:)::T
    complex*16, allocatable, dimension(:,:,:)::V
    !DVR wave function
    real*8::population
    complex*16, allocatable, dimension(:, :)::wfn
    !transmission and reflection
    integer::left_index, right_index
    complex*16, allocatable, dimension(:,:)::p ! DVR momentum
    complex*16, allocatable, dimension(:)::p_left, p_right ! rows in p
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
    !This -1 prevents the grids from touching the boundary
    !since absorbing potential diverges there
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
    !Build DVR Hamiltonian
    allocate(T(NGrids, NGrids))
    call compute_kinetic(dq, T, NGrids)
    allocate(V(NGrids, NStates, NStates))
    call compute_potential(grids, V, NGrids, NStates)
    !Initialize transmission and reflection
    left_index = NAbsorbGrids + 1
    right_index = NAbsorbGrids + NUsualGrids
    allocate(p(NGrids, NGrids))
    call compute_momentum(dq, p, NGrids)
    allocate(p_left(NGrids))
    p_left = p(left_index, :)
    allocate(p_right(NGrids))
    p_right = p(right_index, :)
    transmission = 0d0
    reflection   = 0d0
    !Set initial condition
    allocate(wfn(NGrids, NStates))
    do i = 1, NGrids
        call init_wfn(grids(i), wfn(i, :))
    end do
    !Propagate wave function
    do i = 1, NSnapshots
        call RK4(T, V, wfn, wfn, dt, NGrids, NStates)
        !Calculate transmission and reflection
        forall (j = 1 : NStates)
            transmission(j) = transmission(j) - dble(wfn(right_index, j) * dot_product(wfn(:,j), p_right))
            reflection  (j) = reflection  (j) + dble(wfn(left_index , j) * dot_product(wfn(:,j), p_left ))
        end forall
        !Check population
        call compute_population(population) 
        if (population > 1.0001) then
            write(*,*)"Error: total population = ", population
            write(*,*)"Maybe grid spacing and/or time step is not sufficiently small"
            write(*,*)"Stop propagation at time ", i * dt
            exit
        else if (population < 0.0001) then
            write(*,*)"All population has been absorbed"
            write(*,*)"Stop propagation at time ", i * dt
            exit
        end if
    end do
    if (i > NSnapshots) then
        write(*,*)"Warning: not all population was absorbed within total propagation time"
        write(*,*)population, " population remains"
        write(*,*)"The scattering process has not reached its end yet"
    end if
    transmission = transmission * dt / mass
    reflection   = reflection   * dt / mass
    population = sum(transmission) + sum(reflection)
    if (population < 0.9999) then
        write(*,*)"Warning: sum of transmission and reflection = ", population
        if (i <= NSnapshots) then
            write(*,*)"Maybe grid spacing and/or time step is not sufficiently small"
        else
            write(*,*)"Probably because not all population was absorbed"
        end if
    else if (population > 1.0001) then
        write(*,*)"Warning: sum of transmission and reflection = ", population
        write(*,*)"Maybe grid spacing and/or time step is not sufficiently small"
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
    deallocate(p)
    deallocate(p_left)
    deallocate(p_right)
    contains
    subroutine compute_population(population)
        real*8, intent(out)::population
        integer::i
        population = 0d0
        do i = 1, NStates
            population = population + dot_product( &
                         wfn(NAbsorbGrids + 1 : NAbsorbGrids + NUsualGrids, i), &
                         wfn(NAbsorbGrids + 1 : NAbsorbGrids + NUsualGrids, i))
        end do
        population = population * dq
    end subroutine compute_population
end subroutine transmit_reflect

end module solver