module solver
    use LinearAlgebra
    use basic
    use DVR
    implicit none

contains
!Perform matrix vector multiplication 2H . wfn
subroutine perform_2Hwfn(output, T, V, wfn, NGrids, NStates)
    integer,intent(in)::NGrids, NStates
    complex*16, dimension(NGrids, NStates), intent(out)::output
    complex*16, dimension(NGrids, NGrids), intent(in)::T
    complex*16, dimension(NGrids, NStates, NStates), intent(in)::V
    complex*16, dimension(NGrids, NStates), intent(in)::wfn
    integer::i, j
    !$OMP PARALLEL DO PRIVATE(i, j)
    do i = 1, NStates
        !V . wfn
        output(:,i) = V(:,i,1) * wfn(:,1)
        do j = 2, NStates
            if (j > i) then
                output(:,i) = output(:,i) + V(:,j,i) * wfn(:,j)
            else
                output(:,i) = output(:,i) + V(:,i,j) * wfn(:,j)
            end if
        end do
        !2T . wfn + 2V . wfn
        call zsymv('L', NGrids, (2d0,0d0), T, NGrids, wfn(:,i), 1, (2d0,0d0), output(:,i), 1)
    end do
    !$OMP END PARALLEL DO
end subroutine perform_2Hwfn

!Recursively determine the new Chebyshev wave function
!V. A. Mandelshtam 1995 J. Chem. Phys.
!H. Guo 2006 Phys. Rev. A
subroutine recurse(new, D, T, V, current, old, NGrids, NStates)
    integer,intent(in)::NGrids, NStates
    complex*16, dimension(NGrids, NStates), intent(out)::new
    complex*16, dimension(NGrids), intent(in)::D
    complex*16, dimension(NGrids, NGrids), intent(in)::T
    complex*16, dimension(NGrids, NStates, NStates), intent(in)::V
    complex*16, dimension(NGrids, NStates), intent(in)::current, old
    integer::i
    call perform_2Hwfn(new, T, V, current, NGrids, NStates)
    forall (i = 1 : NStates)
        new(:,i) = D * (new(:,i) - D * old(:,i))
    end forall
end subroutine recurse

!This is a Chebyshev solver:
!    1. Build Hamiltonian matrix
!    2. Propagate by Chebyshev recursion
!Note that wave function is propagated in Chebyshev order domain
subroutine propagate_Chebyshev()
    !Grid points
    integer::NUsualGrids, NAbsorbGrids, NGrids
    real*8, allocatable, dimension(:)::grids
    !DVR Hamiltonian: kinetic energy T, potential energy V, damping D
    complex*16, allocatable, dimension(:,:)::T
    complex*16, allocatable, dimension(:,:,:)::V
    complex*16, allocatable, dimension(:)::D
    !DVR wave function
    real*8::population
    complex*16, allocatable, dimension(:,:)::wfnnew, wfn, wfnold
    !Chebyshev scalor
    real*8::Hmax, Hmin, Tmax, Tmin, Vmax, Vmin, shift
    real*8, dimension(NStates)::energy
    real*8, dimension(NStates, NStates)::Hd
    !Work variable
    integer::i, j
    !Discretize space
    NUsualGrids = floor((right - left) / dq) + 1
    dq = (right - left) / (NUsualGrids - 1)
    NAbsorbGrids = floor(6.283185307179586d0 / kmin / dq)
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
    allocate(D(NGrids))
    call compute_damping(grids, D, NGrids)
    !Scale Hamiltonian to [-1, 1]
    Tmax = 39.47841760435743d0 / dq / dq / 2d0 / mass ! (2pi/dq)^2 / 2m
    Tmin = Tmax / NGrids / NGrids
    Hd = V(1,:,:)
    call My_dsyev('N', Hd, energy, NStates)
    Vmin = energy(1)
    Vmax = energy(NStates)
    do i = 2, NGrids
        Hd = V(i,:,:)
        call My_dsyev('N', Hd, energy, NStates)
        if (energy(1) < Vmin) Vmin = energy(1)
        if (energy(NStates) > Vmax) Vmax = energy(NStates)
    end do
    Hmax = Tmax + Vmax
    Hmin = Tmin + Vmin
    T = -2d0 / (Hmax - Hmin) * T
    V = -2d0 / (Hmax - Hmin) * V
    shift = (Hmax + Hmin) / (Hmax - Hmin)
    forall (i = 1 : NStates)
        V(:,i,i) = V(:,i,i) + shift
    end forall
    !Propagate wave function in Chebyshev order domain
    open(unit=99, file="Chebyshev.out", form="unformatted", status="replace")
        !Initial condition
        allocate(wfnnew(NGrids, NStates))
        allocate(wfn   (NGrids, NStates))
        allocate(wfnold(NGrids, NStates))
        do i = 1, NGrids
            call init_wfn(grids(i), wfnold(i, :))
        end do
        write(99)wfnold
        !The 1st recursion
        call perform_2Hwfn(wfn, T, V, wfnold, NGrids, NStates)
        forall (i = 1 : NStates)
            wfn(:,i) = D * wfn(:,i) / 2d0
        end forall
        write(99)wfn
        !Propagate by Chebyshev recursion
        do i = 2, order
            call recurse(wfnnew, D, T, V, wfn, wfnold, NGrids, NStates)
            !Write trajectory
            write(99)wfnnew
            !Check population
            call compute_population(population) 
            if (population < 0.0001) then
                write(*,*)"All population has been absorbed"
                write(*,*)"Stop propagation at order ", i
                order = i
                exit
            end if
            !Get ready for the next loop
            wfnold = wfn
            wfn = wfnnew
        end do
        if (i > order) then
            write(*,*)"Warning: not all population was absorbed within total order"
            write(*,*)population, " population remains"
            write(*,*)"The scattering process has not reached its end yet"
        end if
    close(99)
    !Leave a checkpoint
    open(unit=99, file="checkpoint-Chebyshev.txt", status="replace")
        write(99,*)"Number of electronic states:"
        write(99,*)NStates
        write(99,*)"Order of Chebyshev propagation:"
        write(99,*)order
        write(99,*)"Number of grid points:"
        write(99,*)NGrids
        write(99,*)"Minimum possible energy:"
        write(99,*)Hmin
        write(99,*)"Maximum possible energy:"
        write(99,*)Hmax
    close(99)
    !Output the grids used in calculation
    open(unit=99, file="grids.out", form="unformatted", status="replace")
        write(99)grids
    close(99)
    !Clean up
    deallocate(grids)
    deallocate(T)
    deallocate(V)
    deallocate(D)
    deallocate(wfnnew)
    deallocate(wfn)
    deallocate(wfnold)
    contains
    subroutine compute_population(population)
        real*8, intent(out)::population
        integer::i
        population = 0d0
        do i = 1, NStates
            population = population + dot_product( &
                         wfnnew(NAbsorbGrids + 1 : NAbsorbGrids + NUsualGrids, i), &
                         wfnnew(NAbsorbGrids + 1 : NAbsorbGrids + NUsualGrids, i))
        end do
        population = population * dq
    end subroutine compute_population
end subroutine propagate_Chebyshev

end module solver