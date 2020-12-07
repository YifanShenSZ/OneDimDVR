!Calculate transmission and reflection from Chebyshev order domain wave function
program main
    implicit none
    !input
    real*8::mass, dt, left, right
    !transmission and reflection
    integer::left_index, right_index
    real*8, allocatable, dimension(:)::transmission, reflection
    complex*16, allocatable, dimension(:,:)::p ! DVR momentum
    complex*16, allocatable, dimension(:)::p_left, p_right ! a row in p
    !space grid points and Chebyshev order domain wave function
    integer::NStates, order, NGrids
    real*8::Hmin, Hmax, dq
    real*8, allocatable, dimension(:)::grids
    complex*16, allocatable, dimension(:,:,:)::Chebyshev
    !transformation to time domain
    real*8::R
    complex*16::scalor
    complex*16, allocatable, dimension(:,:)::wfn
    !parallelization
    integer, external::omp_get_max_threads
    integer::NThreads
    integer, allocatable, dimension(:)::chunk_start, chunk_stop
    complex*16, allocatable, dimension(:,:,:)::wfn_parallel
    !work variable
    logical::early_stop
    integer::i, j, k
    real*8::time, coeff, population

    write(*,*)"Calculate transmission and reflection from Chebyshev order domain wave function"
    write(*,*)"Yifan Shen 2020"
    write(*,*)
    !Read input
    open(unit=99, file="Che2tran.in", status="old", iostat=i)
    if (i == 0) then
        read(99,*); read(99,*)mass
        read(99,*); read(99,*)dt
        read(99,*); read(99,*)left
        read(99,*); read(99,*)right
    else
        close(99)
        open(unit=99, file="Che2tran.in", status="replace")
            write(99,*)"Mass:"
            write(99,*)
            write(99,*)"Time step:"
            write(99,*)
            write(99,*)"Left boundary:"
            write(99,*)
            write(99,*)"Right boundary:"
            write(99,*)
        close(99)
        stop "Please fill in Che2tran.in, a template has been provided"
    end if
    close(99)

    !Read space grid points and Chebyshev order domain wave function
    open(unit=99, file="checkpoint-Chebyshev.txt")
        read(99,*); read(99,*)NStates
        read(99,*); read(99,*)order
        read(99,*); read(99,*)NGrids
        read(99,*); read(99,*)Hmin
        read(99,*); read(99,*)Hmax
    close(99)
    allocate(grids(NGrids))
    open(unit=99, file="grids.out", form="unformatted")
        read(99)grids
    close(99)
    dq = grids(2) - grids(1)
    left_index = 0
    do i = 2, NGrids
        if (grids(i) > left) then
            left_index = i - 1
            exit
        end if
    end do
    if (left_index == 0) stop "Your space grids seem not to cover the left boundary"
    right_index = 0
    do i = NGrids - 1, 1, -1
        if (grids(i) < right) then
            right_index = i + 1
            exit
        end if
    end do
    if (right_index == 0) stop "Your space grids seem not to cover the right boundary"
    allocate(Chebyshev(NGrids, NStates, 0:order))
    open(unit=99, file="Chebyshev.out", form="unformatted")
        do i = 0, order
            read(99)Chebyshev(:,:,i)
        end do
    close(99)

    !Initialize parallelization
    NThreads = omp_get_max_threads()
    allocate(chunk_start(NThreads))
    allocate(chunk_stop(NThreads))
    allocate(wfn_parallel(NGrids, NStates, NThreads))
    j = order / NThreads
    chunk_start(1) = 1
    do i = 2, NThreads
        chunk_start(i) = chunk_start(i - 1) + j
    end do
    do i = 1, NThreads - 1
        chunk_stop(i) = chunk_start(i + 1) - 1
    end do
    chunk_stop(NThreads) = order

    !Initialize transmission and reflection
    allocate(transmission(NStates))
    transmission = 0d0
    allocate(reflection(NStates))
    reflection   = 0d0
    allocate(p(NGrids, NGrids))
    call compute_momentum(dq, p, NGrids)
    allocate(p_left(NGrids))
    p_left = p(left_index, :)
    allocate(p_right(NGrids))
    p_right = p(right_index, :)

    !Transform to time domain
    allocate(wfn(NGrids, NStates))
    early_stop = .false.
    time = 0d0
    do
        time = time + dt
        scalor = exp((0d0,-1d0) * (Hmax + Hmin) * time / 2d0)
        R = (Hmax - Hmin) * time / 2d0
        wfn_parallel = (0d0,0d0)
        !$OMP PARALLEL DO PRIVATE(j, k)
        do j = 1, NThreads
            do k = chunk_start(j), chunk_stop(j)
                coeff = 2d0 * bessel_jn(k, R)
                if (abs(coeff) > 1d-8) &
                wfn_parallel(:,:,j) = wfn_parallel(:,:,j) + (0d0,1d0)**k * coeff * Chebyshev(:,:,k)
            end do
        end do
        !$OMP END PARALLEL DO
        wfn = bessel_j0(R) * Chebyshev(:,:,0)
        do j = 1, NThreads
            wfn = wfn + wfn_parallel(:,:,j)
        end do
        !Calculate absorption
        do j = 1, NStates
            transmission(j) = transmission(j) &
                            - dble(wfn(right_index, j) * dot_product(p_right, conjg(wfn(:,j))) &
                            - conjg(wfn(right_index, j)) * dot_product(p_right, wfn(:,j)))
            reflection(j) = reflection(j) &
                          + dble(wfn(left_index, j) * dot_product(p_left, conjg(wfn(:,j))) &
                          - conjg(wfn(left_index, j)) * dot_product(p_left, wfn(:,j)))
        end do
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
            write(*,*)"Stop propagation at time ", time
            exit
        end if
    end do
    transmission = transmission * dt / 2d0 / mass
    reflection   = reflection   * dt / 2d0 / mass
    population = sum(transmission) + sum(reflection)
    if (population < 0.9999) then
        write(*,*)"Warning: sum of transmission and reflection = ", population
        if (early_stop) then
            write(*,*)"Maybe grid spacing and/or output interval is not sufficiently small"
        else
            write(*,*)"Probably because not all population was absorbed"
        end if
    else if (population > 1.0001) then
        write(*,*)"Warning: sum of transmission and reflection = ", population
        write(*,*)"Maybe grid spacing and/or output interval is not sufficiently small"
    end if

    !Output
    open(unit=99, file="transmission.txt", status="replace")
        write(99,*)"state"//char(9)//"transmission"//char(9)//"reflection"
        do i = 1, NStates
            write(99,*)i, char(9), transmission(i), char(9), reflection(i)
        end do
    close(99)

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

subroutine compute_population(population)
    real*8, intent(out)::population
    integer::i
    population = 0d0
    do i = 1, NStates
        population = population + dot_product(wfn(left_index+1:right_index-1, i), wfn(left_index+1:right_index-1, i))
    end do
    population = population * dq
end subroutine compute_population

end program main