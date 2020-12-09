!Convert OneDimDVR Chebyshev order domain wave function to time domain
program main
    implicit none
    !time grids
    real*8::total_time, output_interval
    integer::NSnapshots
    real*8, allocatable, dimension(:)::snapshots
    !Chebyshev order domain wave function
    integer::NStates, order, NGrids
    real*8::Hmin, Hmax
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
    integer::i, j, k
    real*8::coeff

    write(*,*)"Convert OneDimDVR Chebyshev order domain wave function to time domain"
    write(*,*)"Yifan Shen 2020"
    write(*,*)
    !Read input
    open(unit=99, file="Che2wfn.in", status="old", iostat=i)
    if (i == 0) then
        read(99,*); read(99,*)total_time
        read(99,*); read(99,*)output_interval
    else
        close(99)
        open(unit=99, file="Che2wfn.in", status="replace")
            write(99,*)"Total propagation time:"
            write(99,*)
            write(99,*)"Output interval:"
            write(99,*)
        close(99)
        stop "Please fill in Che2wfn.in, a template has been provided"
    end if
    close(99)
    !Discretize time
    NSnapshots = floor(total_time / output_interval) + 1
    allocate(snapshots(NSnapshots))
    snapshots(1) = 0d0
    do i = 2, NSnapshots
        snapshots(i) = snapshots(i - 1) + output_interval
    end do

    !Read Chebyshev order domain wave function
    open(unit=99, file="checkpoint-Chebyshev.txt")
        read(99,*); read(99,*)NStates
        read(99,*); read(99,*)order
        read(99,*); read(99,*)NGrids
        read(99,*); read(99,*)Hmin
        read(99,*); read(99,*)Hmax
    close(99)
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

    !Transform to time domain
    allocate(wfn(NGrids, NStates))
    open(unit=99, file="wfn.out", form="unformatted", status="replace")
    do i = 1, NSnapshots
        scalor = exp((0d0,-1d0) * (Hmax + Hmin) * snapshots(i) / 2d0)
        R = (Hmax - Hmin) * snapshots(i) / 2d0
        wfn_parallel = (0d0,0d0)
        !$OMP PARALLEL DO PRIVATE(j, k, coeff)
        do j = 1, NThreads
            do k = chunk_start(j), chunk_stop(j)
                coeff = 2d0 * bessel_jn(k, R)
                !Only |coeff| > 1d-8 terms are taken into account
                if (abs(coeff) > 1d-8) &
                wfn_parallel(:,:,j) = wfn_parallel(:,:,j) + (0d0,1d0)**k * coeff * Chebyshev(:,:,k)
            end do
        end do
        !$OMP END PARALLEL DO
        wfn = bessel_j0(R) * Chebyshev(:,:,0) + sum(wfn_parallel, 3)
        write(99)wfn
    end do
    close(99)

    !Output
    open(unit=99, file="checkpoint.txt", status="replace")
        write(99,*)"Number of electronic states:"
        write(99,*)NStates
        write(99,*)"Number of time snapshots:"
        write(99,*)NSnapshots
        write(99,*)"Number of grid points:"
        write(99,*)NGrids
    close(99)
    open(unit=99, file="snapshots.out", form="unformatted", status="replace")
        write(99)snapshots
    close(99)

end program main