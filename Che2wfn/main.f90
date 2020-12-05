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
    !work variable
    integer::i, j

    write(*,*)"Convert OneDimDVR Chebyshev order domain wave function to time domain"
    write(*,*)"Yifan Shen 2020"
    write(*,*)
    !Read input
    open(unit=99, file="Che2wfn.in")
        read(99,*); read(99,*)total_time
        read(99,*); read(99,*)output_interval
    close(99)
    !Discretize time
    NSnapshots = floor(total_time / output_interval) + 1
    allocate(snapshots(NSnapshots))
    snapshots(1) = 0d0
    do i = 2, NSnapshots
        snapshots(i) = snapshots(i - 1) + output_interval
    end do

    !Read Chebyshev order domain wave function
    open(unit=99, file="Chebyshev-checkpoint.txt")
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

    !Transform to time domain
    allocate(wfn(NGrids, NStates))
    open(unit=99, file="wfn.out", form="unformatted", status="replace")
    do i = 1, NSnapshots
        scalor = exp((0d0,-1d0) * (Hmax + Hmin) * snapshots(i) / 2d0)
        R = (Hmax - Hmin) * snapshots(i) / 2d0
        wfn = bessel_j0(R) * Chebyshev(:,:,0)
        do j = 1, order
            wfn = wfn + 2d0 * (0d0,1d0)**j * bessel_jn(j, R) * Chebyshev(:,:,j)
        end do
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