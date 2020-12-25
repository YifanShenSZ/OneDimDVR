!Calculate time dependence of population of each state from OneDimDVR wave function
program main
    implicit none
    !grid points and wave function from wave function calculation
    integer::NStates, NSnapshots, NGrids
    real*8::dq
    real*8, allocatable, dimension(:)::snapshots, grids
    complex*16, allocatable, dimension(:, :)::wfn
    !work variable
    integer::i, j

    write(*,'(A)')"Calculate time dependence of population of each state from OneDimDVR wave function"
    write(*,'(A)')"Yifan Shen 2020"
    write(*,*)
    !Read space grid points from wave function calculation
    open(unit=99, file="checkpoint.txt")
        read(99,*); read(99,*)NStates    
        read(99,*); read(99,*)NSnapshots
        read(99,*); read(99,*)NGrids
    close(99)
    allocate(snapshots(NSnapshots))
    open(unit=99, file="snapshots.out", form="unformatted", status="old")
        read(99)snapshots
    close(99)
    allocate(grids(NGrids))
    open(unit=99, file="grids.out", form="unformatted")
        read(99)grids
    close(99)
    dq = grids(2) - grids(1)
    allocate(wfn(NGrids, NStates))

    !Calculate density
    open(unit=99,  file="wfn.out", form="unformatted")
    open(unit=100, file="population.txt", status="replace")
    do i = 1, NSnapshots
        read(99)wfn
        write(100,"(f14.8)",advance="no")snapshots(i)
        do j = 1, NStates
            write(100,"(f14.8)",advance="no")dq * dble(dot_product(wfn(:, j), wfn(:, j)))
        end do
        write(100,*)
    end do
    close(100); close(99)

end program main