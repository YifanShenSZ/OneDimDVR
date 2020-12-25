!Convert OneDimDVR wave function to density
program main
    implicit none
    !space grid points and wave function from wave function calculation
    integer::NStates, NSnapshots, NGrids
    real*8, allocatable, dimension(:)::grids
    complex*16, allocatable, dimension(:, :)::wfn
    !work variable
    character*128::count
    integer::i, j

    write(*,*)"Convert OneDimDVR wave function to density"
    write(*,*)"Yifan Shen 2020"
    write(*,*)
    !Read space grid points from wave function calculation
    open(unit=99, file="checkpoint.txt")
        read(99,*); read(99,*)NStates    
        read(99,*); read(99,*)NSnapshots
        read(99,*); read(99,*)NGrids
    close(99)
    allocate(grids(NGrids))
    open(unit=99, file="grids.out", form="unformatted")
        read(99)grids
    close(99)
    allocate(wfn(NGrids, NStates))

    !Calculate density
    open(unit=99, file="wfn.out", form="unformatted")
    do j = 1, NStates
        write(count,*)j
        open(unit=99+j, file="density"//trim(adjustl(count))//".out", form="unformatted", status="replace")
    end do
    do i = 1, NSnapshots
        read(99)wfn
        do j = 1, NStates
            write(99+j)dble(conjg(wfn(:, j)) * wfn(:, j))
        end do
    end do
    close(99)
    do j = 1, NStates
        close(99+j)
    end do

end program main