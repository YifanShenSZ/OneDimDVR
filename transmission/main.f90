!Calculate transmission and reflection from OneDimDVR wave function
program main
    implicit none
    !mass and interaction region boundary
    real*8::mass, left, right
    !space grid points and wave function from wave function calculation
    integer::NStates, NSnapshots, NGrids
    real*8::dt, dq
    real*8, allocatable, dimension(:)::snapshots, grids
    complex*16, allocatable, dimension(:, :)::wfn
    !transmission and reflection
    integer::left_index, right_index
    complex*16, allocatable, dimension(:,:)::p
    complex*16, allocatable, dimension(:)::pwfn, pwfnstar
    real*8, allocatable, dimension(:)::transmission, reflection
    !work variable
    integer::i, j
    real*8::check

    write(*,*)"Calculate transmission and reflection from OneDimDVR wave function"
    write(*,*)"Yifan Shen 2020"
    write(*,*)
    !Read input
    open(unit=99, file="transmission.in", status="old", iostat=i)
    if (i == 0) then
        read(99,*); read(99,*)mass
        read(99,*); read(99,*)left
        read(99,*); read(99,*)right
    else
        close(99)
        open(unit=99, file="transmission.in", status="replace")
            write(99,*)"Mass:"
            write(99,*)
            write(99,*)"Left boundary:"
            write(99,*)
            write(99,*)"Right boundary:"
            write(99,*)
        close(99)
        stop "Please fill in transmission.in, a template has been provided"
    end if
    close(99)
    !Read space grid points from wave function calculation
    open(unit=99, file="checkpoint.txt")
        read(99,*); read(99,*)NStates    
        read(99,*); read(99,*)NSnapshots
        read(99,*); read(99,*)NGrids
    close(99)
    allocate(snapshots(NSnapshots))
    open(unit=99, file="snapshots.out", form="unformatted")
        read(99)snapshots
    close(99)
    dt = snapshots(2) - snapshots(1)
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
    allocate(wfn(NGrids, NStates))

    !Calculate transmission and reflection
    write(*,*)"Calculating transmission and reflection..."
    allocate(p(NGrids, NGrids))
    call compute_momentum(dq, p, NGrids)
    allocate(pwfn    (NGrids))
    allocate(pwfnstar(NGrids))
    allocate(transmission(NStates))
    allocate(reflection  (NStates))
    transmission = 0d0
    reflection   = 0d0
    open(unit=99,  file="wfn.out", form="unformatted")
        do i = 1, NSnapshots
            read(99)wfn
            do j = 1, NStates
                call zhemv('L', NGrids, (1d0,0d0), p, NGrids,       wfn(:,j) , 1, (0d0,0d0), pwfn    , 1)
                call zhemv('L', NGrids, (1d0,0d0), p, NGrids, conjg(wfn(:,j)), 1, (0d0,0d0), pwfnstar, 1)
                transmission(j) = transmission(j) - dt / 2d0 / mass &
                                * dble(wfn(right_index, j) * pwfnstar(right_index) &
                                - conjg(wfn(right_index, j)) * pwfn(right_index))
                reflection(j) = reflection(j) + dt / 2d0 / mass &
                                * dble(wfn(left_index, j) * pwfnstar(left_index) &
                                - conjg(wfn(left_index, j)) * pwfn(left_index))
            end do
        end do
    close(99)

    !Output
    check = sum(transmission) + sum(reflection)
    if (check < 0.9999) then
        write(*,*)"Warning: sum of transmission and reflection = ", check
        write(*,*)"Maybe not all population has left interaction region"
        write(*,*)"Please try increasing total propagation time"
    else if (check > 1.0001) then
        write(*,*)"Warning: sum of transmission and reflection = ", check
        write(*,*)"Maybe time step is not sufficiently small"
        write(*,*)"Please try decreasing output interval"
    end if
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
end subroutine compute_momentum

end program main