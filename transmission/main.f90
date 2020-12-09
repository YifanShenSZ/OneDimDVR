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
    complex*16, allocatable, dimension(:,:)::p ! DVR momentum
    complex*16, allocatable, dimension(:)::p_left, p_right ! rows in p
    real*8, allocatable, dimension(:)::transmission, reflection
    !work variable
    integer::i, j
    real*8::population

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
    allocate(p(NGrids, NGrids))
    call compute_momentum(dq, p, NGrids)
    allocate(p_left(NGrids))
    p_left = p(left_index, :)
    allocate(p_right(NGrids))
    p_right = p(right_index, :)
    allocate(transmission(NStates))
    allocate(reflection  (NStates))
    transmission = 0d0
    reflection   = 0d0
    open(unit=99,  file="wfn.out", form="unformatted")
        do i = 1, NSnapshots
            read(99)wfn
            forall (j = 1 : NStates)
                transmission(j) = transmission(j) - dble(wfn(right_index, j) * dot_product(wfn(:,j), p_right))
                reflection  (j) = reflection  (j) + dble(wfn(left_index , j) * dot_product(wfn(:,j), p_left ))
            end forall
        end do
    close(99)
    transmission = transmission * dt / mass
    reflection   = reflection   * dt / mass
    population = sum(transmission) + sum(reflection)
    if (population < 0.9999) then
        write(*,*)"Warning: sum of transmission and reflection = ", population
        write(*,*)"Maybe not all population has left interaction region"
        write(*,*)"Maybe grid spacing and/or output interval is not sufficiently small"
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

end program main