!Convert OneDimDVR wave function to Wigner distribution
program main
    implicit none
    !input
    real*8::pmin, pmax, dp
    !space grid points and wave function from wave function calculation
    integer::NStates, NSnapshots, NGrids
    real*8, allocatable, dimension(:)::grids
    complex*16, allocatable, dimension(:,:)::wfn
    !momentum grid points
    integer::NMomenta
    real*8, allocatable, dimension(:)::momenta
    !Wigner distribution
    real*8, allocatable, dimension(:,:)::Wigner
    !work variable
    character*128::count
    integer::i, j

    write(*,*)"Convert OneDimDVR wave function to Wigner distribution"
    write(*,*)"Yifan Shen 2020"
    write(*,*)
    !Read input
    open(unit=99, file="Wigner.in", status="old", iostat=i)
    if (i == 0) then
        read(99,*); read(99,*)pmin
        read(99,*); read(99,*)pmax
        read(99,*); read(99,*)dp
    else
        close(99)
        open(unit=99, file="Wigner.in", status="replace")
            write(99,*)"Minimum momentum:"
            write(99,*)
            write(99,*)"Maximum momentum:"
            write(99,*)
            write(99,*)"Momentum step:"
            write(99,*)
        close(99)
        stop "Please fill in Wigner.in, a template has been provided"
    end if
    close(99)
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

    !Discretize momentum space
    NMomenta = floor((pmax - pmin) / dp) + 1
    allocate(momenta(NMomenta))
    momenta(1) = pmin
    do i = 2, NMomenta
        momenta(i) = momenta(i - 1) + dp
    end do

    !Calculate Wigner distribution
    write(*,*)"Calculating Wigner distribution..."
    allocate(Wigner(NGrids, NMomenta))
    open(unit=99, file="wfn.out", form="unformatted")
    do j = 1, NStates
        write(count,*)j
        open(unit=99+j, file="Wigner"//trim(adjustl(count))//".out", form="unformatted", status="replace")
    end do
    do i = 1, NSnapshots
        read(99)wfn
        do j = 1, NStates
            call wfn2Wigner(grids, momenta, wfn(:, j), Wigner, NGrids, NMomenta)
            write(99+j)Wigner
        end do
    end do
    close(99)
    do j = 1, NStates
        close(99+j)
    end do

    !Output grid information
    open(unit=99, file="momenta.out", form="unformatted", status="replace")
        write(99)momenta
    close(99)
    open(unit=99, file="checkpoint.txt", access="append")
        write(99,*)"Number of momentum grid points:"
        write(99,*)NMomenta
    close(99)

    write(*,*)
    write(*,*)"Mission success"

contains
!Wigner(q, p) = WignerTransform[wfn(q)]
subroutine wfn2Wigner(grids, momenta, wfn, Wigner, NGrids, NMomenta)
    integer, intent(in)::NGrids, NMomenta
    real*8, dimension(NGrids), intent(in)::grids
    real*8, dimension(NMomenta), intent(in)::momenta
    complex*16, dimension(NGrids), intent(in)::wfn
    real*8, dimension(NGrids, NMomenta), intent(out)::Wigner
    complex*16, parameter::cim2 = (0d0,2d0)
    integer::i, j, k
    real*8::dq
    dq = grids(2) - grids(1)
    !$OMP PARALLEL DO PRIVATE(i, j, k)
    do i = 1, NMomenta
        do j = 1, NGrids
            Wigner(j, i) = dble(conjg(wfn(j)) * wfn(j))
            do k = 1, min(NGrids - j, j - 1)
                Wigner(j, i) = Wigner(j, i) + 2d0 * dble( &
                          + exp(cim2 * momenta(i) * dble(k) * dq) &
                          * conjg(wfn(j + k)) * wfn(j - k) )
            end do
        end do
    end do
    !$OMP END PARALLEL DO
    Wigner = Wigner * dq / 3.141592653589793d0
end subroutine wfn2Wigner

end program main