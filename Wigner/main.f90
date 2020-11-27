!Convert OneDimDVR wave function to Wigner distribution
program main
    implicit none
    !input
    real*8::pmin, pmax, dp
    integer::grids_increment
    !space grid points and wave function from wave function calculation
    integer::NStates, NSnapshots, NGrids
    real*8, allocatable, dimension(:)::grids
    complex*16, allocatable, dimension(:)::wfn
    !momentum grid points
    integer::NMomenta
    real*8, allocatable, dimension(:)::momenta
    !Wigner space grid points
    integer::NWignerGrids
    !Wigner distribution
    real*8, allocatable, dimension(:, :)::Wigner
    !work variable
    character*128::count
    integer::i, j

    write(*,*)"Convert OneDimDVR wave function to Wigner distribution"
    write(*,*)"Yifan Shen 2020"
    write(*,*)
    !Read input
    open(unit=99, file="OneDimDVR-Wigner.in")
        read(99,*); read(99,*)pmin
        read(99,*); read(99,*)pmax
        read(99,*); read(99,*)dp
        read(99,*); read(99,*)grids_increment
    close(99)
    !Read space grid points from wave function calculation
    open(unit=99, file="checkpoint.out")
        read(99,*); read(99,*)NStates    
        read(99,*); read(99,*)NSnapshots
        read(99,*); read(99,*)NGrids
    close(99)
    allocate(grids(NGrids))
    open(unit=99, file="grids.out", form="unformatted")
        read(99)grids
    close(99)
    allocate(wfn(NGrids))

    !Discretize momentum space
    NMomenta = floor((pmax - pmin) / dp) + 1
    allocate(momenta(NMomenta))
    momenta(1) = pmin
    do i = 2, NMomenta
        momenta(i) = momenta(i - 1) + dp
    end do
    !Apply fewer space grid points for Wigner distribution
    NWignerGrids = size(grids(1:NGrids:grids_increment))

    !Calculate Wigner distribution
    write(*,*)"Calculating Wigner distribution..."
    allocate(Wigner(NWignerGrids, NMomenta))
    open(unit=99, file="wfn.out", form="unformatted")
    do j = 1, NStates
        write(count,*)j
        open(unit=99+j, file="Wigner"//trim(adjustl(count))//".out", form="unformatted", status="replace")
    end do
    do i = 1, NSnapshots
        do j = 1, NStates
            read(99)wfn
            call wfn2Wigner(grids, momenta, wfn, Wigner, &
            NGrids, NMomenta, NWignerGrids, grids_increment)
            write(99+j)Wigner
        end do
    end do
    close(99)
    do j = 1, NStates
        close(99+j)
    end do

    !Output grid information
    open(unit=99, file="Wigner_q.out", form="unformatted", status="replace")
        write(99)grids(1:NGrids:grids_increment)
    close(99)
    open(unit=99, file="Wigner_p.out", form="unformatted", status="replace")
        write(99)momenta
    close(99)

    write(*,*)
    write(*,*)"Mission success"

contains
!Wigner(q, p) = WignerTransform[wfn(q)]
subroutine wfn2Wigner(grids, momenta, wfn, Wigner, NGrids, NMomenta, NWignerGrids, grids_increment)
    integer, intent(in)::NGrids, NMomenta, NWignerGrids, grids_increment
    real*8, dimension(NGrids)::grids
    real*8, dimension(NMomenta)::momenta
    complex*16, dimension(NGrids), intent(in)::wfn
    real*8, dimension(NWignerGrids, NMomenta), intent(out)::Wigner
    complex*16, parameter::cim2 = (0d0,2d0)
    integer::centre, i, j, k
    real*8::dq
    dq = grids(2) - grids(1)
    !$OMP PARALLEL DO PRIVATE(centre, i, j, k)
    do i = 1, NMomenta
        do j = 1, NWignerGrids
            centre = (j - 1) * grids_increment + 1
            Wigner(j, i) = dble(conjg(wfn(centre)) * wfn(centre))
            do k = 1, min(NGrids - centre, centre - 1)
                Wigner(j, i) = Wigner(j, i) + 2d0 * dble( &
                          + exp(cim2 * momenta(i) * dble(k) * dq) &
                          * conjg(wfn(centre + k)) * wfn(centre - k) )
            end do
        end do
    end do
    !$OMP END PARALLEL DO
    Wigner = Wigner * dq / 3.141592653589793d0
end subroutine wfn2Wigner

end program main