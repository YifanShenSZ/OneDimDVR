!Calculate transmission and reflection from Chebyshev order domain wave function
program main
    implicit none

type pint
    integer, allocatable, dimension(:)::indices
end type pint

    !input
    real*8::mass, dt, left, right
    integer::left_index, right_index
    !space grid points and Chebyshev order domain wave function
    integer::NStates, order, NGrids
    real*8::Hmin, Hmax, dq
    real*8, allocatable, dimension(:)::grids
    complex*16, allocatable, dimension(:,:,:)::Chebyshev
    !transformation to time domain and transmission & reflection
    real*8::R
    complex*16::scalor
    type(pint), allocatable, dimension(:)::order_left, order_right ! order > 0 and wfn norm > 0 orders
    complex*16, allocatable, dimension(:,:)::p ! DVR momentum
    complex*16, allocatable, dimension(:)::p_left, p_right ! rows in p
    complex*16, allocatable, dimension(:,:)::Chebyshev_left, Chebyshev_right, pCheStar_left, pCheStar_right
    real*8, allocatable, dimension(:)::transmission, reflection
    !parallelization
    integer, external::omp_get_max_threads
    integer::NThreads, thread
    complex*16, allocatable, dimension(:,:)::wfn_left, wfn_right, pwfnstar_left, pwfnstar_right
    real*8, allocatable, dimension(:,:)::tran_para, refl_para
    !work variable
    integer::i, j, k, count, timer
    real*8::time, dcoeff, tol
    complex*16::zcoeff

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

    !Initialize time domain convertion and transmission & reflection
    allocate(order_left (NStates))
    allocate(order_right(NStates))
    do j = 1, NStates
        !Only wfn norm > 1d-8 orders are taken into account
        count = 0
        do i = 1, order
            if (abs(Chebyshev(left_index,j,i)) > 1d-8) count = count + 1
        end do
        allocate(order_left(j)%indices(count))
        count = 0
        do i = 1, order
            if (abs(Chebyshev(left_index,j,i)) > 1d-8) then
                count = count + 1
                order_left(j)%indices(count) = i
            end if
        end do
        count = 0
        do i = 1, order
            if (abs(Chebyshev(right_index,j,i)) > 1d-8) count = count + 1
        end do
        allocate(order_right(j)%indices(count))
        count = 0
        do i = 1, order
            if (abs(Chebyshev(right_index,j,i)) > 1d-8) then
                count = count + 1
                order_right(j)%indices(count) = i
            end if
        end do
    end do
    allocate(p(NGrids, NGrids))
    call compute_momentum(dq, p, NGrids)
    allocate(p_left(NGrids))
    p_left = p(left_index, :)
    allocate(p_right(NGrids))
    p_right = p(right_index, :)
    allocate(Chebyshev_left (NStates, 0:order))
    allocate(Chebyshev_right(NStates, 0:order))
    forall(i = 0 : order)
        Chebyshev_left (:,i) = Chebyshev( left_index,:,i)
        Chebyshev_right(:,i) = Chebyshev(right_index,:,i)
    end forall
    allocate(pCheStar_left (NStates, 0:order))
    allocate(pCheStar_right(NStates, 0:order))
    forall(i = 0 : order, j = 1 : NStates)
        pCheStar_left (j,i) = dot_product(Chebyshev(:,j,i), p_left )
        pCheStar_right(j,i) = dot_product(Chebyshev(:,j,i), p_right)
    end forall
    allocate(transmission(NStates))
    transmission = 0d0
    allocate(reflection  (NStates))
    reflection   = 0d0

    !Initialize parallelization
    NThreads = omp_get_max_threads()
    allocate(wfn_left (NStates, NThreads))
    allocate(wfn_right(NStates, NThreads))
    allocate(pwfnstar_left (NStates, NThreads))
    allocate(pwfnstar_right(NStates, NThreads))
    allocate(tran_para(NStates, NThreads))
    allocate(refl_para(NStates, NThreads))

    !Calculate transmission and reflection
    timer = 0
    tol   = 0.9999 / dt * mass
    do
        tran_para = 0d0
        refl_para = 0d0
        !$OMP PARALLEL DO PRIVATE(thread, time, count, scalor, R, dcoeff, j, i, k, zcoeff)
        do thread = 1, NThreads
            !Take 1000 time steps per thread
            time = (timer * NThreads * 1000 + (thread - 1) * 1000) * dt
            do count = 1, 1000
                !Transform to time domain and calculate transmission & reflection
                time = time + dt
                scalor = exp((0d0,-1d0) * (Hmax + Hmin) * time / 2d0)
                R = (Hmax - Hmin) * time / 2d0
                dcoeff = bessel_j0(R)
                wfn_left (:,thread) = dcoeff * Chebyshev_left (:,0)
                wfn_right(:,thread) = dcoeff * Chebyshev_right(:,0)
                pwfnstar_left (:,thread) = dcoeff * pCheStar_left (:,0)
                pwfnstar_right(:,thread) = dcoeff * pCheStar_right(:,0)
                do j = 1, NStates
                    do i = 1, size(order_left(j)%indices)
                        k = order_left(j)%indices(i)
                        dcoeff = 2d0 * bessel_jn(k, R)
                        !Only |coeff| > 1d-8 terms are taken into account
                        if (abs(dcoeff) > 1d-8) then
                            zcoeff = (0d0,1d0)**k * dcoeff
                            wfn_left(j,thread) = wfn_left(j,thread) + zcoeff * Chebyshev_left(j,k)
                            zcoeff = conjg(zcoeff)
                            pwfnstar_left(j,thread) = pwfnstar_left(j,thread) + zcoeff * pCheStar_left(j,k)
                        end if
                    end do
                    do i = 1, size(order_right(j)%indices)
                        k = order_right(j)%indices(i)
                        dcoeff = 2d0 * bessel_jn(k, R)
                        !Only |coeff| > 1d-8 terms are taken into account
                        if (abs(dcoeff) > 1d-8) then
                            zcoeff = (0d0,1d0)**k * dcoeff
                            wfn_right(j,thread) = wfn_right(j,thread) + zcoeff * Chebyshev_right(j,k)
                            zcoeff = conjg(zcoeff)
                            pwfnstar_right(j,thread) = pwfnstar_right(j,thread) + zcoeff * pCheStar_right(j,k)
                        end if
                    end do
                end do
                !Calculate transmission and reflection
                tran_para(:,thread) = tran_para(:,thread) - dble(wfn_right(:,thread) * pwfnstar_right(:,thread))
                refl_para(:,thread) = refl_para(:,thread) + dble(wfn_left (:,thread) * pwfnstar_left (:,thread))
            end do
        end do
        !$OMP END PARALLEL DO
        transmission = transmission + sum(tran_para, 2)
        reflection   = reflection   + sum(refl_para, 2)
        timer = timer + 1
        write(*,*)"time = ", timer * NThreads * 1000 * dt
        write(99,*)"state    transmission    reflection"
        do i = 1, NStates
            write(99,"(I5,4x,f12.8,4x,f12.8)")i, transmission(i) * dt / mass, reflection(i) * dt / mass
        end do
        !Check population
        if (sum(transmission) + sum(reflection) > tol) then
            write(*,*)"All population has been absorbed"
            write(*,*)"Stop propagation at time ", time
            exit
        end if
    end do
    transmission = transmission * dt / mass
    reflection   = reflection   * dt / mass

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