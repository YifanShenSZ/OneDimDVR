!OneDimDVR: Numerically solve the time-dependent Schrodinger equation by discrete variable representation method
!Solution is carried out in Chebyshev order domain rather than usual time domain
!This is a specialized version for 1 dimensional scatter systems
!Unit: atomic unit
program main
    use General
    use basic
    use solver
    implicit none
    character*128::job

    write(*,'(A)')"OneDimDVR: Numerically solve the time-dependent Schrodinger equation by discrete variable representation method"
    write(*,'(A)')"Solution is carried out in Chebyshev order domain rather than usual time domain"
    write(*,'(A)')"This is a specialized version for 1 dimensional scatter systems"
    write(*,'(A)')"Yifan Shen 2020"
    write(*,*)
    call ShowTime()
    call read_input()
    call initialize()

    write(*,*)"Propagating wave function in Chebyshev order domain..."
    call propagate_Chebyshev()
    write(*,*)

    call ShowTime()
    write(*,*)"Mission success"

contains
subroutine read_input()
    integer::i
    open(unit=99, file="OneDimDVR-Chebyshev.in", status="old", iostat=i)
    if (i == 0) then
        read(99,*); read(99,*)mass
        read(99,*); read(99,*)order
        read(99,*); read(99,*)left
        read(99,*); read(99,*)right
        read(99,*); read(99,*)dq
        read(99,*); read(99,*)Hmin
        read(99,*); read(99,*)Hmax
        read(99,*); read(99,*)kmin
    else
        close(99)
        open(unit=99, file="OneDimDVR-Chebyshev.in", status="replace")
            write(99,*)"Mass:"
            write(99,*)
            write(99,*)"Order of Chebyshev propagation:"
            write(99,*)
            write(99,*)"Left boundary:"
            write(99,*)
            write(99,*)"Right boundary:"
            write(99,*)
            write(99,*)"Grid spacing:"
            write(99,*)
            write(99,*)"Energy lower bound: (set both to 0 to let OneDimDVR estimate it)"
            write(99,*)
            write(99,*)"Energy upper bound: (set both to 0 to let OneDimDVR estimate it)"
            write(99,*)
            write(99,*)"Minimum wave number to be absorbed:"
            write(99,*)
        close(99)
        stop "Please fill in OneDimDVR-Chebyshev.in, a template has been provided"
    end if
    close(99)
end subroutine read_input

subroutine initialize()
    integer, external::omp_get_max_threads
    if (NStates < omp_get_max_threads()) then
        call omp_set_num_threads(NStates)
    end if
    call initialize_libHd()
    call initialize_libwfn()
end subroutine initialize

end program main