!OneDimDVR: Numerically solve the time-dependent Schrodinger equation by discrete variable representation method
!This is a specialized version for 1 dimensional scatter systems
!Unit: atomic unit
program main
    use General
    use basic
    use solver
    implicit none
    character*128::job

    write(*,'(A)')"OneDimDVR: Numerically solve the time-dependent Schrodinger equation by discrete variable representation method"
    write(*,'(A)')"This is a specialized version for 1 dimensional scatter systems"
    write(*,'(A)')"Yifan Shen 2020"
    write(*,*)
    call ShowTime()
    call read_input()
    call initialize()

    write(*,*)"Propagating wave function in Chebyshev order domain..."
    call propagate_Chebyshev()

    call ShowTime()
    write(*,*)"Mission success"

contains
subroutine read_input()
    open(unit=99, file="Chebyshev.in")
        read(99,*); read(99,*)mass
        read(99,*); read(99,*)order
        read(99,*); read(99,*)left
        read(99,*); read(99,*)right
        read(99,*); read(99,*)dq
        read(99,*); read(99,*)kmin
        read(99,*); read(99,*)kmax
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