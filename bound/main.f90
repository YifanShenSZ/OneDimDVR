!OneDimDVR: Numerically solve the time-dependent Schrodinger equation by discrete variable representation method
!This is a specialized version for 1 dimensional bound systems
!Unit: atomic unit
program main
    use General
    use basic
    use solver
    implicit none

    write(*,'(A)')"OneDimDVR: Numerically solve the time-dependent Schrodinger equation by discrete variable representation method"
    write(*,'(A)')"This is a specialized version for 1 dimensional bound systems"
    write(*,'(A)')"Yifan Shen 2020"
    write(*,*)
    call ShowTime()
    call read_input()
    call initialize()
    write(*,*)

    write(*,*)"Propagating wave function..."
    call propagate_wavefunction()
    write(*,*)

    call ShowTime()
    write(*,*)"Mission success"

contains
subroutine read_input()
    integer::i
    open(unit=99, file="OneDimDVR-bound.in", status="old", iostat=i)
    if (i == 0) then
        read(99,*); read(99,*)mass
        read(99,*); read(99,*)total_time
        read(99,*); read(99,*)output_interval
        read(99,*); read(99,*)left
        read(99,*); read(99,*)right
        read(99,*); read(99,*)dq
    else
        close(99)
        open(unit=99, file="OneDimDVR-bound.in", status="replace")
            write(99,*)"Mass:"
            write(99,*)
            write(99,*)"Total propagation time:"
            write(99,*)
            write(99,*)"Output interval:"
            write(99,*)
            write(99,*)"Left boundary:"
            write(99,*)
            write(99,*)"Right boundary:"
            write(99,*)
            write(99,*)"Grid spacing:"
            write(99,*)
        close(99)
        stop "Please fill in OneDimDVR-bound.in, a template has been provided"
    end if
    close(99)
end subroutine read_input

subroutine initialize()
    call initialize_libHd()
    call initialize_libwfn()
end subroutine initialize

end program main