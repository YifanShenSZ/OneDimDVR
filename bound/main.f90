!OneDimDVR: Numerically solve the time-dependent Schrodinger equation by discrete variable representation method
!This is a specialized version for 1 dimensional bound systems
!Unit: atomic unit
program main
    use General
    use libHd
    use libwfn
    use solver
    implicit none

    write(*,*)"OneDimDVR: Numerically solve the time-dependent Schrodinger equation by discrete variable representation method"
    write(*,*)"This is a specialized version for 1 dimensional bound systems"
    write(*,*)"Yifan Shen 2020"
    write(*,*)
    call ShowTime()
    call ReadInput()
    call Initialize()
    write(*,*)

    write(*,*)"Propagating wave function..."
    call propagate_wavefunction()
    write(*,*)

    call ShowTime()
    write(*,*)"Mission success"

contains
subroutine ReadInput()
    open(unit=99, file="OneDimDVR-bound.in")
        read(99,*); read(99,*)NStates
        read(99,*); read(99,*)mass
        read(99,*); read(99,*)total_time
        read(99,*); read(99,*)output_interval
        read(99,*); read(99,*)left
        read(99,*); read(99,*)right
        read(99,*); read(99,*)dq
    close(99)
end subroutine ReadInput

subroutine Initialize()
    call initialize_libHd()
    call initialize_libwfn()
end subroutine Initialize

end program main