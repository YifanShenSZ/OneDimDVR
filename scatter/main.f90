!OneDimDVR: Numerically solve the time-dependent Schrodinger equation by discrete variable representation method
!This is a specialized version for 1 dimensional scatter systems
!Unit: atomic unit
program main
    use General
    use basic
    use solver
    implicit none
    character*128::job

    write(*,*)"OneDimDVR: Numerically solve the time-dependent Schrodinger equation by discrete variable representation method"
    write(*,*)"This is a specialized version for 1 dimensional scatter systems"
    write(*,*)"Yifan Shen 2020"
    write(*,*)
    call ShowTime()
    call ReadInput()
    call Initialize()

    select case(job)
    case ("wavefunction")
        write(*,*)"Propagating wave function..."
        call propagate_wavefunction()
    case("transmission")
        write(*,*)"Calculating transmission and reflection..."
        call transmit_reflect()
    end select
    write(*,*)

    call ShowTime()
    write(*,*)"Mission success"

contains
subroutine ReadInput()
    open(unit=99, file="OneDimDVR-scatter.in")
        read(99,*); read(99,*)job; write(*,*)"Job type: "//job
        read(99,*); read(99,*)mass
        read(99,*); read(99,*)total_time
        read(99,*); read(99,*)dt
        read(99,*); read(99,*)output_interval
        read(99,*); read(99,*)left
        read(99,*); read(99,*)right
        read(99,*); read(99,*)dq
        read(99,*); read(99,*)kmin
        read(99,*); read(99,*)kmax
    close(99)
end subroutine ReadInput

subroutine Initialize()
    call initialize_libHd()
    call initialize_libwfn()
end subroutine Initialize

end program main