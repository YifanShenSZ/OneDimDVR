!OneDimDVR: Numerically solve the time-dependent Schrodinger equation by discrete variable representation method
!This is a specialized version for 1 dimensional scatter systems
!Unit: atomic unit
program main
    use General
    use basic
    use solver
    implicit none
    character*128::job, representation

    write(*,'(A)')"OneDimDVR: Numerically solve the time-dependent Schrodinger equation by discrete variable representation method"
    write(*,'(A)')"This is a specialized version for 1 dimensional scatter systems"
    write(*,'(A)')"Yifan Shen 2020"
    write(*,*)
    call ShowTime()
    call read_input()
    call initialize()

    select case(job)
    case ("wavefunction")
        write(*,*)"Propagating wave function"
        select case(representation)
        case("diabatic")
            write(*,*)"in diabatic representation"
            call propagate_wavefunction()
        case("adiabatic")
            write(*,*)"in adiabatic representation"
            call adiabatic_propagate_wavefunction()
        case default
            stop "Unknown representation"
        end select
    case("transmission")
        write(*,*)"Calculating transmission and reflection"
        select case(representation)
        case("diabatic")
            write(*,*)"in diabatic representation"
            call transmit_reflect()
        case("adiabatic")
            write(*,*)"in adiabatic representation"
            call adiabatic_transmit_reflect()
        case default
            stop "Unknown representation"
        end select
    case default
        stop "Unknown job type"
    end select
    write(*,*)

    call ShowTime()
    write(*,*)"Mission success"

contains
subroutine read_input()
    integer::i
    open(unit=99, file="OneDimDVR-scatter.in", status="old", iostat=i)
    if (i == 0) then
        read(99,*); read(99,*)job
        read(99,*); read(99,*)representation
        read(99,*); read(99,*)mass
        read(99,*); read(99,*)total_time
        read(99,*); read(99,*)dt
        read(99,*); read(99,*)output_interval
        read(99,*); read(99,*)left
        read(99,*); read(99,*)right
        read(99,*); read(99,*)dq
        read(99,*); read(99,*)kmin
        read(99,*); read(99,*)kmax
    else
        close(99)
        open(unit=99, file="OneDimDVR-scatter.in", status="replace")
            write(99,*)"Job type: (wavefunction, transmission)"
            write(99,*)
            write(99,*)"Representation: (diabatic, adiabatic)"
            write(99,*)
            write(99,*)"Mass:"
            write(99,*)
            write(99,*)"Total propagation time:"
            write(99,*)
            write(99,*)"Time step:"
            write(99,*)
            write(99,*)"Output interval:"
            write(99,*)
            write(99,*)"Left boundary:"
            write(99,*)
            write(99,*)"Right boundary:"
            write(99,*)
            write(99,*)"Grid spacing:"
            write(99,*)
            write(99,*)"Minimum wave number to be absorbed:"
            write(99,*)
            write(99,*)"Maximum wave number to be absorbed:"
            write(99,*)
        close(99)
        stop "Please fill in OneDimDVR-scatter.in, a template has been provided"
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