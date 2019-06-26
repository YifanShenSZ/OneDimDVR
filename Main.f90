!Numerically solve the time-dependent Schrodinger equation by discrete variable representation method
!IO unit: atomic unit
!Computation unit: atomic unit
program main
    use Basic; use DVR; use Analyzation
    implicit none
    character*32::JobType!Main only accessed input variable
!---------- Initialize ----------
    call ShowTime(); call ReadInput(); call Initialize()
!------------- End --------------

!----------- Run job ------------
    select case(jobtype)
        case('NewTrajectory')
            call NewTrajectory()
        case('TR-p0')
            call TR_p0()
        case('SMD')
            call ComputeSMD()
        case('pRepresentation')
            call pRepresentation()
        case('WignerDistribution')
            call WignerDistribution()
        case default!Throw a warning
            write(*,*)'Program abort: unsupported job type '//trim(adjustl(JobType))
            stop
    end select
!------------- End --------------

!---------- Clean up ------------
    call ShowTime()
    open(unit=99,file='ParametersUsed.out',status='replace')
        write(99,*)NGrid
        write(99,*)dx
        write(99,*)lt
        write(99,*)NState
    close(99)
    write(*,*)'Mission success'
!------------- End --------------

contains
subroutine ReadInput()
    logical::advance
    open(unit=99,file='OneDimDVR.in')
        read(99,*)
        read(99,*)
        read(99,*)
        read(99,*)JobType
            write(*,*)'Job type: '//JobType
        read(99,*)
        read(99,*)NState
        read(99,*)
        read(99,*)mass
        read(99,*)
        read(99,*)x0
        read(99,*)
        read(99,*)left
        read(99,*)
        read(99,*)right
        read(99,*)
        read(99,*)maxdx
        read(99,*)
        read(99,*)maxdt
        read(99,*)
        read(99,*)p0
        read(99,*)
        read(99,*)TotalTime
        read(99,*)
        read(99,*)OutputInterval
        read(99,*)
        read(99,*)ScatteringProblem
        read(99,*)
        read(99,*)pleft
        read(99,*)
        read(99,*)pright
        read(99,*)
        read(99,*)dp
        read(99,*)
        read(99,*)advance
    close(99)
    if(advance) then
        write(*,*)'Advanced input requested, parameters are set to user specification'
        open(unit=99,file='advance.in',status='old')
            namelist /AdvancedInput/ &
                AutoStep,maxtotallength,minpopdev,minpop
            read(99,nml=AdvancedInput)
        close(99)
    end if
end subroutine ReadInput

subroutine ReadTrajectory()
    integer::i,j,k
    real*8::re,im
    open(unit=99,file='ParametersUsed.out')
        read(99,*)NGrid
        read(99,*)dx
        read(99,*)lt
        read(99,*)NState
    close(99)
    allocate(x(NGrid))
    open(unit=99,file='x.out')
        do i=1,NGrid
            read(99,*)x(i)
        end do
    close(99)
    allocate(psy(NGrid,NState,lt))
    open(unit=99,file='Psy.out')
        do k=1,lt
            do j=1,NState
                do i=1,NGrid
                    Read(99,*)re
                    Read(99,*)im
                    psy(i,j,k)=cmplx(re,im)
                end do
            end do
        end do
    close(99)
end subroutine ReadTrajectory

subroutine Initialize()
    select case(jobtype)
        case('NewTrajectory','TR-p0')
        case default
            call ReadTrajectory()
    end select
end subroutine Initialize

end program main