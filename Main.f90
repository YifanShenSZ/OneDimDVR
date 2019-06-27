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
        case('pRepresentation')
            call pRepresentation()
        case('WignerDistribution')
            call WignerDistribution()
        case('SMD')
            call ComputeSMD()
        case default!Throw a warning
            write(*,*)'Program abort: unsupported job type '//trim(adjustl(JobType))
            stop
    end select
!------------- End --------------

!---------- Clean up ------------
    call ShowTime(); write(*,*)'Mission success'
!------------- End --------------

contains
subroutine ReadInput()
    logical::advance
    open(unit=99,file='OneDimDVR.in')
        read(99,*); read(99,*); read(99,*); read(99,*)JobType; write(*,*)'Job type: '//JobType
        read(99,*); read(99,*)NState
        read(99,*); read(99,*)mass
        read(99,*); read(99,*)x0
        read(99,*); read(99,*)left
        read(99,*); read(99,*)right
        read(99,*); read(99,*)maxdx
        read(99,*); read(99,*)maxdt
        read(99,*); read(99,*)p0
        read(99,*); read(99,*)TotalTime
        read(99,*); read(99,*)OutputInterval
        read(99,*); read(99,*)ScatteringProblem
        read(99,*); read(99,*)pleft
        read(99,*); read(99,*)pright
        read(99,*); read(99,*)dp
        read(99,*); read(99,*)skipx
        read(99,*); read(99,*)advance
    close(99)
    if(advance) then
        write(*,*)'Advanced input requested, parameters are set to user specification'
        open(unit=99,file='advance.in',status='old')
            namelist /AdvancedInput/ &
                AutoStep,maxtotallength,maxpopdev,minpop
            read(99,nml=AdvancedInput)
        close(99)
    end if
end subroutine ReadInput

subroutine Initialize()
    integer::i,j,k
    real*8::re,im
    if(jobtype/='NewTrajectory'.and.jobtype/='TR-p0') then!Read old trajectory
        open(unit=99,file='x.out')
            lx=0; do; read(99,*,iostat=i); if(i/=0) exit; lx=lx+1; end do; rewind 99
            allocate(x(lx)); do i=1,lx; read(99,*)x(i); end do
            dx=x(2)-x(1)
        close(99)
        open(unit=99,file='t.out')
            lt=0; do; read(99,*,iostat=i); if(i/=0) exit; lt=lt+1; end do; rewind 99
            allocate(t(lt)); do i=1,lt; read(99,*)t(i); end do
            dt=t(2)-t(1)
        close(99)
        allocate(psy(lx,NState,lt))
        open(unit=99,file='Psy.out')
            do k=1,lt
                do j=1,NState
                    do i=1,lx
                        Read(99,*)re; Read(99,*)im; psy(i,j,k)=cmplx(re,im)
                    end do
                end do
            end do
        close(99)
    end if
end subroutine Initialize

end program main