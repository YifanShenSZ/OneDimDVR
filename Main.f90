!Numerically solve the time dependent schrodinger equation, method:
!DVR: Colbert Miller 1992, Absorb: Manolopoulos 2004
!IO unit: input time in ps or fs (specified in input file), others in atomic unit
!Computation unit: atomic unit
program main
    use Basic
    use DVR
    implicit none
    !Derived type for SMD quantity calculation
        type SMD3PArray
            type(d2PArray),allocatable,dimension(:)::Order
        end type SMD3PArray
        type SMD4PArray
            type(SMD3PArray),allocatable,dimension(:)::NState
        end type SMD4PArray
    !Work variable
    integer::i,j,k,l,m,n,iorder,jorder,mm,nn,nSMD
    real*8::xave,pave,sigmax,sigmap
    real*8,allocatable,dimension(:)::fct
    real*8,allocatable,dimension(:,:)::Transmission,Reflection,pop
    complex*16,allocatable,dimension(:,:,:)::MatrixForm
    type(d2PArray),allocatable,dimension(:)::SMDTemp
    type(SMD4PArray),allocatable,dimension(:)::SMD
!---------- Initialize ----------
    call ShowTime()
    call ReadInput()
!------------- End --------------

!----------- Run job ------------
    select case(jobtype)
        case('NewTrajectory')
            if(.not.AutoStep) write(*,*)'Auto dx dt determination is disabled, adopt user specified upper limit'
            if(ScatteringProblem) then
                write(*,*)'Scattering problem, turn on absorbing potential'
                call SolveAbsorb()
                write(*,*)'Actual evolution time =',ActualTime/fsInAU/1000d0,'ps'
            else
                call Solve()
            end if
            !output
                allocate(t(lt))
                forall(i=1:lt)
                    t(i)=(i-1)*OutputInterval
                end forall
                open(unit=99,file='x.DVR',status='replace')
                    do i=NAbsorbGrid+1,NAbsorbGrid+NGrid
                        write(99,*)x(i)
                    end do
                close(99)
                open(unit=99,file='t.DVR',status='replace')
                    do i=1,lt
                        write(99,*)t(i)
                    end do
                close(99)
                open(unit=99,file='Psy.DVR',status='replace')
                    do k=1,lt
                        do j=1,NState
                            do i=NAbsorbGrid+1,NAbsorbGrid+NGrid
                                write(99,*)real(psy(i,j,k))
                                write(99,*)imag(psy(i,j,k))
                            end do
                        end do
                    end do
                close(99)
            !clean up
                deallocate(x)
                deallocate(t)
                deallocate(psy)
        case('TR-p0')
            if(.not.AutoStep) write(*,*)'Auto dx dt determination is disabled, adopt user specified upper limit'
            !prepare
                lp0=floor((p0right-p0left)/dp0)+1
                allocate(p0scan(lp0))
                    forall(i=1:lp0)
                        p0scan(i)=p0left+(i-1)*dp0
                    end forall
                allocate(Transmission(NState,lp0))
                allocate(Reflection(NState,lp0))
            do i=1,lp0
                p0=p0scan(i)
                call ShowTime()
                write(*,*)'p0 =',p0
                call scan_solve(Transmission(:,i),Reflection(:,i))
                write(*,*)'Actual evolution time =',ActualTime
            end do
            !output
                open(unit=99,file='TR.DVR',status='replace')
                    do i=1,lp0
                        write(99,*)p0scan(i)
                        do j=1,NState
                            write(99,*)Transmission(j,i)
                            write(99,*)Reflection(j,i)
                        end do
                    end do
                close(99)
            !clean up
                deallocate(p0scan)
                deallocate(Transmission)
                deallocate(Reflection)
        case('SMD')
            call ReadTrajectory()
            !prepare
                nSMD=SMDOrder*(SMDOrder+3)/2
                allocate(MatrixForm(NGrid,NGrid,nSMD))
                call SMDMatrix(NGrid,nSMD,MatrixForm)
                allocate(fct(SMDOrder+1))
                do i=0,SMDOrder
                    fct(i+1)=dFactorial(i)
                end do
                allocate(SMDTemp(SMDOrder))
                do i=1,SMDOrder
                    allocate(SMDTemp(i).Array(i+1))
                end do
                allocate(pop(NState,lt))
                allocate(SMD(lt))
                do i=1,lt
                    allocate(SMD(i).NState(NState))
                    do j=1,NState
                        allocate(SMD(i).NState(j).Order(SMDOrder+1))
                        do k=1,SMDOrder
                            allocate(SMD(i).NState(j).Order(k).Array(k+1))
                        end do
                    end do
                end do
            do k=1,lt
                do j=1,NState
                    pop(j,k)=dot_product(psy(:,j,k),psy(:,j,k))*dx
                    i=1
                    do iorder=1,SMDOrder
                        do jorder=1,iorder+1
                            SMD(k).NState(j).Order(iorder).Array(jorder)=dot_product(psy(:,j,k),matmul(matrixform(:,:,i),psy(:,j,k)))*dx
                            i=i+1
                        end do
                    end do
                    !transform to dimensionless central moments
                        xave=SMD(k).NState(j).order(1).Array(1)/pop(j,k)
                        pave=SMD(k).NState(j).order(1).Array(2)/pop(j,k)
                        sigmax=Sqrt(SMD(k).NState(j).order(2).Array(1)/pop(j,k)-xave**2)
                        sigmap=Sqrt(SMD(k).NState(j).order(2).Array(3)/pop(j,k)-pave**2)
                        SMD(k).NState(j).order(2).Array(1)=sigmax
                        SMD(k).NState(j).order(2).Array(3)=sigmap
                        SMD(k).NState(j).order(2).Array(2)=(SMD(k).NState(j).order(2).Array(2)/pop(j,k)-xave*pave)/sigmax/sigmap
                        do iorder=3,SMDOrder
                            do jorder=1,iorder+1
                                m=iorder-jorder+1
                                n=jorder-1
                                SMDTemp(iorder).Array(jorder)=xave**m*pave**n
                                do nn=1,n
                                    SMDTemp(iorder).Array(jorder)=SMDTemp(iorder).Array(jorder)+dCombination(n,nn)*xave**m*pave**(n-nn)*SMD(k).NState(j).order(nn).Array(nn+1)
                                end do
                                do mm=1,m
                                    do nn=0,n
                                        SMDTemp(iorder).Array(jorder)=SMDTemp(iorder).Array(jorder)+dCombination(m,mm)*dCombination(n,nn)*xave**(m-mm)*pave**(n-nn)*SMD(k).NState(j).order(mm+nn).Array(nn+1)
                                    end do
                                end do
                            end do
                        end do
                        do iorder=3,SMDOrder
                            do jorder=1,iorder+1
                                m=iorder-jorder+1
                                SMD(k).NState(j).order(iorder).Array(jorder)=SMDTemp(iorder).Array(jorder)/sigmax**m/sigmap**n/fct(m+1)/fct(jorder)
                            end do
                            SMD(k).NState(j).order(iorder).Array=SMD(k).NState(j).order(iorder).Array/pop(j,k)
                        end do
                end do
            end do
            !output
                open(unit=99,file='SMD.DVR',status='replace')
                    do i=1,lt
                        do j=1,NState
                            do iorder=1,SMDOrder
                                do jorder=1,iorder+1
                                    write(99,*)SMD(i).NState(j).order(iorder).Array(jorder)
                                end do
                            end do
                            write(99,*)pop(j,i)
                        end do
                    end do
                close(99)
            !clean up
                deallocate(MatrixForm)
                deallocate(fct)
                do i=1,SMDOrder
                    deallocate(SMDTemp(i).Array)
                end do
                deallocate(SMDTemp)
                deallocate(pop)
                do i=1,lt
                    do j=1,NState
                        do k=1,SMDOrder
                            deallocate(SMD(i).NState(j).Order(k).Array)
                        end do
                        deallocate(SMD(i).NState(j).Order)
                    end do
                    deallocate(SMD(i).NState)
                end do
                deallocate(SMD)
        case('pRepresentation')
            call ReadTrajectory()
            !prepare
                allocate(p0scan(2*NGrid))
                    dp0=hbar*pi/(right-left)
                    forall(i=1:NGrid)
                        p0scan(NGrid+i)=i*dp0
                        p0scan(NGrid-i+1)=-i*dp0
                    end forall
                allocate(phi(2*NGrid,NState,lt))
            do j=1,lt
                do i=1,NState
                    call Transform2p(psy(:,i,j),phi(:,i,j),NGrid)
                end do
            end do
            !output
                open(unit=99,file='k.DVR',status='replace')
                    do i=1,2*NGrid
                        write(99,*)p0scan(i)
                    end do
                close(99)
                open(unit=99,file='Phi.DVR',status='replace')
                    do k=1,lt
                        do j=1,NState
                            do i=1,2*NGrid
                                write(99,*)real(phi(i,j,k))
                                write(99,*)imag(phi(i,j,k))
                            end do
                        end do
                    end do
                close(99)
            !clean up
               deallocate(x)
               deallocate(psy)
               deallocate(p0scan)
               deallocate(phi)
        case('WignerDistribution')
            call ReadTrajectory()
            !prepare
                allocate(p0scan(NGrid))
                    open(unit=99,file='k.DVR')
                        do i=1,NGrid
                            read(99,*)p0scan(i)
                        end do
                    close(99)
                allocate(wigner(NGrid,NGrid,NState,lt))
                wigner=0d0
            do k=1,lt
                do i=1,NState
                    call Transform2Wigner(psy(:,i,k),NGrid,NGrid,wigner(:,:,i,k))
                end do
            end do
            !output
                open(unit=99,file='Wigner.DVR',status='replace')
                    do k=1,NState
                        do i=1,lt
                            do j=1,NGrid
                                do l=1,NGrid
                                    write(99,*)wigner(l,j,i,k)
                                end do
                            end do
                        end do
                    end do
                close(99)
            !clean up
                deallocate(x)
                deallocate(psy)
                deallocate(p0scan)
                deallocate(wigner)
        case default!Throw a warning
            write(*,*)'Program abort: unsupported job type '//trim(adjustl(JobType))
            stop
    end select
!------------- End --------------

!---------- Clean up ------------
    call ShowTime()
    open(unit=99,file='ParametersUsed.DVR',status='replace')
        write(99,*)NGrid
        write(99,*)dx
        write(99,*)ActualTime
        write(99,*)lt
        write(99,*)lp0
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
        read(99,*)TotalTime
            TotalTime=dAbs(TotalTime)*1000d0*fsInAU!Convert to atomic unit
        read(99,*)
        read(99,*)left
        read(99,*)
        read(99,*)right
        read(99,*)
        read(99,*)maxdx
        read(99,*)
        read(99,*)maxdt
            maxdt=dAbs(maxdt)*fsInAU!Convert to atomic unit
        read(99,*)
        read(99,*)p0
        read(99,*)
        read(99,*)OutputInterval
            OutputInterval=dAbs(OutputInterval)*fsInAU!Convert to atomic unit
        read(99,*)
        read(99,*)ScatteringProblem
        read(99,*)
        read(99,*)p0left
        read(99,*)
        read(99,*)p0right
        read(99,*)
        read(99,*)dp0
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
    open(unit=99,file='ParametersUsed.DVR')
        read(99,*)NGrid
        read(99,*)dx
        read(99,*)ActualTime
        read(99,*)lt
        read(99,*)lp0
        read(99,*)NState
    close(99)
    allocate(x(NGrid))
        open(unit=99,file='x.DVR')
            do i=1,NGrid
                read(99,*)x(i)
            end do
        close(99)
    allocate(psy(NGrid,NState,lt))
        open(unit=99,file='Psy.DVR')
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

end program main