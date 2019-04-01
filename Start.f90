!numerically solve the time dependent schrodinger equation, method:
!DVR: Colbert Miller 1992, Absorb: Manolopoulos 2004
!supports: jobtype: 'NewTrajectory', 'TR-p0', 'QHD', 'pRepresentation', 'WignerRepresentation'
!          steptype: 'Auto', 'Manual'
!atomic unit
program main
    use General
    use DVR
    use V_Psy0
    implicit none
    integer::i,j,k,l,m,n,iorder,jorder,mm,nn,half,nQHD
    real*8::xaveFormer,paveForme,sigmaxForme,sigmapForme,xaveLatter,paveLatter,sigmaxLatter,sigmapLatter,xave,pave,sigmax,sigmap
    real*8,allocatable,dimension(:,:)::Transmission,Reflection
    complex*16,allocatable,dimension(:,:,:)::MatrixForm
    type(double2p),allocatable,dimension(:)::QHDTempFormer,QHDTempLatter,QHDTemp
    type(double4p),allocatable,dimension(:)::QHDFormer,QHDLatter,QHD
    call ShowTime()
    call ReadInput()
    select case(jobtype)
    case('NewTrajectory')
        call solve()
            write(*,*)'actualtime=',actualtime
            lt=floor(actualtime/dt)+1        
            allocate(t(lt))
                forall(i=1:lt)
                    t(i)=(i-1)*dt
                end forall
        open(unit=99,file='x.DVR',status='replace')
            do i=1,lx
                write(99,*)x(i)
            end do
        close(99)
        open(unit=99,file='t.DVR',status='replace')
            do i=1,lt
                write(99,*)t(i)
            end do
        close(99)
        open(unit=99,file='Psy.DVR',status='replace')
            do k=1,surface
                do i=1,lt
                    do j=1,lx
                        write(99,*)real(psy(j,i,k))
                        write(99,*)imag(psy(j,i,k))
                    end do
                end do
            end do
        close(99)
        deallocate(x)
        deallocate(t)
        deallocate(psy)
    case('TR-p0')
        lp0=floor((p0right-p0left)/dp0)+1
        allocate(p0scan(lp0))
            forall(i=1:lp0)
                p0scan(i)=p0left+(i-1)*dp0
            end forall
        allocate(Transmission(surface,lp0))
        allocate(Reflection(surface,lp0))
        do i=1,lp0
            p0=p0scan(i)
            call showtime()
            write(*,*)'p0=',p0
            call scan_solve(surface,Transmission(:,i),Reflection(:,i))
            write(*,*)'actualtime=',actualtime
            deallocate(x)
        end do
        open(unit=99,file='TR.DVR',status='replace')
            do i=1,lp0
                write(99,*)p0scan(i)
                do j=1,surface
                    write(99,*)Transmission(j,i)
                    write(99,*)Reflection(j,i)
                end do
            end do
        close(99)
        deallocate(p0scan)
        deallocate(Transmission)
        deallocate(Reflection)
    case('QHD')
        call ReadTrajectory()
        nQHD=QHDOrder*(QHDOrder+3)/2+1
        allocate(MatrixForm(lx,lx,nQHD-1))
        call QHDMatrix(nQHD-1,MatrixForm)
        do i=1,lx
            if(x(i)>=0d0) then
                exit
            end if
        end do
        half=i
        allocate(QHDFormer(surface))
        allocate(QHDLatter(surface))
        allocate(QHD(surface))
        do i=1,surface
            allocate(QHDFormer(i).time(lt))
            allocate(QHDLatter(i).time(lt))
            allocate(QHD(i).time(lt))
            do j=1,lt
                allocate(QHDFormer(i).time(j).order(QHDOrder+1))
                allocate(QHDLatter(i).time(j).order(QHDOrder+1))
                allocate(QHD(i).time(j).order(QHDOrder+1))
                do k=1,QHDOrder
                    allocate(QHDFormer(i).time(j).order(k).CM(k+1))
                    allocate(QHDLatter(i).time(j).order(k).CM(k+1))
                    allocate(QHD(i).time(j).order(k).CM(k+1))
                end do
                allocate(QHDFormer(i).time(j).order(QHDOrder+1).CM(1))
                allocate(QHDLatter(i).time(j).order(QHDOrder+1).CM(1))
                allocate(QHD(i).time(j).order(QHDOrder+1).CM(1))
            end do
        end do
        allocate(QHDTempFormer(QHDOrder))
        allocate(QHDTempLatter(QHDOrder))
        allocate(QHDTemp(QHDOrder))
        do iorder=1,QHDOrder
            allocate(QHDTemp(i).CM(iorder+1))
            allocate(QHDTempFormer(i).CM(iorder+1))
            allocate(QHDTempLatter(i).CM(iorder+1))
        end do
        do k=1,surface
            do j=1,lt
                i=1
                do iorder=1,QHDOrder
                    do jorder=1,iorder+1
                        QHDFormer(k).time(j).order(iorder).CM(jorder)=dot_product(matmul(matrixform(1:half,1:half,i),psy(1:half,j,k)),psy(1:half,j,k))*dx
                        QHDLatter(k).time(j).order(iorder).CM(jorder)=dot_product(matmul(matrixform(half+1:lx,half+1:lx,i),psy(half+1:lx,j,k)),psy(half+1:lx,j,k))*dx
                        QHD(k).time(j).order(iorder).CM(jorder)=dot_product(matmul(matrixform(:,:,i),psy(:,j,k)),psy(:,j,k))*dx
                    i=i+1
                    end do
                end do
                !transform to dimensionless central moments
                    xaveFormer=QHDFormer(k).time(j).order(1).CM(1)
                    paveFormer=QHDFormer(k).time(j).order(1).CM(2)
                    sigmaxFormer=Sqrt(QHDFormer(k).time(j).order(2).CM(1))
                    sigmapFormer=Sqrt(QHDFormer(k).time(j).order(2).CM(3))
                    xaveLatter=QHDLatter(k).time(j).order(1).CM(1)
                    paveLatter=QHDLatter(k).time(j).order(1).CM(2)
                    sigmaxLatter=Sqrt(QHDLatter(k).time(j).order(2).CM(1))
                    sigmapLatter=Sqrt(QHDLatter(k).time(j).order(2).CM(3))
                    xave=QHD(k).time(j).order(1).CM(1)
                    pave=QHD(k).time(j).order(1).CM(2)
                    sigmax=Sqrt(QHD(k).time(j).order(2).CM(1))
                    sigmap=Sqrt(QHD(k).time(j).order(2).CM(3))
                    do iorder=1,QHDOrder
                        do jorder=1,iorder+1
                            QHDTempFormer(iorder).CM(jorder)=0d0
                            QHDTempLatter(iorder).CM(jorder)=0d0
                            QHDTemp(iorder).CM(jorder)=0d0
                            m=iorder-jorder+1
                            n=jorder-1
                            do mm=0,m
                                do nn=0,n
                                    QHDTempFormer(iorder).CM(jorder)=QHDTempFormer(iorder).CM(jorder)+Combination(m,mm)*Combination(n,nn)*xaveFormer**(m-mm)*paveFormer**(n-nn)*QHDFormer(k).time(j).order(mm+nn).CM(nn+1)
                                    QHDTempLatter(iorder).CM(jorder)=QHDTempLatter(iorder).CM(jorder)+Combination(m,mm)*Combination(n,nn)*xaveLatter**(m-mm)*paveLatter**(n-nn)*QHDLatter(k).time(j).order(mm+nn).CM(nn+1)
                                    QHDTemp(iorder).CM(jorder)=QHDTemp(iorder).CM(jorder)+Combination(m,mm)*Combination(n,nn)*xave**(m-mm)*pave**(n-nn)*QHD(k).time(j).order(mm+nn).CM(nn+1)
                                end do
                            end do
                        end do
                    end do
                    do iorder=1,QHDOrder
                        do jorder=1,iorder+1
                            QHDFormer(k).time(j).order(iorder).CM(jorder)=QHDTempFormer(iorder).CM(jorder)/sigmaxFormer**(iorder-jorder+1)/sigmapFormer**(jorder-1)
                            QHDLatter(k).time(j).order(iorder).CM(jorder)=QHDTempLatter(iorder).CM(jorder)/sigmaxLatter**(iorder-jorder+1)/sigmapLatter**(jorder-1)
                            QHD(k).time(j).order(iorder).CM(jorder)=QHDTemp(iorder).CM(jorder)/sigmax**(iorder-jorder+1)/sigmap**(jorder-1)
                        end do
                    end do
                QHDFormer(k).surface(k).order(QHDOrder+1).CM(1)=dot_product(psy(1:half,j,k),psy(1:half,j,k))*dx
                QHDLatter(k).surface(k).order(QHDOrder+1).CM(1)=dot_product(psy(half+1:lx,j,k),psy(half+1:lx,j,k))*dx
                QHD(j).surface(k).order(QHDOrder+1).CM(1)=QHDFormer(nQHD,j,k)+QHDLatter(nQHD,j,k)*dx
            end do
        end do
        open(unit=99,file='QHDFormer.DVR',status='replace')
            do k=1,surface
                do i=1,lt
                    do j=1,nQHD
                        write(99,*)QHDFormer(j,i,k)
                    end do
                end do
            end do
        close(99)
        open(unit=99,file='QHDLatter.DVR',status='replace')
            do k=1,surface
                do i=1,lt
                    do j=1,nQHD
                        write(99,*)QHDLatter(j,i,k)
                    end do
                end do
            end do
        close(99)
        open(unit=99,file='QHD.DVR',status='replace')
            do k=1,surface
                do i=1,lt
                    do j=1,nQHD
                        write(99,*)QHD(j,i,k)
                    end do
                end do
            end do
        close(99)
        deallocate(QHDFormer)
        deallocate(QHDLatter)
        deallocate(QHD)
        deallocate(matrixform)
    case('pRepresentation')
        call ReadTrajectory()
        allocate(p0scan(NGrid))
            dp0=hbar*pim2/(right-left)
            forall(i=1:NGrid)
                p0scan(i)=i*dp0
            end forall
            p0scan=p0scan-dp0*NGrid/2d0
        allocate(phi(NGrid,lt,surface))
        do k=1,surface
            do i=1,lt
                call Transform2p(psy(:,i,k),lx,phi(:,i,k),NGrid)
            end do
        end do
        open(unit=99,file='k.DVR',status='replace')
            do i=1,NGrid
                write(99,*)p0scan(i)
            end do
        close(99)
        open(unit=99,file='Phi.DVR',status='replace')
            do k=1,surface
                do i=1,lt
                    do j=1,NGrid
                        write(99,*)real(phi(j,i,k))
                        write(99,*)imag(phi(j,i,k))
                    end do
                end do
            end do
        close(99)
        deallocate(x)
        deallocate(psy)
        deallocate(p0scan)
        deallocate(phi)
    case('WignerRepresentation')
        call ReadTrajectory()
        allocate(p0scan(NGrid))
            open(unit=99,file='k.DVR')
                do i=1,NGrid
                    read(99,*)p0scan(i)
                end do
            close(99)
        allocate(wigner(NGrid,NGrid,lt,surface))
        wigner=0d0
        do k=1,surface
            do i=1,lt
                call Transform2Wigner(psy(:,i,k),lx,NGrid,wigner)
            end do
        end do
        open(unit=99,file='Wigner.DVR',status='replace')
            do k=1,surface
                do i=1,lt
                    do j=1,NGrid
                        do l=1,lx
                            write(99,*)wigner(l,j,i,k)
                        end do
                    end do
                end do
            end do
        close(99)
        deallocate(x)
        deallocate(psy)
        deallocate(p0scan)
        deallocate(wigner)
    case default
        write(*,*)'Unkown job type'
    end select
    open(unit=99,file='ParametersUsed.DVR',status='replace')
        write(99,*)NGrid
        write(99,*)NAbsorbGrid
        write(99,*)lx
        write(99,*)dx
        write(99,*)actualtime
        write(99,*)lt
        write(99,*)dt
        write(99,*)lp0
        write(99,*)surface
    close(99)
    write(*,*)'Mission success'
end program main