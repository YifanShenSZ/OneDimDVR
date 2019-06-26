module Analyzation
    use Basic
    use DVR
    implicit none

!Parameter
    integer,parameter::SMDOrder=2!Higher order SMD operators in DVR have not been derived

!Derived type
    type SMD3PArray!For SMD quantity calculation
        type(d2PArray),allocatable,dimension(:)::Order
    end type SMD3PArray
    type SMD4PArray
        type(SMD3PArray),allocatable,dimension(:)::NState
    end type SMD4PArray

!Global variable
    integer::lp
    real*8,allocatable,dimension(:)::p
    real*8,allocatable,dimension(:,:,:,:)::wigner
    complex*16,allocatable,dimension(:,:,:)::phi

contains
subroutine ComputeSMD()
    integer::i,j,k,m,n,iorder,jorder,mm,nn,nSMD
    real*8::xave,pave,sigmax,sigmap
    real*8,allocatable,dimension(:)::fct
    real*8,allocatable,dimension(:,:)::pop
    complex*16,allocatable,dimension(:,:,:)::MatrixForm
    type(d2PArray),allocatable,dimension(:)::SMDTemp
    type(SMD4PArray),allocatable,dimension(:)::SMD
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
    open(unit=99,file='SMD.out',status='replace')!output
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
end subroutine ComputeSMD

!Get the matrix form of x, p, xx, xp, pp (Higher order p terms have not been derived)
!Note that DVR can't be considered as quantum Hilbert space:
!    < A(t) > = dx * < psy(:,t) | A | psy(:,t) >
subroutine SMDMatrix(NGrid,nSMD,MatrixForm)
    integer,intent(in)::NGrid,nSMD
    complex*16,dimension(NGrid,NGrid,nSMD),intent(inout)::MatrixForm
    integer::i,j,k
    MatrixForm=(0d0,0d0)
    do i=1,SMDOrder
        k=(i-1)*(i+2)/2+1
        forall(j=1:NGrid)
            MatrixForm(j,j,k)=x(j)**i
        end forall
    end do
    do i=1,NGrid
        do j=1,NGrid
            if(i==j) then
                MatrixForm(i,j,2)=0d0
                MatrixForm(i,j,5)=pisqd3
            else
                k=i-j
                MatrixForm(i,j,2)=1d0/k
                MatrixForm(i,j,5)=2d0/k**2
                if(mod(k,2)) then
                    MatrixForm(i,j,2)=-MatrixForm(i,j,2)
                    MatrixForm(i,j,5)=-MatrixForm(i,j,5)
                end if
            end if
        end do
    end do
    MatrixForm(:,:,2)=-ci*hbar/dx*MatrixForm(:,:,2)
    MatrixForm(:,:,5)=hbar**2/dx**2*MatrixForm(:,:,5)
    forall(i=1:NGrid,j=1:NGrid)
        MatrixForm(i,j,4)=MatrixForm(i,j,2)*(x(i)+x(j))/2d0
    end forall
end subroutine SMDMatrix

subroutine pRepresentation()
    integer::i,j,k
    lp=floor((pright-pleft)/dp)+1!Prepare
    allocate(p(lp))
    forall(i=1:lp)
        p(i)=pleft+(i-1)*dp
    end forall
    allocate(phi(lp,NState,lt))
    do j=1,lt
        do i=1,NState
            call Transform2p(psy(:,i,j),phi(:,i,j))
        end do
    end do
    open(unit=99,file='k.out',status='replace')!Output
        do i=1,lp
            write(99,*)p(i)
        end do
    close(99)
    open(unit=99,file='Phi.out',status='replace')
        do k=1,lt
            do j=1,NState
                do i=1,lp
                    write(99,*)real(phi(i,j,k))
                    write(99,*)imag(phi(i,j,k))
                end do
            end do
        end do
    close(99)
end subroutine pRepresentation

!phi(p) = Fourier[psy(x)]
subroutine Transform2p(psy,phi)
    complex*16,dimension(NGrid),intent(in)::psy
    complex*16,dimension(lp),intent(inout)::phi
    integer::k,i
    do k=1,lp
        phi(k)=exp(-ci/hbar*p(k)*x(1))*psy(1)
        do i=2,NGrid
            phi(k)=phi(k)+exp(-ci/hbar*p(k)*x(i))*psy(i)
        end do
    end do
    phi=phi*dx/sqrtpim2/dSqrt(hbar)
end subroutine Transform2p

subroutine WignerDistribution()
    integer::i,j,k,l
    lp=floor((pright-pleft)/dp)+1!Prepare
    allocate(p(lp))
    forall(i=1:lp)
        p(i)=pleft+(i-1)*dp
    end forall
    allocate(wigner(NGrid,lp,NState,lt))
    wigner=0d0
    do k=1,lt
        do i=1,NState
            call Transform2Wigner(psy(:,i,k),NGrid,NGrid,wigner(:,:,i,k))
        end do
    end do
    open(unit=99,file='Wigner.out',status='replace')!Output
        do i=1,lt
            do j=1,NState
                do k=1,lp
                    do l=1,NGrid
                        write(99,*)wigner(l,k,j,i)
                    end do
                end do
            end do
        end do
    close(99)
end subroutine WignerDistribution

!p(x,p) = WignerTransform[psy(x)]
subroutine Transform2Wigner(psy,lx,lk,wigner)
    integer::lx,lk,k,i,j,width
    real*8,dimension(lx)::ones
    real*8,dimension(lx,lk)::wigner
    complex*16,dimension(lx)::psy,temp
    ones=1d0
    do k=1,lk
        do j=2,lx-1
            temp=(0d0,0d0)
            width=min(j,lx-j)-1
            forall(i=-width:width)
                temp(i+width+1)=exp(2d0*ci/hbar*p(k)*i*dx)*psy(j+i)*psy(j-i)
            end forall
            wigner(j,k)=dot_product(ones,temp)
        end do
    end do
    wigner=wigner*dx/pi/hbar
end subroutine Transform2Wigner

end module Analyzation