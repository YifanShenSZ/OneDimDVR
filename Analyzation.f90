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
    integer::lp,lwx
    real*8,allocatable,dimension(:)::p
    real*8,allocatable,dimension(:,:,:,:)::wigner
    complex*16,allocatable,dimension(:,:,:)::phi

contains
subroutine ComputeSMD()
    integer::i,j,k,m,n,mm,nn,nSMD
    real*8::xave,pave,sigmax,sigmap
    real*8,allocatable,dimension(:)::fct
    real*8,allocatable,dimension(:,:)::pop
    complex*16,allocatable,dimension(:,:,:)::MatrixForm
    type(d2PArray),allocatable,dimension(:)::SMDTemp
    type(SMD4PArray),allocatable,dimension(:)::SMD
    !prepare
        nSMD=SMDOrder*(SMDOrder+3)/2
        allocate(MatrixForm(lx,lx,nSMD))
        call SMDMatrix(lx,nSMD,MatrixForm)
        allocate(fct(0:SMDOrder))
        do i=0,SMDOrder
            fct(i)=dFactorial(i)
        end do
        allocate(SMDTemp(SMDOrder))
        do i=1,SMDOrder
            allocate(SMDTemp(i).Array(0:i))
        end do
        allocate(pop(NState,lt))
        allocate(SMD(lt))
        do i=1,lt
            allocate(SMD(i).NState(NState))
            do j=1,NState
                allocate(SMD(i).NState(j).Order(SMDOrder))
                do k=1,SMDOrder
                    allocate(SMD(i).NState(j).Order(k).Array(0:k))
                end do
            end do
        end do
    do k=1,lt
        do j=1,NState
            pop(j,k)=dot_product(psy(:,j,k),psy(:,j,k))*dx
            m=1
            do i=1,SMDOrder
                do n=0,i
                    SMD(k).NState(j).Order(i).Array(n)=dot_product(psy(:,j,k),matmul(matrixform(:,:,m),psy(:,j,k)))*dx
                    m=m+1
                end do
            end do
            !x, p, xx, xp, pp have been obtained, transform to dimensionless central moments
            xave=SMD(k).NState(j).order(1).Array(0)/pop(j,k); SMD(k).NState(j).order(1).Array(0)=xave
            pave=SMD(k).NState(j).order(1).Array(1)/pop(j,k); SMD(k).NState(j).order(1).Array(1)=pave
            sigmax=Sqrt(SMD(k).NState(j).order(2).Array(0)/pop(j,k)-xave**2); SMD(k).NState(j).order(2).Array(0)=sigmax
            sigmap=Sqrt(SMD(k).NState(j).order(2).Array(2)/pop(j,k)-pave**2); SMD(k).NState(j).order(2).Array(2)=sigmap
            SMD(k).NState(j).order(2).Array(1)=(SMD(k).NState(j).order(2).Array(1)/pop(j,k)-xave*pave)/sigmax/sigmap
            do i=3,SMDOrder
                do n=0,i
                    m=i-n
                    SMDTemp(i).Array(n)=xave**m*pave**n
                    do nn=1,n
                        SMDTemp(i).Array(n)=SMDTemp(i).Array(n)+dCombination(n,nn)*xave**m*pave**(n-nn)*SMD(k).NState(j).order(nn).Array(nn)
                    end do
                    do mm=1,m
                        do nn=0,n
                            SMDTemp(i).Array(n)=SMDTemp(i).Array(n)+dCombination(m,mm)*dCombination(n,nn)*xave**(m-mm)*pave**(n-nn)*SMD(k).NState(j).order(mm+nn).Array(nn)
                        end do
                    end do
                end do
            end do
            do i=3,SMDOrder
                do n=0,i
                    m=i-n
                    SMD(k).NState(j).order(i).Array(n)=SMDTemp(i).Array(n)/sigmax**m/sigmap**n/fct(m)/fct(n)
                end do
                SMD(k).NState(j).order(i).Array=SMD(k).NState(j).order(i).Array/pop(j,k)
            end do
        end do
    end do
    open(unit=99,file='SMD.out',status='replace')!output
        do k=1,lt
            do j=1,NState
                do i=1,SMDOrder
                    do n=0,i
                        write(99,*)SMD(k).NState(j).order(i).Array(n)
                    end do
                end do
                write(99,*)pop(j,k)
            end do
        end do
    close(99)
end subroutine ComputeSMD

!Get the matrix form of x, p, xx, xp, pp (Higher order p terms have not been derived)
!Note that DVR can't be considered as quantum Hilbert space:
!    < A(t) > = dx * < psy(:,t) | A | psy(:,t) >
subroutine SMDMatrix(lx,nSMD,MatrixForm)
    integer,intent(in)::lx,nSMD
    complex*16,dimension(lx,lx,nSMD),intent(inout)::MatrixForm
    integer::i,j,k
    MatrixForm=(0d0,0d0)
    do i=1,SMDOrder
        k=(i-1)*(i+2)/2+1
        forall(j=1:lx)
            MatrixForm(j,j,k)=x(j)**i
        end forall
    end do
    do i=1,lx
        do j=1,lx
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
    forall(i=1:lx,j=1:lx)
        MatrixForm(i,j,4)=MatrixForm(i,j,2)*(x(i)+x(j))/2d0
    end forall
end subroutine SMDMatrix

subroutine pRepresentation()
    integer::i,j,k
    lp=floor((pright-pleft)/dp)+1; allocate(p(lp)); forall(i=1:lp); p(i)=pleft+(i-1)*dp; end forall
    allocate(phi(lp,NState,lt))
    do j=1,lt
        do i=1,NState
            call psy2phi(psy(:,i,j),lx,phi(:,i,j),lp)
        end do
    end do
    open(unit=99,file='p.out',status='replace')!Output
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

subroutine psy2phi(psy,lx,phi,lp)!phi(p) = FourierTransform[psy(x)]
    integer,intent(in)::lx,lp
    complex*16,dimension(lx),intent(in)::psy
    complex*16,dimension(lp),intent(inout)::phi
    integer::k,i
    do k=1,lp
        phi(k)=exp(-ci/hbar*p(k)*x(1))*psy(1)
        do i=2,lx
            phi(k)=phi(k)+exp(-ci/hbar*p(k)*x(i))*psy(i)
        end do
    end do
    phi=phi*dx/sqrtpim2/dSqrt(hbar)
end subroutine psy2phi

subroutine WignerDistribution()
    integer::i,j,k,l
    lp=floor((pright-pleft)/dp)+1; allocate(p(lp)); forall(i=1:lp); p(i)=pleft+(i-1)*dp; end forall
    if(mod(lx,1+skipx)==0) then
        lwx=lx/(1+skipx)
    else
        lwx=lx/(1+skipx)+1
    end if
    allocate(wigner(lwx,lp,NState,lt)); wigner=0d0
    do k=1,lt
        do i=1,NState
            call psy2Wigner(psy(:,i,k),lx,wigner(:,:,i,k),lwx,lp)
        end do
    end do
    open(unit=99,file='Wigner_x.out',status='replace')!Output
        do i=1,lx,1+skipx; write(99,*)x(i); end do
    close(99)
    open(unit=99,file='Wigner_p.out',status='replace')
        do i=1,lp; write(99,*)p(i); end do
    close(99)
    open(unit=99,file='Wigner.out',status='replace')
        do i=1,lt
            do j=1,NState
                do k=1,lp
                    do l=1,lwx
                        write(99,*)wigner(l,k,j,i)
                    end do
                end do
            end do
        end do
    close(99)
end subroutine WignerDistribution

subroutine psy2Wigner(psy,lx,wigner,lwx,lp)!rho(x,p) = WignerTransform[psy(x)]
    integer,intent(in)::lx,lwx,lp
    complex*16,dimension(lx),intent(in)::psy
    real*8,dimension(lwx,lp),intent(out)::wigner
    integer::position,i,j,k
    do i=1,lp
        do j=1,lwx
            position=(j-1)*(1+skipx)+1; wigner(j,i)=0d0
            do k=max(1-position,position-lx),min(lx-position,position-1)
                wigner(j,i)=wigner(j,i)+exp(2d0*ci/hbar*p(i)*dble(k)*dx)*conjg(psy(position+k))*psy(position-k)
            end do
        end do
    end do
    wigner=wigner*dx/pi/hbar
end subroutine psy2Wigner

end module Analyzation