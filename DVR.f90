module DVR
    use V_Psy0
    use General
    implicit none

type Double2P
    real*8,allocatable,dimension(:)::CM
end type Double2P

type Double3P
    type(Double2P),allocatable,dimension(:)::order
end type Double3P

type Double4P
    type(Double3P),allocatable,dimension(:)::time
end type Double4P

!parameters
    !to prevent stackoverflow during auto dt
    integer::maxlx=1300
    !stop evolution when population no longer equals to 1
    real*8::minpopdev=1d-6,minpop=1d-6

!global variables
    complex*16,allocatable,dimension(:,:)::hamiltonian,hamiltonian_absorb

contains
!AbsorbedPotential_ij(x)
function AbsorbedPotential(x,i,j)
    integer::i,j
    real*8::x,y
    complex*16::AbsorbedPotential
    AbsorbedPotential=potential(x,i,j)
    if(i==j) then
        if(x>right) then
            y=cs*kmins/(pim2**2)*(x-right)**2
            AbsorbedPotential=AbsorbedPotential-ci*hbar**2*kmins*4d0/mass*((cs+y)/(cs-y)**2-1d0/cs)
        else if(x<left) then
            y=cs*kmins/(pim2**2)*(x-left)**2
            AbsorbedPotential=AbsorbedPotential-ci*hbar**2*kmins*4d0/mass*((cs+y)/(cs-y)**2-1d0/cs)
        end if
    end if
end function AbsorbedPotential

!x is the vector of grid points
!HamiltonianDVR store as a whole surface*size(x) order matrix
!can be considered as consisting of surface*surface blocks, ij block is surface i x surface j 
function HamiltonianDVR(x,dx,surface)
    integer::surface,i,j,k,lx,ii,jj,kk,kkk,totallength
    real*8::dx
    real*8,allocatable,dimension(:)::x
    real*8,allocatable,dimension(:,:)::tran
    complex*16,allocatable,dimension(:,:)::HamiltonianDVR
    lx=size(x)
    allocate(tran(lx,lx))
    totallength=surface*lx
    allocate(HamiltonianDVR(totallength,totallength))
    HamiltonianDVR=(0d0,0d0)
    do i=1,lx
        do j=1,lx
            if(i==j) then
                tran(i,j)=pipid3
            else
                k=i-j
                tran(i,j)=2d0/k**2
                if(mod(k,2)) then
                    tran(i,j)=-tran(i,j)
                end if
            end if
        end do
    end do
    tran=hbar**2/2d0/mass/dx**2*tran
    j=1
    do i=1,surface
        k=j+lx-1
        HamiltonianDVR(j:k,j:k)=tran
        do ii=1,surface
            kk=(i-1)*lx
            kkk=(ii-1)*lx
            do jj=1,lx
                kk=kk+1
                kkk=kkk+1
                HamiltonianDVR(kk,kkk)=HamiltonianDVR(kk,kkk)+potential(x(jj),i,ii)
            end do
        end do
        j=k+1
    end do
end function HamiltonianDVR

function HamiltonianDVR_absorb(x,dx,surface)
    integer::surface,i,j,k,lx,ii,jj,kk,kkk,totallength
    real*8::dx
    real*8,allocatable,dimension(:)::x
    real*8,allocatable,dimension(:,:)::tran
    complex*16,allocatable,dimension(:,:)::HamiltonianDVR_absorb
    lx=size(x)
    allocate(tran(lx,lx))
    totallength=surface*lx
    allocate(HamiltonianDVR_absorb(totallength,totallength))
    HamiltonianDVR_absorb=(0d0,0d0)
    do i=1,lx
        do j=1,lx
            if(i==j) then
                tran(i,j)=pipid3
            else
                k=i-j
                tran(i,j)=2d0/k**2
                if(mod(k,2)) then
                    tran(i,j)=-tran(i,j)
                end if
            end if
        end do
    end do
    tran=hbar**2/2d0/mass/dx**2*tran
    j=1
    do i=1,surface
        k=j+lx-1
        HamiltonianDVR_absorb(j:k,j:k)=tran
        do ii=1,surface
            kk=(i-1)*lx
            kkk=(ii-1)*lx
            do jj=1,lx
                kk=kk+1
                kkk=kkk+1
                HamiltonianDVR_absorb(kk,kkk)=HamiltonianDVR_absorb(kk,kkk)+AbsorbedPotential(x(jj),i,ii)
            end do
        end do
        j=k+1
    end do
end function HamiltonianDVR_absorb

subroutine evolve(d,u,dim)
    integer::dim
    complex*16,dimension(dim)::d,u
    d=matmul(hamiltonian,u)
end subroutine evolve

subroutine evolve_absorb(d,u,dim)
    integer::dim
    complex*16,dimension(dim)::d,u
    d=matmul(hamiltonian_absorb,u)
end subroutine evolve_absorb

!numerically solve the evolution of a psy0, and save the trajectory
subroutine solve()
    logical::success
    integer::i,j,k,totallength
    real*8::al
    complex*16,allocatable,dimension(:)::eigval
    complex*16,allocatable,dimension(:,:)::psyevolve,eigvec
    call initialize()
    if(kminabs0) then
        stop 'parameter error: min(|k|)=0'
    end if
    if(steptype=='Auto') then
        dx=min(pim2/5d0/kmax,maxdx)
        call dScientificNotation(dx,i)
        dx=floor(dx)*10d0**i
    end if
    NGrid=floor((right-left)/dx)+1
    dx=(right-left)/(NGrid-1)
    write(*,*)'dx=',dx
    al=pim2/kmin
    NAbsorbGrid=ceiling(al/dx)-1
    lx=NGrid+2*NAbsorbGrid
    allocate(x(NGrid+2*NAbsorbGrid))
        j=NAbsorbGrid+NGrid
        forall(i=1:NAbsorbGrid)
            x(i)=left-(NAbsorbGrid-i+1)*dx
            x(i+j)=right+i*dx
        end forall
        forall(i=1:NGrid)
            x(i+NAbsorbGrid)=left+(i-1)*dx
        end forall
    totallength=surface*lx
    allocate(hamiltonian_absorb(totallength,totallength))
    hamiltonian_absorb=-ci*HamiltonianDVR_absorb(x,dx,surface)
    if(steptype=='Auto'.and.lx<maxlx) then
        allocate(eigval(totallength))
        allocate(eigvec(totallength,totallength))
        call My_zgeev('N',hamiltonian_absorb,eigval,eigvec,totallength,success)
        dt=min(1d0/maxval(abs(eigval)),pi*4d0*mass/kmax**2,maxdt)
        deallocate(eigval)
        deallocate(eigvec)
        call dScientificNotation(dt,i)
        dt=floor(dt)*10d0**i
    end if
    lt=floor(totaltime/dt)+1
    write(*,*)'dt=',dt
    allocate(psy(lx,lt,surface))
    psy=0d0
    allocate(psyevolve(totallength,lt))
    psyevolve=0d0
    do j=1,surface
        do i=1,lx
            psy(i,1,j)=psy0(x(i),j)
        end do
    end do
    forall(j=1:surface)
        psyevolve((j-1)*lx+1:j*lx,1)=psy(:,1,j)
    end forall
    do i=2,lt
        call zRK4(psyevolve(:,i-1),psyevolve(:,i),evolve_absorb,dt,totallength)
        if(abs(dot_product(psyevolve(:,i),psyevolve(:,i))*dx-1d0)>minpopdev) then
            exit
        end if
    end do
    actualtime=(min(i,lt)-1)*dt
    forall (i=2:lt)
        forall(j=1:surface)
            psy(:,i,j)=psyevolve((j-1)*lx+1:j*lx,i)
        end forall
    end forall
    deallocate(hamiltonian_absorb)
    deallocate(psyevolve)
end subroutine solve

!numerically solve the evolution of different psy0 with difference p0, and save the transmission and reflection
subroutine scan_solve(Surface,Transmission,Reflection)
    logical::success
    integer::surface,i,j,k,totallength,nleft,nright
    real*8::al
    real*8,dimension(surface)::transmission,reflection,reigval
    complex*16,allocatable,dimension(:)::eigval,psyevolve,psyevolve_absorb
    complex*16,allocatable,dimension(:,:)::eigvec,psyleft,psyright,psyleftab,psyrightab
    call initialize()
    transmission=0d0
    reflection=0d0
    if(kminabs0) then
        stop 'parameter error: kmin<0 .and. kmax>0'
    end if
    if(steptype=='Auto') then
        dx=min(pim2/5d0/kmax,maxdx)
        call dScientificNotation(dx,i)
        dx=floor(dx)*10d0**i
    end if
    NGrid=floor((right-left)/dx)+1
    dx=(right-left)/(NGrid-1)
    write(*,*)'dx=',dx
    al=pim2/kmin
    NAbsorbGrid=ceiling(al/dx)-1
    lx=NGrid+2*NAbsorbGrid
    allocate(x(NGrid+2*NAbsorbGrid))
        j=NAbsorbGrid+NGrid
        forall(i=1:NAbsorbGrid)
            x(i)=left-(NAbsorbGrid-i+1)*dx
            x(i+j)=right+i*dx
        end forall
        forall(i=1:NGrid)
            x(i+NAbsorbGrid)=left+(i-1)*dx
        end forall
        nleft=NAbsorbGrid
        nright=j+1
    allocate(psyleft(NAbsorbGrid,surface))
    allocate(psyright(NAbsorbGrid,surface))
    allocate(psyleftab(NAbsorbGrid,surface))
    allocate(psyrightab(NAbsorbGrid,surface))
    totallength=surface*lx
    allocate(hamiltonian_absorb(totallength,totallength))
    hamiltonian=-ci*HamiltonianDVR(x,dx,surface)
    allocate(hamiltonian_absorb(totallength,totallength))
    hamiltonian_absorb=-ci*HamiltonianDVR_absorb(x,dx,surface)
    if(steptype=='Auto'.and.lx<maxlx) then
        allocate(eigval(lx))
        allocate(eigvec(lx,lx))
        call my_zgeev('N',hamiltonian,eigval,eigvec,totallength,success)
        dt=min(1d0/maxval(abs(eigval)),pi*4d0*mass/kmax**2,maxdt)
        call my_zgeev('N',hamiltonian_absorb,eigval,eigvec,totallength,success)
        dt=min(dt,1d0/maxval(abs(eigval)))
        deallocate(eigval)
        deallocate(eigvec)
        call dScientificNotation(dt,i)
        dt=floor(dt)*10d0**i
    end if
    lt=floor(totaltime/dt)+1
    dt=totaltime/(lt-1)
    write(*,*)'dt=',dt
    allocate(psyevolve(totallength))
    psyevolve=0d0
    allocate(psyevolve_absorb(totallength))
    psyevolve_absorb=0d0
    do j=1,surface
        do i=1,lx
            psyevolve_absorb(i)=psy0(x(i),j)
        end do
    end do
    do i=2,lt
        call zRK4(psyevolve_absorb,psyevolve,evolve,dt,totallength)
        call zRK4(psyevolve_absorb,psyevolve_absorb,evolve_absorb,dt,totallength)
        forall(j=1:surface)
            psyleft(:,j)=psyevolve((j-1)*lx+1:(j-1)*lx+nleft)
            psyright(:,j)=psyevolve((j-1)*lx+nright:j*lx)
            psyleftab(:,j)=psyevolve_absorb((j-1)*lx+1:(j-1)*lx+nleft)
            psyrightab(:,j)=psyevolve_absorb((j-1)*lx+nright:j*lx)
            reflection(j)=reflection(j)&
                +dx*(dot_product(psyleft(:,j),psyleft(:,j))&
                -dot_product(psyleftab(:,j),psyleftab(:,j)))
            transmission(j)=transmission(j)&
                +dx*(dot_product(psyright(:,j),psyright(:,j))&
                -dot_product(psyrightab(:,j),psyrightab(:,j)))
        end forall
        if(real(dot_product(psyevolve_absorb,psyevolve_absorb))*dx<minpop) then
            exit
        end if
    end do
    actualtime=(min(i,lt)-1)*dt
    deallocate(hamiltonian)
    deallocate(hamiltonian_absorb)
    deallocate(psyevolve)
    deallocate(psyevolve_absorb)
    deallocate(psyleft)
    deallocate(psyright)
    deallocate(psyleftab)
    deallocate(psyrightab)
end subroutine scan_solve

!get the matrix form of x, p, xx, xp, pp...
!note that DVR can't be considered as quantum Hilbert space, 
!the expectation of A(t) is dx*dot_product(psy(:,t),matmul(A,psy(:,t)))
subroutine QHDMatrix(nQHD,MatrixForm)
    integer::nQHD,i,j,k
    complex*16,dimension(lx,lx,nQHD)::MatrixForm
    MatrixForm=(0d0,0d0)
    do i=1,QHDOrder
        k=(i-1)*(i+2)/2+1
       forall(j=1:lx)
           MatrixForm(j,j,k)=x(j)**i
       end forall
    end do
    do i=1,lx
        do j=1,lx
            if(i==j) then
                MatrixForm(i,j,2)=0d0
                MatrixForm(i,j,5)=pipid3
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
    forall(i=1:lx)
        MatrixForm(:,i,4)=MatrixForm(:,i,2)*x(i)
    end forall
    MatrixForm(:,:,4)=MatrixForm(:,:,4)-ci/2d0*hbar
    if(QHDOrder>2) then
        forall(i=1:lx)
            MatrixForm(:,i,7)=MatrixForm(:,i,2)*x(i)**2
            MatrixForm(:,i,8)=MatrixForm(:,i,5)*x(i)
        end forall
        MatrixForm(:,:,7)=MatrixForm(:,:,7)-ci*hbar*MatrixForm(:,:,1)
        MatrixForm(:,:,8)=MatrixForm(:,:,8)-ci*hbar*MatrixForm(:,:,2)
    end if
end subroutine QHDMatrix

!phi(p)=F[psy(x)]
subroutine Transform2p(psy,lx,phi,lk)
    integer::lx,lk,k,i
    real*8,dimension(lk)::ones
    complex*16,dimension(lx)::psy
    complex*16,dimension(lk)::phi,temp
    ones=1d0
    do k=1,lk
        forall(i=NAbsorbGrid+1:NAbsorbGrid+NGrid)
            temp(i-NAbsorbGrid)=exp(-ci/hbar*p0scan(k)*x(i))*psy(i)
        end forall
        phi(k)=dot_product(ones,temp)*dx/pim2
    end do
end subroutine Transform2p

!p(x,p)=WignerTransform[psy(x)]
subroutine Transform2Wigner(psy,lx,lk,wigner)
    integer::lx,lk,k,i,j,width
    real*8,dimension(lk)::ones
    real*8,dimension(lx,lk)::wigner
    complex*16,dimension(lx)::psy
    complex*16,dimension(lk)::temp
    ones=1d0
    do k=1,lk
        do j=2,lk-1
            temp=(0d0,0d0)
            width=min(j,lk-j)-1
            forall(i=-width:width)
                temp(i+width+1)=exp(2d0*ci/hbar*p0scan(k)*i*dx)*psy(j+i+NAbsorbGrid)*psy(j-i+NAbsorbGrid)
            end forall
            wigner(j,k)=dot_product(ones,temp)
        end do
    end do
end subroutine Transform2Wigner

end module DVR