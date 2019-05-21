module DVR
    use Basic
    implicit none

!Parameter
    logical::AutoStep=.true.!Whether determine dx and dt automatically
    integer::maxtotallength=2147483647!To prevent too slow generation of auto dt
    real*8::minpopdev=1d-6,minpop=1d-6!Stop evolution when population no longer equals to 1

!DVR module only variable
    complex*16,allocatable,dimension(:,:)::Hamiltonian,HamiltonianAbsorb

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
!HamiltonianDVR store as a whole NState*size(x) order matrix
!can be considered as consisting of NState*NState blocks, ij block is NState i x NState j 
function HamiltonianDVR(x,dx,NState,lx)
    real*8,dimension(lx),intent(in)::x
    real*8,intent(in)::dx
    integer,intent(in)::NState,lx
    integer::i,j,k,ii,jj,kk,kkk,totallength
    real*8,dimension(lx,lx)::tran
    complex*16,allocatable,dimension(:,:)::HamiltonianDVR
    totallength=NState*lx
    allocate(HamiltonianDVR(totallength,totallength))
    HamiltonianDVR=(0d0,0d0)
    do i=1,lx
        do j=1,lx
            if(i==j) then
                tran(i,j)=pisqd3
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
    do i=1,NState
        k=j+lx-1
        HamiltonianDVR(j:k,j:k)=tran
        do ii=1,NState
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
function HamiltonianDVRAbsorb(x,dx,NState,lx)
    real*8,dimension(lx),intent(in)::x
    real*8,intent(in)::dx
    integer,intent(in)::NState,lx
    integer::i,j,k,ii,jj,kk,kkk,totallength
    real*8,dimension(lx,lx)::tran
    complex*16,allocatable,dimension(:,:)::HamiltonianDVRAbsorb
    totallength=NState*lx
    allocate(HamiltonianDVRAbsorb(totallength,totallength))
    HamiltonianDVRAbsorb=(0d0,0d0)
    do i=1,lx
        do j=1,lx
            if(i==j) then
                tran(i,j)=pisqd3
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
    do i=1,NState
        k=j+lx-1
        HamiltonianDVRAbsorb(j:k,j:k)=tran
        do ii=1,NState
            kk=(i-1)*lx
            kkk=(ii-1)*lx
            do jj=1,lx
                kk=kk+1
                kkk=kkk+1
                HamiltonianDVRAbsorb(kk,kkk)=HamiltonianDVRAbsorb(kk,kkk)+AbsorbedPotential(x(jj),i,ii)
            end do
        end do
        j=k+1
    end do
end function HamiltonianDVRAbsorb

subroutine Evolve(d,u,N)
    complex*16,dimension(N),intent(inout)::d
    complex*16,dimension(N),intent(in)::u
    integer,intent(in)::N
    d=matmul(Hamiltonian,u)
end subroutine Evolve
subroutine EvolveAbsorb(d,u,N)
    complex*16,dimension(N),intent(inout)::d
    complex*16,dimension(N),intent(in)::u
    integer,intent(in)::N
    d=matmul(HamiltonianAbsorb,u)
end subroutine EvolveAbsorb

!Numerically solve the evolution of psy0, and output the trajectory
subroutine Solve()
    integer::i,j,index,totallength,OutputFreq
    real*8::al
    complex*16,allocatable,dimension(:)::eigval,psyevolveold,psyevolvenew
    complex*16,allocatable,dimension(:,:)::psyevolve,eigvec
    !discretize space
        if(AutoStep) then
            dx=min(pim2/5d0/kmax,maxdx)
            call dScientificNotation(dx,i)
            dx=floor(dx)*10d0**i
        end if
        write(*,*)'dx =',dx
        NGrid=floor((right-left)/dx)+1
        NAbsorbGrid=0
        lx=NGrid
        allocate(x(NGrid))
        forall(i=1:NGrid)
            x(i)=left+(i-1)*dx
        end forall
        totallength=NState*lx
        allocate(Hamiltonian(totallength,totallength))
        Hamiltonian=-ci*HamiltonianDVR(x,dx,NState,lx)
    !discretize time
        if(totallength>maxtotallength) then
            write(*,*)'To prevent too slow generation of auto dt, auto dt disabled'
            write(*,*)'The dimension of Hamiltonian is',totallength
        elseif(AutoStep) then
            allocate(eigval(totallength))
            allocate(eigvec(totallength,totallength))
            call My_zgeev('N',Hamiltonian,eigval,eigvec,totallength)
            dt=min(1d0/maxval(abs(eigval)),pim2*mass/hbar/kmax**2,maxdt)
            deallocate(eigval)
            deallocate(eigvec)
        end if
        OutputFreq=Ceiling(OutputInterval/dt)
        dt=OutputInterval/OutputFreq
        write(*,*)'dt =',dt
    !prepare
        lt=floor(totaltime/OutputInterval)+1
        allocate(psy(lx,NState,lt))
        allocate(psyevolveold(totallength))
        allocate(psyevolvenew(totallength))
        index=1
        do j=1,NState
            do i=1,lx
                psy(i,j,1)=psy0(x(i),j)
                psyevolveold(index)=psy(i,j,1)
                index=index+1
            end do
        end do
    !solve
        write(*,*)'Total snap shots =',(lt-1)*OutputFreq
        write(*,*)'Evolving...'
        do i=1,(lt-1)*OutputFreq
            call zRK4(psyevolveold,psyevolvenew,Evolve,dt,totallength)
            !write trajectory
                if(mod(i,OutputFreq)==0) then
                    forall(j=1:NState)
                        psy(:,j,1+i/OutputFreq)=psyevolvenew((j-1)*lx+1:j*lx)
                    end forall
                end if
            !get ready for next loop
                psyevolveold=psyevolvenew
        end do
        ActualTime=(lt-1)*OutputInterval
    !clean up
        deallocate(Hamiltonian)
        deallocate(psyevolveold)
        deallocate(psyevolvenew)
end subroutine Solve
subroutine SolveAbsorb()
    integer::i,j,index,totallength,OutputFreq
    real*8::al
    complex*16,allocatable,dimension(:)::eigval,psyevolveold,psyevolvenew
    complex*16,allocatable,dimension(:,:)::eigvec
    !discretize space
        if(AutoStep) then
            dx=min(pim2/5d0/kmax,maxdx)
            call dScientificNotation(dx,i)
            dx=floor(dx)*10d0**i
        end if
        write(*,*)'dx =',dx
        NGrid=floor((right-left)/dx)+1
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
        totallength=NState*lx
        allocate(HamiltonianAbsorb(totallength,totallength))
        HamiltonianAbsorb=HamiltonianDVRAbsorb(x,dx,NState,lx)
        HamiltonianAbsorb=-ci*HamiltonianAbsorb
    !discretize time
        if(totallength>maxtotallength) then
            write(*,*)'to prevent stackoverflow, auto dt disabled'
            write(*,*)'The dimension of Hamiltonian is',totallength
        elseif(AutoStep) then
            allocate(eigval(totallength))
            allocate(eigvec(totallength,totallength))
            call My_zgeev('N',HamiltonianAbsorb,eigval,eigvec,totallength)
            dt=min(1d0/maxval(abs(eigval)),pim2*mass/hbar/kmax**2,maxdt)
            deallocate(eigval)
            deallocate(eigvec)
        end if
        OutputFreq=Ceiling(OutputInterval/dt)
        dt=OutputInterval/OutputFreq
        write(*,*)'dt =',dt
    !prepare
        lt=floor(totaltime/OutputInterval)+1
        allocate(psy(lx,NState,lt))
        allocate(psyevolveold(totallength))
        allocate(psyevolvenew(totallength))
        index=1
        do j=1,NState
            do i=1,lx
                psy(i,j,1)=psy0(x(i),j)
                psyevolveold(index)=psy(i,j,1)
                index=index+1
            end do
        end do
    !solve
        write(*,*)'Total snap shots =',(lt-1)*OutputFreq
        write(*,*)'Evolving...'
        do i=1,(lt-1)*OutputFreq
            call zRK4(psyevolveold,psyevolvenew,EvolveAbsorb,dt,totallength)
            !write trajectory
                if(mod(i,OutputFreq)==0) then
                    forall(j=1:NState)
                        psy(:,j,1+i/OutputFreq)=psyevolvenew((j-1)*lx+1:j*lx)
                    end forall
                end if
            !check whether the population has been absorbed too much
                if(abs(dot_product(psyevolvenew,psyevolvenew)*dx-1d0)>minpopdev) then
                    lt=1+i/OutputFreq
                    exit
                end if
            !get ready for next loop
                psyevolveold=psyevolvenew
        end do
        ActualTime=(lt-1)*OutputInterval
    !clean up
        deallocate(HamiltonianAbsorb)
        deallocate(psyevolveold)
        deallocate(psyevolvenew)
end subroutine SolveAbsorb

!Numerically solve the evolution of different psy0 with difference p0, and save the transmission and reflection
subroutine scan_solve(Transmission,Reflection)
    real*8,dimension(NState),intent(inout)::transmission,reflection
    integer::i,j,index,totallength,nleft,nright
    real*8::al
    complex*16,allocatable,dimension(:)::eigval,psyevolve,psyevolve_absorb
    complex*16,allocatable,dimension(:,:)::eigvec
    transmission=0d0
    reflection=0d0
    !discretize space
        if(AutoStep) then
            dx=min(pim2/5d0/kmax,maxdx)
            call dScientificNotation(dx,i)
            dx=floor(dx)*10d0**i
        end if
        write(*,*)'dx =',dx
        NGrid=floor((right-left)/dx)+1
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
        totallength=NState*lx
        allocate(Hamiltonian(totallength,totallength))
        Hamiltonian=-ci*HamiltonianDVR(x,dx,NState,lx)
        allocate(HamiltonianAbsorb(totallength,totallength))
        HamiltonianAbsorb=-ci*HamiltonianDVRAbsorb(x,dx,NState,lx)
    !discretize time
        if(totallength>maxtotallength) then
            write(*,*)'to prevent stackoverflow, auto dt disabled'
            write(*,*)'The dimension of Hamiltonian is',totallength
        elseif(AutoStep) then
            allocate(eigval(totallength))
            allocate(eigvec(totallength,totallength))
            call My_zgeev('N',Hamiltonian,eigval,eigvec,totallength)
            dt=min(1d0/maxval(abs(eigval)),maxdt)
            call My_zgeev('N',HamiltonianAbsorb,eigval,eigvec,totallength)
            dt=min(dt,1d0/maxval(abs(eigval)))
            deallocate(eigval)
            deallocate(eigvec)
            call dScientificNotation(dt,i)
            dt=floor(dt)*10d0**i
        end if
        write(*,*)'dt =',dt
        lt=floor(totaltime/dt)+1
    !prepare
        allocate(psyevolve(totallength))
        allocate(psyevolve_absorb(totallength))
        index=1
        do j=1,NState
            do i=1,lx
                psyevolve_absorb(index)=psy0(x(i),j)
                index=index+1
            end do
        end do
    !solve
        do i=2,lt
            call zRK4(psyevolve_absorb,psyevolve,Evolve,dt,totallength)
            call zRK4(psyevolve_absorb,psyevolve_absorb,EvolveAbsorb,dt,totallength)
            forall(j=1:NState)
                reflection(j)=reflection(j)&
                    +dx*(dot_product(psyevolve((j-1)*lx+1:(j-1)*lx+nleft),psyevolve((j-1)*lx+1:(j-1)*lx+nleft))&
                    -dot_product(psyevolve_absorb((j-1)*lx+1:(j-1)*lx+nleft),psyevolve_absorb((j-1)*lx+1:(j-1)*lx+nleft)))
                transmission(j)=transmission(j)&
                    +dx*(dot_product(psyevolve((j-1)*lx+nright:j*lx),psyevolve((j-1)*lx+nright:j*lx))&
                    -dot_product(psyevolve_absorb((j-1)*lx+nright:j*lx),psyevolve_absorb((j-1)*lx+nright:j*lx)))
            end forall
            if(real(dot_product(psyevolve_absorb,psyevolve_absorb))*dx<minpop) then
                lt=i
                exit
            end if
        end do
        ActualTime=(lt-1)*dt
    !clean up
        deallocate(Hamiltonian)
        deallocate(HamiltonianAbsorb)
        deallocate(psyevolve)
        deallocate(psyevolve_absorb)
end subroutine scan_solve

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

!phi(p)=F[psy(x)]
subroutine Transform2p(psy,phi,NGrid)
    complex*16,dimension(NGrid),intent(in)::psy
    complex*16,dimension(NGrid),intent(inout)::phi
    integer,intent(in)::NGrid
    integer::k,i
    real*8,dimension(NGrid)::ones
    complex*16,dimension(NGrid)::temp
    ones=1d0
    do k=1,NGrid
        forall(i=1:NGrid)
            temp(i)=exp(-ci/hbar*p0scan(k)*x(i))*psy(i)
        end forall
        phi(k)=dot_product(ones,temp)
    end do
    phi=phi*dx/Sqrt(pim2*hbar)
end subroutine Transform2p

!p(x,p)=WignerTransform[psy(x)]
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
                temp(i+width+1)=exp(2d0*ci/hbar*p0scan(k)*i*dx)*psy(j+i)*psy(j-i)
            end forall
            wigner(j,k)=dot_product(ones,temp)
        end do
    end do
    wigner=wigner*dx/pi/hbar
end subroutine Transform2Wigner

end module DVR