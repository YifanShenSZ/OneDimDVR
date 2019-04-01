module General
    implicit none    

!parameters
    !for Factorial
        logical::FactorialWarning=1

!global variables
    logical::kminabs0
    integer::NGrid,NAbsorbGrid,lx,lt,lp0
    real*8::d,dx,dt,kmin,kmax,kmins,actualtime
    real*8,allocatable,dimension(:)::x,t,p0scan
    real*8,allocatable,dimension(:,:,:,:)::wigner
    complex*16,allocatable,dimension(:,:,:)::psy,phi

!inputs
    character*32::jobtype,steptype
    integer::surface,QHDOrder
    real*8::mass,x0,totaltime,left,right,maxdx,maxdt,p0,p0left,p0right,dp0

!constants
    real*8::hbar=1d0,pi=3.141592653589793d0,enature=2.718281828459045d0,c=2.62206d0,&
            pim2=6.283185307179586d0,pipid3=3.2898681336964524d0,cs=6.87519864356d0
    complex*16::ci=(0d0,1d0)

contains
subroutine ReadInput()
    integer::i
    open(unit=99,file='input')
        read(99,*)
        read(99,*)jobtype
        read(99,*)
        read(99,*)steptype
        read(99,*)
        read(99,*)surface
        read(99,*)
        read(99,*)mass
        read(99,*)
        read(99,*)x0
        read(99,*)
        read(99,*)totaltime
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
        read(99,*)QHDOrder
        read(99,*)
        read(99,*)p0left
        read(99,*)
        read(99,*)p0right
        read(99,*)
        read(99,*)dp0
    close(99)
    write(*,*)'Job type: '//jobtype
    write(*,*)'Step type: '//steptype
end subroutine ReadInput

subroutine ReadTrajectory()
    integer::i,j,k
    real*8::re,im
    open(unit=99,file='ParametersUsed.DVR')
        read(99,*)NGrid
        read(99,*)NAbsorbGrid
        read(99,*)lx
        read(99,*)dx
        read(99,*)actualtime
        read(99,*)lt
        read(99,*)dt
        read(99,*)lp0
        read(99,*)surface
    close(99)
    allocate(x(lx))
        open(unit=99,file='x.DVR')
            do i=1,lx
                read(99,*)x(i)
            end do
        close(99)
    allocate(psy(lx,lt,surface))
        open(unit=99,file='Psy.DVR')
            do k=1,surface
                do i=1,lt
                    do j=1,lx
                        Read(99,*)re
                        Read(99,*)im
                        psy(j,i,k)=cmplx(re,im)
                    end do
                end do
            end do
        close(99)
end subroutine ReadTrajectory

!jobtype: 'N' eigen values only, 'V' eigen vectors as well
!N order general matrix A
!eigval harvests the eigen values, eigvec harvest the eigen vectors
!only for right eigen: matmul(A,eigvec(:,j))=eigval(j)*eigvec(:,j)
!success harvest whether the iteration succeeded
subroutine My_zgeev(jobtype,A,eigval,eigvec,N,success)
    logical::success
    character::jobtype
    integer::n,lwork,info
    real*8,dimension(2*n)::rwork
    complex*16,dimension(3*n)::work
    complex*16,dimension(n)::eigval
    complex*16,dimension(1,n)::vl
    complex*16,dimension(n,n)::a,temp,eigvec
    temp=a
    info=0
    call zgeev('N',jobtype,n,temp,n,eigval,vl,1,eigvec,n,work,3*n,rwork,info)
    success=info==0
end subroutine My_zgeev

!jobtype: 'N' eigen values only, 'V' eigen vectors as well
!N order hermitian matrix A
!eigval harvests the eigen values, eigvec harvest the eigen vectors
!success harvest whether the iteration succeeded
subroutine My_zheev(jobtype,A,eigval,eigvec,N,success)
    logical::success
    character::jobtype
    integer::N,info
    real*8,dimension(N)::eigval
    real*8,dimension(3*N-2)::rwork
    complex*16,dimension(3*N)::work
    complex*16,dimension(N,N)::A,eigvec
    info=0
    eigvec=A
    call zheev(jobtype,'U',n,eigvec,n,eigval,work,3*N,rwork,info)
    success=info==0
end subroutine My_zheev

!f has the form of: subroutine f(derivative,u,dim)
!Runge Kutta 4 order 
subroutine zRK4(old,new,f,dt,dim)
    external::f
    integer::dim
    complex*16,dimension(dim)::old,new,k1,k2,k3,k4
    real*8::dt,dtd2
    dtd2=dt/2d0
    call f(k1,old,dim)
    call f(k2,old+k1*dtd2,dim)
    call f(k3,old+k2*dtd2,dim)
    call f(k4,old+k3*dt,dim)
    new=old+dt/6d0*(k1+2d0*k2+2d0*k3+k4)
end subroutine zRK4

!8 significant digits for n>=41
function Factorial(n)
    integer::n,i
    real*8::Factorial,temp
    select case(n)
        case(0)
            Factorial=1d0
        case(1)
            Factorial=1d0
        case(2)
            Factorial=2d0
        case(3)
            Factorial=6d0
        case(4)
            Factorial=24d0
        case(5)
            Factorial=120d0
        case(6)
            Factorial=720d0
        case(7)
            Factorial=5040d0
        case(8)
            Factorial=40320d0
        case(9)
            Factorial=362880d0
        case(10)
            Factorial=3628800d0
        case(11)
            Factorial=39916800d0
        case(12)
            Factorial=479001600d0
        case(13)
            Factorial=6227020800d0
        case(14)
            Factorial=87178291200d0
        case(15)
            Factorial=1307674368d3
        case(16)
            Factorial=20922789888d3
        case(17)
            Factorial=355687428096d3
        case(18)
            Factorial=6402373705728d3
        case(19)
            Factorial=121645100408832d3
        case(20)
            Factorial=243290200817664d4
        case default
            if(FactorialWarning) then
                write(*,*)'Not accurate factorial'
                FactorialWarning=0
            end if
            temp=(9d0*n)**pi
            Factorial=sqrt(2*pi*n)*(n/enature)**n*Exp(1d0/12d0/n-Log(9d0*n)/(temp-1d0/temp))
    end select
end function Factorial
function Factorial2(n)
    integer::n,i
    real*8::Factorial2
    select case(n)
        case(-1)
            Factorial2=1d0
        case(0)
            Factorial2=1d0
        case(1)
            Factorial2=1d0
        case(2)
            Factorial2=2d0
        case(3)
            Factorial2=3d0
        case(4)
            Factorial2=8d0
        case(5)
            Factorial2=15d0
        case(6)
            Factorial2=48d0
        case(7)
            Factorial2=105d0
        case(8)
            Factorial2=384d0
        case(9)
            Factorial2=945d0
        case(10)
            Factorial2=3840d0
        case(11)
            Factorial2=10395d0
        case(12)
            Factorial2=46080d0
        case(13)
            Factorial2=135135d0
        case(14)
            Factorial2=645120d0
        case(15)
            Factorial2=2027025d0
        case(16)
            Factorial2=10321920d0
        case(17)
            Factorial2=34459425d0
        case(18)
            Factorial2=185794560d0
        case(19)
            Factorial2=654729075d0
        case(20)
            Factorial2=3715891200d0
        case(21)
            Factorial2=13749310575d0
        case(22)
            Factorial2=81749606400d0
        case(23)
            Factorial2=316234143225d0
        case(24)
            Factorial2=1961990553600d0
        case(25)
            Factorial2=7905853580625d0
        case(26)
            Factorial2=51011754393600d0
        case(27)
            Factorial2=213458046676875d0
        case(28)
            Factorial2=1428329123020800d0
        case(29)
            Factorial2=6190283353629375d0
        case(30)
            Factorial2=42849873690624000d0
        case(31)
            Factorial2=191898783962510625d0
        case(32)
            Factorial2=1371195958099968000d0
        case(33)
            Factorial2=6332659870762850625d0
        case(34)
            Factorial2=46620662575398912000d0
        case(35)
            Factorial2=221643095476699771875d0
        case(36)
            Factorial2=1678343852714360832000d0
        case(37)
            Factorial2=8200794532637891559375d0
        case(38)
            Factorial2=63777066403145711616000d0
        case(39)
            Factorial2=319830986772877770815625d0
        case(40)
            Factorial2=2551082656125828464640000d0
        case default
            if(mod(n,2)) then
                Factorial2=Factorial(n+1)/(2**((n+1)/2)*Factorial((n+1)/2))
            else
                Factorial2=2**(n/2)*Factorial(n/2)
            end if
    end select
end function Factorial2
function Arrangement(m,n)
    integer::m,n,i
    real*8::Arrangement
    if(m<2.or.n==0) then
        Arrangement=1d0
    else if(n==1) then
        Arrangement=m
    else if(n==m.or.n==m-1) then
        Arrangement=Factorial(m)
    else
        select case(m)
            case(4)
                Arrangement=12d0
            case(5)
                select case(n)
                    case(2)
                        Arrangement=20d0
                    case(3)
                        Arrangement=60d0
                end select
            case(6)
                select case(n)
                    case(2)
                        Arrangement=30d0
                    case(3)
                        Arrangement=120d0
                    case(4)
                        Arrangement=360d0
                end select
            case(7)
                select case(n)
                    case(2)
                        Arrangement=42d0
                    case(3)
                        Arrangement=210d0
                    case(4)
                        Arrangement=840d0
                    case(5)
                        Arrangement=2520d0
                end select
            case(8)
                select case(n)
                    case(2)
                        Arrangement=56d0
                    case(3)
                        Arrangement=336d0
                    case(4)
                        Arrangement=1680d0
                    case(5)
                        Arrangement=6720d0
                    case(6)
                        Arrangement=20160d0
                end select
            case(9)
                select case(n)
                    case(2)
                        Arrangement=72d0
                    case(3)
                        Arrangement=504d0
                    case(4)
                        Arrangement=3024d0
                    case(5)
                        Arrangement=15120d0
                    case(6)
                        Arrangement=60480d0
                    case(7)
                        Arrangement=181440d0
                end select
            case(10)
                select case(n)
                    case(2)
                        Arrangement=90d0
                    case(3)
                        Arrangement=720d0
                    case(4)
                        Arrangement=5040d0
                    case(5)
                        Arrangement=30240d0
                    case(6)
                        Arrangement=151200d0
                    case(7)
                        Arrangement=604800d0
                    case(8)
                        Arrangement=1814400d0
                end select
            case default
                if(m>20.and.n<10) then
                    Arrangement=m*(m-1)
                    do i=m-2,m-n+1,-1
                        Arrangement=Arrangement*i
                    end do
                else
                    Arrangement=Factorial(m)/Factorial(m-n)
                end if
        end select
    end if
end function Arrangement
function Combination(m,n)
    integer::m,n
    real*8::Combination
    if(m<2.or.n==0.or.n==m) then
        Combination=1d0
    else if(n==1.or.n==(m-1)) then
        Combination=m
    else 
        select case(m)
            case(4)
                Combination=6d0
            case(5)
                Combination=10d0
            case default
                if(n<m/2d0) then
                    n=m-n
                end if
                select case(m)
                    case(6)
                        select case(n)
                            case(4)
                                Combination=15d0
                            case(3)
                                Combination=20d0
                        end select
                    case(7)
                        select case(n)
                            case(5)
                                Combination=21d0
                            case(4)
                                Combination=35d0
                        end select
                    case(8)
                        select case(n)
                            case(6)
                                Combination=28d0
                            case(5)
                                Combination=56d0
                            case(4)
                                Combination=70d0
                        end select
                    case(9)
                        select case(n)
                            case(7)
                                Combination=36d0
                            case(6)
                                Combination=84d0
                            case(5)
                                Combination=126d0
                        end select
                    case(10)
                        select case(n)
                            case(8)
                                Combination=45d0
                            case(7)
                                Combination=120d0
                            case(6)
                                Combination=210d0
                            case(5)
                                Combination=252d0
                        end select
                    case(11)
                        select case(n)
                            case(9)
                                Combination=55d0
                            case(8)
                                Combination=165d0
                            case(7)
                                Combination=330d0
                            case(6)
                                Combination=462d0
                        end select
                    case default
                        Combination=Arrangement(m,n)/Factorial(n)
                end select
        end select
    end if
end function Combination

!x_input=x_output*1di
subroutine dScientificNotation(x,i)
    real*8::x
    integer::i
    i=0
    do
        if(x<1d0) then
            x=x*10d0
            i=i-1
        else
            exit
        end if
    end do
    do
        if(x>=10d0) then
            x=x/10d0
            i=i+1
        else
            exit
        end if
    end do
end subroutine dScientificNotation

!show date hour minute second
subroutine ShowTime()
    integer::value(1:8)
    call date_and_time(values=value)
    write(*,*)value(3),'d',value(5),':',value(6),':',value(7)
end subroutine ShowTime

end module General