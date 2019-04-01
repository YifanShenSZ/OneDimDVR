!modify this module for different psy0 and potential
module V_Psy0
    use General
    implicit none

contains
!gaussian potential
    function potential(x,i,j)
        integer::i,j
        real*8::potential,x,&
                p_peak=8d0,miu=0d0,sigma=3d0
        potential=-p_peak**2*0.5d0/mass*exp(-(x-miu)**2/2d0/sigma**2)
    end function potential
!simple avoided crossing
    !function potential(x,i,j)
    !    integer::i,j
    !    real*8::potential,x,&
    !            p_peak=10d0
    !    if(i==1.and.j==1) then
    !        if(x>0d0) then
    !            potential=0.01d0*(1d0-Exp(-1.6d0*x))
    !        else
    !            potential=0.01d0*(Exp(1.6d0*x)-1d0)
    !        end if
    !    else if(i==2.and.j==2) then
    !        if(x>0d0) then
    !            potential=0.01d0*(Exp(-1.6d0*x)-1d0)
    !        else
    !            potential=0.01d0*(1d0-Exp(1.6d0*x))
    !        end if
    !    else
    !        potential=0.005d0*exp(-x**2)
    !    end if
    !end function potential
!dual avoided crossing
    !function potential(x,i,j)
    !    integer::i,j
    !    real*8::potential,x,&
    !            p_peak=10d0
    !    if(i==1.and.j==1) then
    !        potential=0d0
    !    else if(i==2.and.j==2) then
    !        potential=0.05d0-0.1d0*exp(-0.28d0*x**2)
    !    else
    !        potential=0.015d0*exp(-0.06d0*x**2)
    !    end if
    !end function potential
!extended coupling with reflection
    !function potential(x,i,j)
    !    integer::i,j
    !    real*8::potential,x,&
    !            p_peak=10d0
    !    if(i==1.and.j==1) then
    !        potential=6d-4
    !    else if(i==2.and.j==2) then
    !        potential=-6d-4
    !    else
    !        if(x>0d0) then
    !            potential=0.1d0*(2d0-exp(-0.9d0*x))
    !        else
    !            potential=0.1d0*exp(0.9d0*x)
    !        end if
    !    end if
    !end function potential

!initial wave function
!default is gaussian wave packet
function psy0(x,i)
    integer::i
    real*8::x
    complex*16::psy0
    if(i==1) then
        psy0=(pi*2d0*d)**(-0.25d0)*exp(-(x-x0)**2/4d0/d+ci*p0*x)
    end if
end function psy0
!call when start solving
!suits only gaussian wigner distribution .and. 0 potential out [left,right]
subroutine initialize()
    real*8::temp
    dx=maxdx
    dt=maxdt
    !sigma_p=p0/20d0, equal to d=100d0/p0**2
    d=hbar**2*100d0/p0**2
    !for gaussian wave packet psy0 only
    kmin=p0-3d0*p0/20d0
    kmax=p0+3d0*p0/20d0
    kminabs0=kmin<0.and.kmax>0
    temp=min(abs(kmin),abs(kmax))
    kmax=max(abs(kmin),abs(kmax))
    kmin=temp/3d0
    kmins=kmin**2
end subroutine initialize

end module V_Psy0