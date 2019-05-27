!Library interface, data storage, basic routine
!Modify this module for different initial wave function (psy0) and potential:
!    REMEMBER to modify InitializeDVRParameter accordingly,
!    different models have different assumptions on initial condition
module Basic
    use General
    use Mathematics
    use LinearAlgebra
    implicit none

!Parameter
    integer,parameter::SMDOrder=2!Higher order SMD operators in DVR have not been derived

!Programwide accessed input variable
    !General control
    integer::NState
    real*8::mass,x0,left,right,maxdx,maxdt
    !NewTrajectory only
    real*8::p0,TotalTime,OutputInterval
    logical::ScatteringProblem
    !TR-p0 only
    real*8::p0left,p0right,dp0

!Global variable
    integer::NGrid,NAbsorbGrid,lx,lt,lp0
    real*8::d,dx,dt,kmin,kmax,kmins
    real*8,allocatable,dimension(:)::x,t,p0scan
    real*8,allocatable,dimension(:,:,:,:)::wigner
    complex*16,allocatable,dimension(:,:,:)::psy,phi

!Constant
    real*8,parameter::hbar=1d0,c=2.62206d0,cs=6.87519864356d0

contains
!Initial wave function
function psy0(x,i)!Gaussian wave packet on 1st state, taking x0 and p0 as arguments
    integer::i
    real*8::x
    complex*16::psy0
    if(i==1) then
        psy0=(pi*2d0*d)**(-0.25d0)*exp(-(x-x0)**2/4d0/d+ci*p0*x)
    end if
end function psy0

!Under such psy0, d (the initial variance of x) is the only thing user has to worry about:
!    For scattering problem, here is coded to use sigma_p=p0/20d0, equivalent to d=100d0/p0**2
!    For bounded problem, user has to specify d
subroutine InitializeDVRParameter()
    real*8::temp
    if(ScatteringProblem) then
        dx=maxdx
        dt=maxdt
        !For gaussian wave packet psy0 only
            kmin=p0-3d0*p0/20d0
            kmax=p0+3d0*p0/20d0
            d=hbar**2*100d0/p0**2!Automatically use sigma_p=p0/20d0, equivalent to d=100d0/p0**2
        if(kmin<0.and.kmax>0) stop 'Program abort: wavepacket contains 0 momentum component, please check your input'
        temp=min(abs(kmin),abs(kmax))
        kmax=max(abs(kmin),abs(kmax))
        kmin=temp/3d0
        kmins=kmin*kmin
    else
        dx=maxdx
        dt=maxdt
        kmax=0d0
        !User have to specify the initial variance of x
        d=0.1776604496914905d0**2!6-th order double well,x0=-0.463950899364933
    end if
end subroutine InitializeDVRParameter

real*8 function potential(x,i,j)!Some scattering models take mass as argument
    real*8,intent(in)::x
    integer,intent(in)::i,j
    !Scattering model
        !Gaussian well
            !real*8,parameter::p_peak=10d0,miu=0d0,sigma=3d0
            !potential=-p_peak**2*0.5d0/mass*exp(-(x-miu)**2/2d0/sigma**2)
        !Gaussian barrier
            !real*8,parameter::p_peak=10d0,miu=0d0,sigma=3d0
            !potential=p_peak**2*0.5d0/mass*exp(-(x-miu)**2/2d0/sigma**2)
        !Erfc down step
            !real*8,parameter::p_peak=10d0,halfwidth=5d0
            !potential=p_peak**2*0.25d0/mass*Erfc(x/halfwidth)
        !Erfc up step
            !real*8,parameter::p_peak=10d0,halfwidth=5d0
            !potential=-p_peak**2*0.25d0/mass*Erfc(x/halfwidth)
        !Simple avoided crossing
            if(i==1.and.j==1) then
                if(x>0d0) then
                    potential=0.01d0*(1d0-Exp(-1.6d0*x))
                else
                    potential=0.01d0*(Exp(1.6d0*x)-1d0)
                end if
            else if(i==2.and.j==2) then
                if(x>0d0) then
                    potential=0.01d0*(Exp(-1.6d0*x)-1d0)
                else
                    potential=0.01d0*(1d0-Exp(1.6d0*x))
                end if
            else
                potential=0.005d0*exp(-x**2)
            end if
        !Dual avoided crossing
            !if(i==1.and.j==1) then
            !    potential=0d0
            !else if(i==2.and.j==2) then
            !    potential=0.05d0-0.1d0*exp(-0.28d0*x**2)
            !else
            !    potential=0.015d0*exp(-0.06d0*x**2)
            !end if
        !Extended coupling with reflection
            !if(i==1.and.j==1) then
            !    potential=6d-4
            !else if(i==2.and.j==2) then
            !    potential=-6d-4
            !else
            !    if(x>0d0) then
            !        potential=0.1d0*(2d0-exp(-0.9d0*x))
            !    else
            !        potential=0.1d0*exp(0.9d0*x)
            !    end if
            !end if
    !Bounded model
        !4th order double well
            !real*8::x2
            !x2=x*x
            !potential=x2*(2.4d-5*x2-3d-4)
        !6th order double well
            !real*8::x2,x4
            !x2=x*x
            !x4=x2*x2
            !potential=3d-2*x2*(x4+2d0*x2-1d0)
end function potential

end module Basic