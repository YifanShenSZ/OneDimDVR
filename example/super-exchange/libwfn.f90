!We set coordinate standard deviation = 10 / momentum expectation
module libwfn
    implicit none

real*8::qmean, pmean, qsigma

contains
subroutine initialize_libwfn()
    open(unit=99, file="initial.in", status="old")
        read(99,*); read(99,*)qmean
        read(99,*); read(99,*)pmean
    close(99)
    qsigma = 10d0 / pmean
end subroutine initialize_libwfn

subroutine init_wfn(q, psy)
    real*8, intent(in)::q
    complex*16, dimension(3), intent(out)::psy
    real*8::q_dimless
    q_dimless = (q - qmean) / qsigma
    psy(1) = 0.6316187777460647d0 / Sqrt(qsigma) & !Sqrt{1/[Sqrt(2pi)*qsigma]}
           * exp(-q_dimless * q_dimless / 4d0 + (0d0,1d0) * pmean * q)
    psy(2:3) = (0d0,0d0)
end subroutine init_wfn

end module libwfn