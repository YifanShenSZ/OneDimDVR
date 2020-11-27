module libwfn
    implicit none

contains
subroutine initialize_libwfn()
end subroutine initialize_libwfn

subroutine init_wfn(q, psy)
    real*8, intent(in)::q
    complex*16, dimension(1), intent(out)::psy
    real*8, parameter::qsigma = 0.17766d0, qmean = -0.463951d0, pmean = 3d0
    complex*16, parameter::ci = (0d0,1d0)
    real*8::q_dimless
    q_dimless = (q - qmean) / qsigma
    psy(1) = 0.6316187777460647d0 / Sqrt(qsigma) & !Sqrt{1/[Sqrt(2pi)*qsigma]}
           * exp(- q_dimless*q_dimless / 4d0 + ci * pmean * q)
end subroutine init_wfn

end module libwfn