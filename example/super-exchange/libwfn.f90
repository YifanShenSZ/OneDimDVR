module libwfn
    implicit none

contains
subroutine initialize_libwfn()
end subroutine initialize_libwfn

complex*16 function init_wfn(q,i)
    real*8,intent(in)::q
    integer,intent(in)::i
    real*8, parameter::qsigma = 1d0, qmean = -10d0, pmean = 0d0
    complex*16, parameter::ci = (0d0,1d0)
    real*8::q_dimless
    q_dimless = (q - qmean) / qsigma
    if (i == 1) then
        !init_wfn = (pi*2d0*d) **(-0.25d0)*exp(-(q-q0)**2/4d0/d+ci*p0*q)
        init_wfn = 0.6316187777460647d0 / Sqrt(qsigma) & !Sqrt{1/[Sqrt(2pi)*qsigma]}
                 * exp(- q_dimless*q_dimless / 4d0 + ci * pmean * q)
    end if
end function init_wfn

end module libwfn