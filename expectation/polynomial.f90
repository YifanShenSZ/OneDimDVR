module polynomial
    use StringUtility
    implicit none

!Polynomial definition
type PolDef
    integer, allocatable, dimension(:)::q, p
    contains
    procedure::construct => PolDef_construct
    procedure::evaluate => PolDef_evaluate
end type PolDef

contains
subroutine PolDef_construct(this, order, input_line)
    class(PolDef), intent(inout)::this
    integer, intent(in)::order
    character(*), dimension(order), intent(in)::input_line
    integer::qcount, pcount, i, j
    !Count the number of coordinate and momentum monomials
    qcount = 0
    pcount = 0
    do i = 1, order
        read(input_line(i),*)j
        if (j > 0) then
            qcount = qcount + 1
        else
            pcount = pcount + 1
        end if
    end do
    !Read the definition
    allocate(this%q(qcount))
    allocate(this%p(pcount))
    qcount = 0
    pcount = 0
    do i = 1, order
        read(input_line(i),*)j
        if (j > 0) then
            qcount = qcount + 1
            this%q(qcount) = j
        else
            pcount = pcount + 1
            this%p(pcount) = j
        end if
    end do
end subroutine PolDef_construct

!This is specialized for 1-dimensional case
real*8 function PolDef_evaluate(this, q, p)
    class(PolDef), intent(in)::this
    real*8, intent(in)::q, p
    integer::i
    PolDef_evaluate = 1d0
    do i = 1, size(this%q)
        PolDef_evaluate = PolDef_evaluate * q
    end do
    do i = 1, size(this%p)
        PolDef_evaluate = PolDef_evaluate * p
    end do
end function PolDef_evaluate

end module polynomial