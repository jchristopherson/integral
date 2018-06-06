! integral_ode_auto.f90

submodule (integral_core) integral_ode_auto
contains
! ------------------------------------------------------------------------------
    module function oa_step(this, fcn, x, y, xout, err) result(brk)
        ! Arguments
        class(ode_auto), intent(inout) :: this
        class(ode_helper), intent(in) :: fcn
        real(real64), intent(inout) :: x
        real(real64), intent(inout), dimension(:) :: y
        real(real64), intent(in) :: xout
        class(errors), intent(inout), optional, target :: err
        logical :: brk

        ! Local Variables
    end function
end submodule
