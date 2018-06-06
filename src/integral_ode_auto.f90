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
        integer(int32) :: neq, itask, iopt, lrw, liw, jt, ncnst, lrn, lrs
        logical :: useConstraints
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        procedure(ode_fcn), pointer :: fptr

        ! Initialization
        neq = fcn%get_equation_count()
        ncnst = fcn%get_constraint_count()
        useConstraints = fcn%get_constraints_defined()
        fptr => fcn%get_ode_fcn()
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        if (this%get_allow_overshoot()) then
            itask = 1
            if (this%get_provide_all_output()) itask = 2
        else
            itask = 4
            if (this%get_provide_all_output()) itask = 5
        end if
        if (fcn%get_jacobian_defined()) then
            jt = 1
        else
            jt = 2
        end if
        iopt = 0

        ! Ensure workspace arrays are sized appropriately
        lrn = 20 + 16 * neq
        lrs = 22 + 9 * neq + neq**2
        lrw = max(lrn, lrs)
        liw = 20 + neq

        ! Determine if we are to call DLSODA or DLSODAR
        if (useConstraints .and. ncnst > 0) then
        else
        end if

        ! Check output condition

    contains
        subroutine odepackfcn(n, x, y, dydx)
            integer(int32), intent(in) :: n
            real(real64), intent(in) :: x, y(n)
            real(real64), intent(out) :: dydx(n)
            call fptr(x, y, dydx)
        end subroutine
    end function
end submodule
