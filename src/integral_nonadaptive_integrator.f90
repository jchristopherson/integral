! integral_nonadaptive_integrator.f90

submodule (integral_core) integral_nonadaptive_integrator
contains
! ------------------------------------------------------------------------------
    module function ni_integrate(this, fcn, a, b, info, err) result(rst)
        ! Arguments
        class(nonadaptive_integrator), intent(inout) :: this
        procedure(integrand), intent(in), pointer :: fcn
        real(real64), intent(in) :: a, b
        type(integration_behavior), intent(out), optional :: info
        class(errors), intent(inout), optional, target :: err
        real(real64) :: rst

        ! Local Variables
        real(real64) :: abserr, resabs, resasc
        integer(int32) :: n

        ! Determine the correct rule, and perform the integration
        select case (this%get_rule())
        case (INT_15_POINT_RULE)
            n = 15
            call dqk15(local_fcn, a, b, rst, abserr, resabs, resasc)
        case (INT_21_POINT_RULE)
            n = 21
            call dqk21(local_fcn, a, b, rst, abserr, resabs, resasc)
        case (INT_31_POINT_RULE)
            n = 31
            call dqk31(local_fcn, a, b, rst, abserr, resabs, resasc)
        case (INT_41_POINT_RULE)
            n = 41
            call dqk41(local_fcn, a, b, rst, abserr, resabs, resasc)
        case (INT_51_POINT_RULE)
            n = 51
            call dqk51(local_fcn, a, b, rst, abserr, resabs, resasc)
        case (INT_61_POINT_RULE)
            n = 61
            call dqk61(local_fcn, a, b, rst, abserr, resabs, resasc)
        end select

        ! Report out the integration statistics
        if (present(info)) then
            info%error_estimate = abserr
            info%evaluation_count = n
            info%subinterval_count = 0
        end if
    contains
        ! This is required as the F77 code doesn't deal with procedure pointers
        function local_fcn(xarg) result(frst)
            real(real64), intent(in) :: xarg
            real(real64) :: frst
            frst = fcn(xarg)
        end function
    end function

! ------------------------------------------------------------------------------
    pure module function ni_get_rule(this) result(x)
        class(nonadaptive_integrator), intent(in) :: this
        integer(int32) :: x
        x = this%m_rule
    end function

! ------------------------------------------------------------------------------
    module subroutine ni_set_rule(this, x, err)
        class(nonadaptive_integrator), intent(inout) :: this
        integer(int32), intent(in) :: x
        class(errors), intent(inout), optional, target :: err
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        if (x /= INT_15_POINT_RULE .and. x /= INT_21_POINT_RULE .and. &
            x /= INT_31_POINT_RULE .and. x /= INT_41_POINT_RULE .and. &
            x /= INT_51_POINT_RULE .and. x /= INT_61_POINT_RULE) then
                call errmgr%report_error("ni_set_rule", &
                    "Unknown integration rule requested.", &
                    INT_INVALID_INPUT_ERROR)
                return
        end if
        this%m_rule = x
    end subroutine

! ------------------------------------------------------------------------------
end submodule
