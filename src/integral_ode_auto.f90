! integral_ode_auto.f90

submodule (integral_core) integral_ode_auto
contains
! ------------------------------------------------------------------------------
    module function oa_step(this, fcn, x, y, xout, rtol, atol, err) result(brk)
        ! Arguments
        class(ode_auto), intent(inout) :: this
        class(ode_helper), intent(inout) :: fcn
        real(real64), intent(inout) :: x
        real(real64), intent(inout), dimension(:) :: y
        real(real64), intent(in) :: xout
        real(real64), intent(in), dimension(:) :: rtol, atol
        class(errors), intent(inout), optional, target :: err
        logical :: brk

        ! Local Variables
        integer(int32) :: neq, itask, iopt, lrw, liw, jt, ncnst, lrn, lrs, itol
        logical :: useConstraints
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        brk = .false.
        neq = fcn%get_equation_count()
        ncnst = fcn%get_constraint_count()
        useConstraints = fcn%get_constraints_defined()
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Ensure workspace arrays are sized appropriately
        lrn = 20 + 16 * neq + 3 * ncnst
        lrs = 22 + 9 * neq + neq**2 + 3 * ncnst
        lrw = max(lrn, lrs)
        liw = 20 + neq

        ! Define workspace arrays
        call this%init_workspace(liw, lrw, ncnst, errmgr)
        if (errmgr%has_error_occurred()) return

        ! Additional Initialization
        iopt = 0
        if (this%get_allow_overshoot()) then
            itask = 1
            if (this%get_provide_all_output()) itask = 2
        else
            itask = 4
            iopt = 1
            this%m_rwork(1) = this%get_integration_limit()
            if (this%get_provide_all_output()) itask = 5
        end if
        if (fcn%get_jacobian_defined()) then
            jt = 1
        else
            jt = 2
        end if
        itol = 4

        if (this%get_limit_step_size()) then
            iopt = 1
            this%m_rwork(6) = this%get_max_step_size()
        end if

        if (this%get_iteration_limit() /= 500) then
            iopt = 1
            this%m_iwork(6) = this%get_iteration_limit()
        end if

        ! Determine if we are to call DLSODA or DLSODAR
        if (useConstraints .and. ncnst > 0) then
            call DLSODAR(odepackfcn, neq, y, x, xout, itol, rtol, atol, itask, &
                this%m_istate, iopt, this%m_rwork, lrw, this%m_iwork, liw, &
                odepackjac, jt, odepackconstraints, ncnst, this%m_rststats)

            ! Check output condition
            select case (this%m_istate)
            case (3)
                ! Root Found
                brk = .true.
            case (-1)
                call errmgr%report_error("oa_step", &
                    "Excessive work has been done during this " // &
                    "integration step.", INT_EXCESSIVE_WORK_ERROR)
            case (-2)
                call errmgr%report_error("oa_step", &
                    "Unachievable tolerances have been requested.", &
                    INT_IMPRACTICAL_TOLERANCE_ERROR)
            case (-3)
                call errmgr%report_error("oa_step", &
                    "An invalid input was encountered.", &
                    INT_INVALID_INPUT_ERROR)
            case (-4)
                call errmgr%report_error("oa_step", &
                    "Error testing has failed repeatadly.  " // &
                    "Ensure inputs are valid.", INT_REPEATED_ERROR_TEST_FAILURE)
            case (-5)
                call errmgr%report_error("oa_step", "Repeated convergence " // &
                    "failures have been encountered.", INT_CONVERGENCE_ERROR)
            case (-6)
                call errmgr%report_error("oa_step", "A solution component " // &
                    "has effectively vanished.", INT_INTEGRAND_BEHAVIOR_ERROR)
            case (-7)
                call errmgr%report_error("oa_step", "An insufficiently " // &
                    "sized workspace array was provided.", &
                    INT_INVALID_INPUT_ERROR)
            end select
        else
            call DLSODA(odepackfcn, neq, y, x, xout, itol, rtol, atol, itask, &
                this%m_istate, iopt, this%m_rwork, lrw, this%m_iwork, liw, &
                odepackjac, jt)

            ! Check output condition
            select case (this%m_istate)
            case (-1)
                call errmgr%report_error("oa_step", &
                    "Excessive work has been done during this " // &
                    "integration step.", INT_EXCESSIVE_WORK_ERROR)
            case (-2)
                call errmgr%report_error("oa_step", &
                    "Unachievable tolerances have been requested.", &
                    INT_IMPRACTICAL_TOLERANCE_ERROR)
            case (-3)
                call errmgr%report_error("oa_step", &
                    "An invalid input was encountered.", &
                    INT_INVALID_INPUT_ERROR)
            case (-4)
                call errmgr%report_error("oa_step", &
                    "Error testing has failed repeatadly.  " // &
                    "Ensure inputs are valid.", INT_REPEATED_ERROR_TEST_FAILURE)
            case (-5)
                call errmgr%report_error("oa_step", "Repeated convergence " // &
                    "failures have been encountered.", INT_CONVERGENCE_ERROR)
            case (-6)
                call errmgr%report_error("oa_step", "A solution component " // &
                    "has effectively vanished.", INT_INTEGRAND_BEHAVIOR_ERROR)
            case (-7)
                call errmgr%report_error("oa_step", "An insufficiently " // &
                    "sized workspace array was provided.", &
                    INT_INVALID_INPUT_ERROR)
            end select
        end if

    contains
        subroutine odepackfcn(n, x, y, dydx)
            integer(int32), intent(in) :: n
            real(real64), intent(in) :: x, y(n)
            real(real64), intent(out) :: dydx(n)
            call fcn%evaluate_ode(x, y, dydx)
        end subroutine

        subroutine odepackjac(n, x, y, ml, mu, pd, nrowpd)
            integer(int32), intent(in) :: n, ml, mu, nrowpd
            real(real64), intent(in) :: x, y(n)
            real(real64), intent(out) :: pd(nrowpd, n)
            call fcn%evaluate_jacobian(x, y, pd)
        end subroutine

        subroutine odepackconstraints(neq, x, y, ng, gout)
            integer(int32), intent(in) :: neq, ng
            real(real64), intent(in) :: x, y(neq)
            real(real64), intent(out) :: gout(ng)
            call fcn%evaluate_constraints(x, y, gout)
        end subroutine
    end function

! ------------------------------------------------------------------------------
    module subroutine oa_reset_integrator(this)
        class(ode_auto), intent(inout) :: this
        this%m_istate = 1
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine oa_init_workspace(this, liw, lrw, ncnsts, err)
        ! Arguments
        class(ode_auto), intent(inout) :: this
        integer(int32), intent(in) :: liw, lrw, ncnsts
        class(errors), intent(inout) :: err

        ! Local Variables
        integer(int32) :: flag

        ! Process
        flag = 0
        if (allocated(this%m_iwork)) then
            if (size(this%m_iwork) < liw) then
                deallocate(this%m_iwork)
                allocate(this%m_iwork(liw), stat = flag)
            end if
        else
            allocate(this%m_iwork(liw), stat = flag)
        end if

        if (flag == 0) then
            if (allocated(this%m_rwork)) then
                if (size(this%m_rwork) < lrw) then
                    deallocate(this%m_rwork)
                    allocate(this%m_rwork(lrw), stat = flag)
                end if
            else
                allocate(this%m_rwork(lrw), stat = flag)
            end if
        end if

        if (flag == 0 .and. ncnsts > 0) then
            if (allocated(this%m_rststats)) then
                if (size(this%m_rststats) /= ncnsts) then
                    deallocate(this%m_rststats)
                    allocate(this%m_rststats(ncnsts), stat = flag)
                end if
            else
                allocate(this%m_rststats(ncnsts), stat = flag)
            end if
        end if

        ! Error Checking
        if (flag /= 0) then
            call err%report_error("oa_init_workspace", &
                "Insufficient memory available.", INT_OUT_OF_MEMORY_ERROR)
            return
        end if
    end subroutine

! ------------------------------------------------------------------------------
    module function oa_get_constraint_info(this) result(x)
        ! Arguments
        class(ode_auto), intent(in) :: this
        logical, allocatable, dimension(:) :: x

        ! Local Variables
        integer(int32) :: i, n

        ! Process
        if (.not.allocated(this%m_rststats)) return
        n = size(this%m_rststats)
        allocate(x(n))
        do i = 1, n
            x(i) = this%m_rststats(i) == 1
        end do
    end function

! ------------------------------------------------------------------------------
end submodule
