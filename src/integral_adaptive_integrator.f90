! integral_adaptive_integrator.f90

submodule (integral_core) integral_adaptive_integrator
contains
! ------------------------------------------------------------------------------
    module subroutine ai_init(this, err)
        ! Arguments
        class(adaptive_integrator), intent(inout) :: this
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        integer(int32) :: flag, liwork, lwork, n
        logical :: use_breakpoints

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        use_breakpoints = .false.
        if (this%get_use_breakpoints()) then
            if (allocated(this%m_breakpoints)) use_breakpoints = .true.
        end if

        ! Allocate memory
        if (use_breakpoints) then
            n = size(this%m_breakpoints) + 2
            liwork = max(2 * this%get_max_subintervals() + n, 3 * n - 2)
            lwork = liwork * 2 - n
        else
            liwork = this%get_max_subintervals()
            lwork = 4 * liwork
        end if

        if (.not.allocated(this%m_iwork)) then
            allocate(this%m_iwork(liwork), stat = flag)
        else
            if (size(this%m_iwork) < liwork) then
                deallocate(this%m_iwork)
                allocate(this%m_iwork(liwork), stat = flag)
            end if
        end if
        if (flag == 0) then
            if (.not.allocated(this%m_work)) then
                allocate(this%m_work(lwork), stat = flag)
            else
                if (size(this%m_work) < lwork) then
                    deallocate(this%m_work)
                    allocate(this%m_work(lwork), stat = flag)
                end if
            end if
        end if

        ! Report any errors
        if (flag /= 0) then
            call errmgr%report_error("ai_init", &
                "Insufficient memory available.", INT_OUT_OF_MEMORY_ERROR)
            return
        end if
    end subroutine

! ------------------------------------------------------------------------------
    module function ai_integrate(this, fcn, a, b, info, err) result(rst)
        ! Arguments
        class(adaptive_integrator), intent(inout) :: this
        procedure(integrand), intent(in), pointer :: fcn
        real(real64), intent(in) :: a, b
        type(integration_behavior), intent(out), optional :: info
        class(errors), intent(inout), optional, target :: err
        real(real64) :: rst

        ! Local Variables
        integer(int32) :: ier, neval, limit, lenw, last, leniw, npts, flag
        real(real64) :: abserr, epsabs, epsrel
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg
        logical :: usedqagp
        real(real64), allocatable, dimension(:) :: pts

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        usedqagp = .false.
        if (this%get_use_breakpoints()) then
            if (allocated(this%m_breakpoints)) then
                npts = size(this%m_breakpoints) + 2
                allocate(pts(npts), stat = flag)
                if (flag /= 0) then
                    call errmgr%report_error("ai_integrate", &
                        "Insufficient memory available.", &
                        INT_OUT_OF_MEMORY_ERROR)
                    return
                end if
                pts(2:npts-1) = this%m_breakpoints
                pts(1) = a
                pts(npts) = b
                usedqagp = .true.
            end if
        end if

        call this%initialize(errmgr)
        if (errmgr%has_error_occurred()) return

        limit = this%get_max_subintervals()
        leniw = size(this%m_iwork)
        lenw = size(this%m_work)
        epsabs = this%get_abs_tol()
        epsrel = this%get_rel_tol()

        ! Compute the integral
        if (usedqagp) then
            call dqagp(local_fcn, a, b, npts, pts, epsabs, epsrel, &
                rst, abserr, neval, ier, leniw, lenw, last, iwork, work)
        else
            call dqags(local_fcn, a, b, epsabs, epsrel, rst, abserr, neval, &
                ier, limit, lenw, last, this%m_iwork, this%m_work)
        end if

        ! Collect info regarding the integration process
        if (present(info)) then
            info%error_estimate = abserr
            info%evaluation_count = neval
            info%subinterval_count = last
        end if

        ! Check for errors
        select case (ier)
        case (1)
            write(errmsg, '(AI0A)') "The maximum number (", limit, &
                ") of subdivisions allowed has been achieved."
            call errmgr%report_error("ai_integrate", trim(errmsg), &
                INT_COUNT_EXCEEDED_ERROR)
            return
        case (2)
            call errmgr%report_error("ai_integrate", &
                "The occurence of roundoff error has been detected " // &
                "thereby preventing the requested tolerance from being " // &
                "achieved.", INT_ROUND_OFF_ERROR)
            return
        case (3)
            call errmgr%report_error("ai_integrate", &
                "The integrand appears to be behaving too poorly to proceed.", &
                INT_INTEGRAND_BEHAVIOR_ERROR)
            return
        case (4)
            call errmgr%report_error("ai_integrate", &
                "The algorithm cannot converge with the assigned tolerances.", &
                INT_CONVERGENCE_ERROR)
            return
        case (5)
            call errmgr%report_error("ai_integrate", &
                "The integral is likely divergent", INT_DIVERGENT_ERROR)
            return
        case (6)
            call errmgr%report_error("ai_integrate", &
                "An invalid input has been provided.", INT_INVALID_INPUT_ERROR)
            return
        end select

    contains
        ! This is required as the F77 code doesn't deal with procedure pointers
        function local_fcn(xarg) result(frst)
            real(real64), intent(in) :: xarg
            real(real64) :: frst
            frst = fcn(xarg)
        end function
    end function

! ------------------------------------------------------------------------------
    pure module function ai_get_use_brkpnts(this) result(x)
        class(adaptive_integrator), intent(in) :: this
        logical :: x
        x = this%m_userDefinedBreaks
    end function

! --------------------
    module subroutine ai_set_use_brkpnts(this, x)
        class(adaptive_integrator), intent(inout) :: this
        logical, intent(in) :: x
        this%m_userDefinedBreaks = x
    end subroutine

! ------------------------------------------------------------------------------
    module function ai_get_breakpoints(this) result(x)
        class(adaptive_integrator), intent(in) :: this
        real(real64), allocatable, dimension(:) :: x
        if (allocated(this%m_breakpoints)) then
            x = this%m_breakpoints
        end if
    end function

! --------------------
    module subroutine ai_set_breakpoints(this, x)
        class(adaptive_integrator), intent(inout) :: this
        real(real64), intent(in), dimension(:) :: x
        if (allocated(this%m_breakpoints)) deallocate(this%m_breakpoints)
        this%m_breakpoints = x
        call this%set_use_breakpoints(.true.)
    end subroutine

! ------------------------------------------------------------------------------
end submodule
