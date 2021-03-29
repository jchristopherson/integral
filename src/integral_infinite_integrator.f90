! integral_infinite_integrator.f90

submodule (integral_core) integral_infinite_integrator
contains
! ------------------------------------------------------------------------------
    module subroutine iii_initialize(this, err)
        ! Arguments
        class(infinite_interval_integrator), intent(inout) :: this
        class(errors), intent(inout), optional, target :: err

        ! Local Variables
        integer(int32) :: flag, lwork, liwork
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Set up the error handling
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Initialization
        liwork = this%get_max_subintervals()
        lwork = 4 * this%get_max_subintervals()

        ! Memory Allocation
        if (allocated(this%m_work)) deallocate(this%m_work)
        if (allocated(this%m_iwork)) deallocate(this%m_iwork)
        allocate(this%m_work(lwork), stat = flag)
        if (flag == 0) allocate(this%m_iwork(liwork), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("iii_initialize", &
                "Insufficient memory available.", INT_OUT_OF_MEMORY_ERROR)
            return
        end if
    end subroutine

! ------------------------------------------------------------------------------
    module function iii_integrate(this, fcn, bound, upper, info, err) result(rst)
        ! Arguments
        class(infinite_interval_integrator), intent(inout) :: this
        procedure(integrand), intent(in), pointer :: fcn
        real(real64), intent(in), optional :: bound
        logical, intent(in), optional :: upper
        type(integration_behavior), intent(out), optional :: info
        class(errors), intent(inout), optional, target :: err
        real(real64) :: rst

        ! Local Variables
        integer(int32) :: inf, neval, ier, limit, lwork, last
        real(real64) :: bnd, epsabs, epsrel, abserr
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg
        
        ! Error Handling Set-Up
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Initialization
        limit = this%get_max_subintervals()
        lwork = 4 * this%get_max_subintervals()
        bnd = 0.0d0
        inf = 2
        if (present(bound)) then
            bnd = bound
            if (present(upper)) then
                if (upper) then
                    ! bound is the upper limit
                    inf = -1
                else
                    ! bound is the lower limit
                    inf = 1
                end if
            else
                ! ERROR: upper must be supplied if bound is supplied
                call errmgr%report_error("iii_integrate", &
                    "The parameter upper must be supplied when the " // &
                    "parameter bound is supplied.", &
                    INT_INVALID_OPERATION_ERROR)
                return
            end if
        end if
        
        ! Ensure the object is initialized properly
        if (.not.allocated(this%m_work) .or. .not.allocated(this%m_iwork)) then
            call this%initialize(errmgr)
            if (errmgr%has_error_occurred()) return
        end if
        if (size(this%m_work) < lwork) then
            call this%initialize(errmgr)
            if (errmgr%has_error_occurred()) return
        end if

        ! Get the tolerance info
        epsrel = this%get_rel_tol()
        epsabs = this%get_abs_tol()

        ! Compute the integral
        call dqagi(local_fcn, bnd, inf, epsabs, epsrel, rst, abserr, neval, &
            ier, limit, lwork, last, this%m_iwork, this%m_work)

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
            call errmgr%report_error("iii_integrate", trim(errmsg), &
                INT_COUNT_EXCEEDED_ERROR)
            return
        case (2)
            call errmgr%report_error("iii_integrate", &
                "The occurence of roundoff error has been detected " // &
                "thereby preventing the requested tolerance from being " // &
                "achieved.", INT_ROUND_OFF_ERROR)
            return
        case (3)
            call errmgr%report_error("iii_integrate", &
                "The integrand appears to be behaving too poorly to proceed.", &
                INT_INTEGRAND_BEHAVIOR_ERROR)
            return
        case (4)
            call errmgr%report_error("iii_integrate", &
                "The algorithm cannot converge with the assigned tolerances.", &
                INT_CONVERGENCE_ERROR)
            return
        case (5)
            call errmgr%report_error("iii_integrate", &
                "The integral is likely divergent", INT_DIVERGENT_ERROR)
            return
        case (6)
            call errmgr%report_error("iii_integrate", &
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
end submodule
