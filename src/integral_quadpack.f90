! integral_quadpack.f90

! REF:
! - http://www.netlib.org/quadpack/
! - https://people.sc.fsu.edu/~jburkardt/f_src/quadpack/quadpack.html
submodule (integral_core) integral_quadpack
contains
    module function integrate_adapt(f, a, b, cntrls, bhvr, err) result(rst)
        ! Arguments
        procedure(integrand), pointer, intent(in) :: f
        real(real64), intent(in) :: a, b
        type(integration_controls), intent(in), optional :: cntrls
        type(integration_behavior), intent(out), optional :: bhvr
        class(errors), intent(inout), optional, target :: err
        real(real64) :: rst

        ! Local Variables
        integer(int32) :: ier, neval, limit, lenw, last, flag
        integer(int32), allocatable, dimension(:) :: iwork
        real(real64) :: epsabs, epsrel, abserr, epsmch
        real(real64), allocatable, dimension(:) :: work
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        epsmch = sqrt(epsilon(epsmch))
        if (present(cntrls)) then
            epsabs = cntrls%absolute_tolerance
            epsrel = cntrls%relative_tolerance
            limit = cntrls%max_subinterval_count
        else
            epsabs = epsmch
            epsrel = epsmch
            limit = 50
        end if

        ! Allocate the workspaces
        lenw = 4 * limit
        allocate(iwork(limit), stat = flag)
        if (flag == 0) allocate(work(lenw), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("integrate_adapt", &
                "Insufficient memory available.", INT_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Call the integrator
        call dqags(local_fcn, a, b, epsabs, epsrel, rst, abserr, neval, ier, &
            limit, lenw, last, iwork, work)

        ! Collect the integration results
        if (present(bhvr)) then
            bhvr%error_estimate = abserr
            bhvr%evaluation_count = neval
            bhvr%subinterval_count = last
        end if

        ! Check for errors
        select case (ier)
        case (1)
            write(errmsg, '(AI0A)') "The maximum number (", limit, &
                ") of subdivisions allowed has been achieved."
            call errmgr%report_error("integrate_adapt", trim(errmsg), &
                INT_COUNT_EXCEEDED_ERROR)
            return
        case (2)
            call errmgr%report_error("integrate_adapt", &
                "The occurence of roundoff error has been detected " // &
                "thereby preventing the requested tolerance from being " // &
                "achieved.", INT_ROUND_OFF_ERROR)
            return
        case (3)
            call errmgr%report_error("integrate_adapt", &
                "The integrand appears to be behaving too poorly to proceed.", &
                INT_INTEGRAND_BEHAVIOR_ERROR)
            return
        case (4)
            call errmgr%report_error("integrate_adapt", &
                "The algorithm cannot converge with the assigned tolerances.", &
                INT_CONVERGENCE_ERROR)
            return
        case (5)
            call errmgr%report_error("integrate_adapt", &
                "The integral is likely divergent", INT_DIVERGENT_ERROR)
            return
        case (6)
            call errmgr%report_error("integrate_adapt", &
                "An invalid input has been provided.", INT_INVALID_INPUT_ERROR)
            return
        end select

    contains
        ! This is required as the F77 code doesn't deal with procedure pointers
        function local_fcn(xarg) result(frst)
            real(real64), intent(in) :: xarg
            real(real64) :: frst
            frst = f(xarg)
        end function
    end function
end submodule
