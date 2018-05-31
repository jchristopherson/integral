! integral_quadpack.f90

! REF:
! - http://www.netlib.org/quadpack/
! - https://people.sc.fsu.edu/~jburkardt/f_src/quadpack/quadpack.html
submodule (integral_core) integral_quadpack
contains
    module function integrate_adaptive(f, a, b, cntrls, bhvr, err) result(rst)
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
            call errmgr%report_error("qags", "Insufficient memory available.", &
                INT_OUT_OF_MEMORY_ERROR)
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
        case (2)
        case (3)
        case (4)
        case (5)
        case (6)
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
