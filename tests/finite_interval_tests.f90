! finite_interval_tests.f90

module finite_interval_tests
    use iso_fortran_env
    use integral_core
contains
    ! This test computes the integral of f(x) = exp(x) / (x**2 + 1) over the inteval
    ! of 0 -> 1.  The result should be 1.27072413983362.
    function integral_test_1() result(rst)
        ! Parameters
        real(real64), parameter :: a = 0.0d0
        real(real64), parameter :: b = 1.0d0
        real(real64), parameter :: ans = 1.27072413983362d0
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        logical :: rst
        procedure(integrand), pointer :: fcn
        type(adaptive_integrator) :: integrator
        real(real64) :: y

        ! Initialization
        rst = .true.
        fcn => integrand_1

        ! Compute the integral
        y = integrator%integrate(fcn, a, b)

        ! Check the answer
        if (abs(y - ans) > tol) then
            rst = .false.
            print '(AF8.5AF8.5A)', "TEST 1 FAILED: Expected: ", ans, &
                ", but computed: ", y, "."
        end if
    end function

    ! This test computes the integral of f(x) = exp(x) / (x**2 + 1) over the
    ! interval of 0 -> 1.  The result should be 1.27072413983362.  The
    ! difference in this test is that there are defined breakpoints in the
    ! interval.
    function integral_test_2() result(rst)
        ! Parameters
        real(real64), parameter :: a = 0.0d0
        real(real64), parameter :: b = 1.0d0
        real(real64), parameter :: ans = 1.27072413983362d0
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        logical :: rst
        procedure(integrand), pointer :: fcn
        type(adaptive_integrator) :: integrator
        real(real64) :: y, pts(3)

        ! Initialization
        rst = .true.
        fcn => integrand_1
        pts = [0.25d0, 0.5d0, 0.75d0]

        call integrator%set_use_breakpoints(.true.)
        call integrator%set_breakpoints(pts)

        ! Compute the integral
        y = integrator%integrate(fcn, a, b)

        ! Check the answer
        if (abs(y - ans) > tol) then
            rst = .false.
            print '(AF8.5AF8.5A)', "TEST 1 FAILED: Expected: ", ans, &
                ", but computed: ", y, "."
        end if
    end function

! ------------------------------------------------------------------------------
    function integrand_1(x) result(f)
        use iso_fortran_env
        real(real64), intent(in) :: x
        real(real64) :: f

        f = exp(x) / (x**2 + 1.0d0)
    end function
end module
