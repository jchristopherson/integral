! infinite_interval_tests.f90

module infinite_interval_tests
    use iso_fortran_env
    use integral_core
contains
    ! This test computes the integral of f(x) = 1 / x**2 over the interval
    ! of [1, inf).  The result should be 1.
    function inf_integral_test_1() result(rst)
        ! Parameters
        real(real64), parameter :: ans = 1.0d0
        real(real64), parameter :: tol = 1.0d-8
        real(real64), parameter :: bounds = 1.0d0

        ! Local Variables
        logical :: rst
        procedure(integrand), pointer :: fcn
        real(real64) :: y
        type(infinite_interval_integrator) :: integrator

        ! Initialization
        rst = .true.
        fcn => integrand_1

        ! Compute the integral
        y = integrator%integrate(fcn, bounds, .false.)

        ! Check the answer
        if (abs(y - ans) > tol) then
            rst = .false.
            print '(AF0.5AF0.5A)', &
                "INFINITE INTERVAL TEST 1 FAILED: Expected: ", ans, &
                ", but computed: ", y, "."
        end if
    end function

    ! This test computes the integral of f(x) = 1 / (1 + x**2) over an interval
    ! of (-inf, inf).  The result should be pi.
    function inf_integral_test_2() result(rst)
        ! Parameters
        real(real64), parameter :: ans = 2.0d0 * acos(0.0d0)
        real(real64), parameter :: tol = 1.0d-8

        ! Local Variables
        logical :: rst
        procedure(integrand), pointer :: fcn
        real(real64) :: y
        type(infinite_interval_integrator) :: integrator

        ! Initialization
        rst = .true.
        fcn => integrand_2

        ! Compute the integral
        y = integrator%integrate(fcn)

        ! Check the answer
        if (abs(y - ans) > tol) then
            rst = .false.
            print '(AF0.5AF0.5A)', &
                "INFINITE INTERVAL TEST 2 FAILED: Expected: ", ans, &
                ", but computed: ", y, "."
        end if
    end function

! ------------------------------------------------------------------------------
    function integrand_1(x) result(f)
        real(real64), intent(in) :: x
        real(real64) :: f
        f = 1.0d0 / x**2
    end function

    function integrand_2(x) result(f)
        real(real64), intent(in) :: x
        real(real64) :: f
        f = 1.0d0 / (1.0d0 + x**2)
    end function
end module
