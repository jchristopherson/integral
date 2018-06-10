! integral_example_1.f90

! This example is from: http://tutorial.math.lamar.edu/Classes/CalcI/ComputingDefiniteIntegrals.aspx#Int_CompDef_Ex3a
program example
    use iso_fortran_env
    use integral_core
    implicit none

    ! Variables
    real(real64) :: ans, y, pi, a, b
    procedure(integrand), pointer :: fcn
    type(adaptive_integrator) :: integrator

    ! Define the integration limits
    pi = 2.0d0 * acos(0.0d0)
    a = pi / 6.0d0
    b = pi / 4.0d0

    ! Evaluate the integral
    fcn => int_fcn
    y = integrator%integrate(fcn, a, b)

    ! Display the results
    ans = 5.0d0 * pi / 12.0d0 - 2.0d0 * sqrt(2.0d0) + 4.0d0 / sqrt(3.0d0)
    print '(AEN13.5AEN13.5A)', "The solution is: ", ans, &
        ", the integrator computed: ", y, "."

contains
    ! This example is from http://tutorial.math.lamar.edu/Classes/CalcI/ComputingDefiniteIntegrals.aspx#Int_CompDef_Ex3a
    ! The integrand is: f(x) = 5 - 2 sec(x) tan(x).
    ! If the integral is considered over the range [pi/6, pi/4], the solution
    ! is 5 pi / 12 - 2 sqrt(2) + 4 / sqrt(3).
    function int_fcn(x) result(f)
        real(real64), intent(in) :: x
        real(real64) :: f
        f = 5.0d0 - 2.0d0 * tan(x) / cos(x) ! Remember, sec(x) = 1 / cos(x)
    end function
end program
