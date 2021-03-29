! infinite_integral_example_2.f90

program main
    use iso_fortran_env
    use integral_core
    implicit none

    ! Variables
    type(infinite_interval_integrator) :: integrator
    procedure(integrand), pointer :: fptr
    real(real64) :: bounds, y
    logical :: isUpperLimit

    ! Define the function to integrate
    fptr => fcn

    ! Define the non-infinite integration limit
    bounds = 1.0d0
    isUpperLimit = .false.

    ! Compute the integral using the limits [1, infinity)
    y = integrator%integrate(fptr, bounds, isUpperLimit)
    print '(AF0.1A)', "Result: ", y, ", Expected: 1.0"

contains
    ! The function to integrate:
    ! f(x) = 1 / x**2
    function fcn(x) result(rst)
        real(real64), intent(in) :: x
        real(real64) :: rst

        rst = 1.0d0 / x**2
    end function
end program
