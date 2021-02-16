! infinite_integral_example_1.f90

program main
    use iso_fortran_env
    use integral_core
    implicit none

    ! Variables
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    type(infinite_interval_integrator) :: integrator
    procedure(integrand), pointer :: fptr
    real(real64) :: y

    ! Define the function to integrate
    fptr => fcn

    ! Compute the integral using the limits (-infinity, infinity)
    y = integrator%integrate(fptr)
    print '(AF0.6AF0.6)', "Result: ", y, ", Expected: ", pi

contains
    ! The function to integrate:
    ! f(x) = 1 / (1 + x**2)
    function fcn(x) result(rst)
        real(real64), intent(in) :: x
        real(real64) :: rst

        rst = 1.0d0 / (1.0d0 + x**2)
    end function
end program
