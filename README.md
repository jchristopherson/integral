# integral
The INTEGRAL library provides routines for the integration of functions of various types.  Additionally, the INTEGRAL library provides routines for the integration of systems of ordinary differential equations (ODEs).

The integration routines are provided by [QUADPACK](http://www.netlib.org/quadpack/), and the ODE routines are provided by [ODEPACK](http://www.netlib.org/odepack/).

## Example 1
The following example illustrates the use of an adaptive integrator to compute the integral of an equation over a finite interval.
```fortran
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
```
```text
The solution is: 789.97089E-03, the integrator computed: 789.97089E-03.
```

## Example 2
The following example illustrates how to compute the solution to a system of ODEs modeling the bouncing of a ball.  The example also utilizes the [FPLOT](https://github.com/jchristopherson/fplot) library in order to plot the solution.
```fortran
program example
    use iso_fortran_env
    use integral_core
    use fplot_core
    implicit none

    ! Parameters
    real(real64), parameter :: g = 9.81d0 ! Gravitational acceleration
    real(real64), parameter :: k = -0.8d0 ! Coefficient of restitution

    ! Local Variables
    procedure(ode_fcn), pointer :: ptr
    procedure(ode_constraint), pointer :: cptr
    type(ode_helper) :: fcn
    type(ode_auto) :: integrator
    integer(int32) :: n
    real(real64) :: ic(2), t(2)
    real(real64), allocatable, dimension(:,:) :: x1, x2, x3, x4
    type(plot_2d) :: plt
    type(plot_data_2d) :: d1, d2, d3, d4
    class(plot_axis), pointer :: xAxis, yAxis
    type(legend), pointer :: lgnd

    ! Set up the integrator
    ptr => ball
    cptr => ground_constraint
    call fcn%define_equations(2, ptr)
    call fcn%define_constraints(1, cptr)
    call integrator%set_max_step_size(1.0d-3)
    call integrator%set_limit_step_size(.true.)

    ! Compute the solution
    t = [0.0d0, 1.0d1]
    ic = [1.0d1, 5.0d0]
    x1 = integrator%integrate(fcn, t, ic)

    ! The integrator stops when the ball first makes contact.  As a result, lets
    ! reset the time limits and initial conditions to continue the integration
    n = size(x1, 1)
    t(1) = x1(n,1)
    ic = [abs(x1(n,2)), k * x1(n,3)]
    call integrator%reset()
    x2 = integrator%integrate(fcn, t, ic)

    ! Again
    n = size(x2, 1)
    t(1) = x2(n,1)
    ic = [abs(x2(n,2)), k * x2(n,3)]
    call integrator%reset()
    x3 = integrator%integrate(fcn, t, ic)

    ! Again
    n = size(x3, 1)
    t(1) = x3(n,1)
    ic = [abs(x3(n,2)), k * x3(n,3)]
    call integrator%reset()
    x4 = integrator%integrate(fcn, t, ic)


    ! Plot the solution
    call plt%initialize()
    call plt%set_font_size(14)

    lgnd => plt%get_legend()
    call lgnd%set_is_visible(.false.)

    xAxis => plt%get_x_axis()
    call xAxis%set_title("t")

    yAxis => plt%get_y_axis()
    call yAxis%set_title("x(t)")

    call d1%set_line_color(CLR_BLUE)
    call d1%set_line_width(2.0)
    call d1%define_data(x1(:,1), x1(:,2))

    call d2%set_line_color(CLR_BLUE)
    call d2%set_line_width(2.0)
    call d2%define_data(x2(:,1), x2(:,2))

    call d3%set_line_color(CLR_BLUE)
    call d3%set_line_width(2.0)
    call d3%define_data(x3(:,1), x3(:,2))

    call d4%set_line_color(CLR_BLUE)
    call d4%set_line_width(2.0)
    call d4%define_data(x4(:,1), x4(:,2))

    call plt%push(d1)
    call plt%push(d2)
    call plt%push(d3)
    call plt%push(d4)
    call plt%draw()

contains
    ! A bouncing ball can be described by the following equation:
    ! x" = -g
    !
    ! Where g = gravitational acceleration
    subroutine ball(t, x, dxdt)
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: dxdt
        dxdt(1) = x(2)
        dxdt(2) = -g
    end subroutine

    ! The constraint function
    subroutine ground_constraint(t, x, f)
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: f
        f(1) = x(1)     ! Find when x == 0
    end subroutine
end program
```
This is the plot resulting from the above program.
![](images/bouncing_ball_example.png?raw=true)

## Documentation
Documentation can be found [here](http://htmlpreview.github.io/?https://github.com/jchristopherson/integral/blob/master/doc/html/index.html).

## Build Instructions
This library utilizes [CMake](https://cmake.org/) to facilitate its build.  Using CMake is as simple as issuing the following commands.
- cmake ...
- make
- make install

## Dependencies
This library depends upon the following libraries.
- [FERROR](https://github.com/jchristopherson/ferror)

See [Running CMake](https://cmake.org/runningcmake/) for more details on the use of CMake.
