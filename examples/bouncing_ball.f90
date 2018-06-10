! bouncing_ball.f90

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
