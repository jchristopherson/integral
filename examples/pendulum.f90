! pendulum.f90

program example
    use iso_fortran_env
    use integral_core
    use fplot_core
    implicit none

    ! Local Variables
    procedure(ode_fcn), pointer :: ptr
    type(ode_helper) :: fcn
    type(ode_auto) :: integrator
    type(plot_2d) :: plt
    type(plot_data_2d) :: d1, d2
    class(plot_axis), pointer :: xAxis, yAxis
    type(legend), pointer :: lgnd
    real(real64) :: t(2), ic(4)
    real(real64), allocatable, dimension(:,:) :: s
    real(real64), allocatable, dimension(:) :: x1, y1, x2, y2

    ! Model Parameters
    real(real64), parameter :: g = 9.81d0
    real(real64), parameter :: m1 = 2.0d-1
    real(real64), parameter :: m2 = 3.0d0
    real(real64), parameter :: l1 = 1.5d0
    real(real64), parameter :: l2 = 0.75d0

    ! Set up the integrator
    ptr => eqns
    call fcn%define_equations(4, ptr)

    ! Compute the solution
    t = [0.0d0, 40.0d0]
    ic = [0.2d0, 0.0d0, 0.8d0, 0.0d0]
    s = integrator%integrate(fcn, t, ic)

    ! Compute the X, Y positions of each mass
    x1 = l1 * cos(s(:,2))
    y1 = l1 * sin(s(:,2))

    x2 = x1 + l2 * cos(s(:,4))
    y2 = y1 + l2 * sin(s(:,4))

    ! Plot the solution
    call plt%initialize()
    call plt%set_font_size(14)

    lgnd => plt%get_legend()
    call lgnd%set_horizontal_position(LEGEND_LEFT)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)

    xAxis => plt%get_x_axis()
    call xAxis%set_title("t")

    yAxis => plt%get_y_axis()
    call yAxis%set_title("Angle [rad]")

    call d1%set_name("Mass 1")
    call d1%set_line_color(CLR_BLUE)
    call d1%define_data(s(:,1), s(:,2))

    call d2%set_name("Mass 2")
    call d2%set_line_color(CLR_GREEN)
    call d2%define_data(s(:,1), s(:,4))

    call plt%push(d1)
    call plt%push(d2)
    call plt%draw()

    ! ----------
    call plt%clear_all()
    call lgnd%set_horizontal_position(LEGEND_RIGHT)
    call lgnd%set_vertical_position(LEGEND_TOP)
    call lgnd%set_draw_inside_axes(.false.)
    call lgnd%set_draw_border(.false.)
    call xAxis%set_title("")
    call yAxis%set_title("")

    call d1%define_data(y1, x1)
    call d2%define_data(y2, x2)

    call plt%push(d1)
    call plt%push(d2)
    call plt%draw()

contains
    ! This is the routine containing the systems of equations describing the
    ! motion of a double pendulum.  See
    ! http://scienceworld.wolfram.com/physics/DoublePendulum.html for the
    ! derivation of these equations.
    !
    ! (m1 + m2) l1 x1" + m2 l2 x2" cos(x1 - x2) + m2 l2 x2'**2 sin(x1 - x2) + g (m1 + m2) sin(x1) = 0
    ! m2 l2 x2" + m2 l1 * x1" cos(x1 - x2) - m2 l1 x1'**2 sin(x1 - x2) + m2 g sin(x2) = 0
    subroutine eqns(t, x, dxdt)
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: dxdt

        ! Variables
        real(real64) :: theta1, theta2, dtheta1, dtheta2

        ! Extract variables
        theta1 = x(1)
        dtheta1 = x(2)
        theta2 = x(3)
        dtheta2 = x(4)

        ! Equations of Motion
        dxdt(1) = dtheta1
        dxdt(2) = -(l2 * m2 * dtheta2**2 * sin(theta2 - theta1) + &
            l1 * m2 * dtheta1**2 * cos(theta2 - theta1) * &
            sin(theta2 - theta1) + g * m2 * sin(theta2) * &
            cos(theta2 - theta1) - g * (m1 + m1) * sin(theta1)) / &
            (l1 * (m2 * (cos(theta2 - theta1)**2 - 1) - m1))
        dxdt(3) = dtheta2
        dxdt(4) = (l2 * m2 * cos(theta2 - theta1) * sin(theta2 - theta1) * &
            dtheta2**2 + l1 * dtheta1**2 * (m2 + m1) * sin(theta2 - theta1) - &
            g * (m2 + m1) * sin(theta1) * cos(theta2 - theta1) + &
            g * (m2 + m1) * sin(theta2)) / &
            (l2 * (m2 * (cos(theta2 - theta1)**2 - 1) - m1))
    end subroutine
end program
