! vanderpol.f90

program example
    use iso_fortran_env
    use integral_core
    use fplot_core
    implicit none

    ! Local Variables
    type(ode_helper) :: fcn
    type(ode_irk) :: integrator1
    type(ode_auto) :: integrator2
    procedure(ode_fcn), pointer :: ptr
    real(real64) :: ic(2), t(2)
    real(real64), allocatable, dimension(:,:) :: x1, x2
    type(plot_2d) :: plt
    type(plot_data_2d) :: d1, d2
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd

    ! Set up the integrator
    ptr => vdp
    call fcn%define_equations(2, ptr)

    ! Define the initial conditions
    t = [0.0d0, 8.0d1]
    ic = [2.0d0, 0.0d0]

    ! Compute the solution
    x1 = integrator1%integrate(fcn, t, ic)
    x2 = integrator2%integrate(fcn, t, ic)

    ! Display the number of solution points in each
    print '(AI0)', "ODE_IRK Solution Point Count: ", size(x1, 1)
    print '(AI0)', "ODE_AUTO Solution Point Count: ", size(x2, 1)

    ! Plot the solution
    call plt%initialize()
    call plt%set_font_size(14)

    xAxis => plt%get_x_axis()
    call xAxis%set_title("t")

    yAxis => plt%get_y_axis()
    call yAxis%set_title("x(t)")

    lgnd => plt%get_legend()
    call lgnd%set_is_visible(.true.)
    call lgnd%set_draw_border(.false.)
    call lgnd%set_draw_inside_axes(.false.)

    call d1%set_name("IRK")
    call d1%set_draw_line(.false.)
    call d1%set_draw_markers(.true.)
    call d1%set_marker_style(MARKER_FILLED_TRIANGLE)
    call d1%set_marker_scaling(1.5)
    call d1%define_data(x1(:,1), x1(:,2))
    call plt%push(d1)

    call d2%set_name("AUTO")
    call d2%set_draw_line(.false.)
    call d2%set_draw_markers(.true.)
    call d2%set_marker_style(MARKER_EMPTY_CIRCLE)
    call d2%set_line_color(CLR_RED)
    call d2%set_line_style(LINE_DASHED)
    call d2%define_data(x2(:,1), x2(:,2))
    call plt%push(d2)

    call plt%draw()

    call plt%clear_all()

    call d1%define_data(x1(:,2), x1(:,3))
    call d2%define_data(x2(:,2), x2(:,3))
    call xAxis%set_title("x(t)")
    call yAxis%set_title("dx/dt")
    call plt%push(d1)
    call plt%push(d2)
    call plt%draw()

contains
    ! Van Der Pol Equation
    ! x" + x - mu * (1 - x**2) * x' = 0
    subroutine vdp(t, x, dxdt)
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: dxdt

        real(real64), parameter :: mu = 20.0d0

        dxdt(1) = x(2)
        dxdt(2) = mu * (1.0d0 - x(1)**2) * x(2) - x(1)
    end subroutine
end program
