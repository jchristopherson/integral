! rk_example_1.f90

program example
    use iso_fortran_env
    use integral_core
    use fplot_core
    implicit none

    ! Local Variables
    procedure(ode_fcn), pointer :: ptr
    type(ode_helper) :: fcn
    type(ode_rk45) :: integrator1
    type(ode_irk) :: integrator2
    type(plot_2d) :: plt
    type(plot_data_2d) :: d1, d2
    class(plot_axis), pointer :: xAxis, yAxis
    type(legend), pointer :: lgnd
    real(real64) :: ts(2), ic(2)
    real(real64), allocatable, dimension(:,:) :: z1, z2

    ! Set up the integrators
    ptr => eqns
    call fcn%define_equations(2, ptr)

    ! Compute the solution
    ts = [0.0d0, 1.0d2]
    ic = [0.0d0, 0.0d0]
    z1 = integrator1%integrate(fcn, ts, ic)
    z2 = integrator2%integrate(fcn, ts, ic)

    ! Display the number of solution points in each
    print '(AI0)', "ODE_RK45 Solution Point Count: ", size(z1, 1)
    print '(AI0)', "ODE_IRK Solution Point Count: ", size(z2, 1)

    ! Plot the solution
    call plt%initialize()
    call plt%set_font_size(14)

    lgnd => plt%get_legend()
    call lgnd%set_is_visible(.true.)
    call lgnd%set_vertical_position(LEGEND_BOTTOM)

    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()

    call xAxis%set_title("t")
    call yAxis%set_title("x(t)")
    call plt%set_title("Runge-Kutta Method Comparison")

    call d1%set_name("4-5 Method")
    call d1%define_data(z1(:,1), z1(:,2))
    call plt%push(d1)

    call d2%set_name("Implicit")
    call d2%define_data(z2(:,1), z2(:,2))
    call d2%set_line_color(CLR_RED)
    call d2%set_line_style(LINE_DASHED)
    call d2%set_line_width(2.0)
    call plt%push(d2)

    call plt%draw()

contains
    ! This is Duffing's equation of the form: 
    ! x" + s*x' + a*x + b*x**3 = g*sin(w * t)
    subroutine eqns(t, x, dxdt)
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: dxdt

        ! Variables
        real(real64), parameter :: s = 0.3d0
        real(real64), parameter :: a = -1.0d0
        real(real64), parameter :: b = 1.0d0
        real(real64), parameter :: g = 0.2d0
        real(real64), parameter :: w = 1.2d0

        ! Derivatives
        dxdt(1) = x(2)
        dxdt(2) = g * sin(w * t) - (s * x(2) + a * x(1) + b * x(1)**3)
    end subroutine
end program
