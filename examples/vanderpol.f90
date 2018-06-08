! vanderpol.f90

program example
    use iso_fortran_env
    use integral_core
    use fplot_core
    implicit none

    ! Local Variables
    type(ode_helper) :: fcn
    type(ode_auto) :: integrator
    procedure(ode_fcn), pointer :: ptr
    real(real64) :: ic(2), t(2)
    real(real64), allocatable, dimension(:,:) :: x
    type(plot_2d) :: plt
    type(plot_data_2d) :: d1
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd

    ! Set up the integrator
    ptr => vdp
    call fcn%define_equations(2, ptr)

    ! Define the initial conditions
    t = [0.0d0, 8.0d1]
    ic = [2.0d0, 0.0d0]

    ! Compute the solution
    x = integrator%integrate(fcn, t, ic)

    ! Plot the solution
    call plt%initialize()
    call plt%set_font_size(14)

    lgnd => plt%get_legend()
    call lgnd%set_is_visible(.false.)

    xAxis => plt%get_x_axis()
    call xAxis%set_title("t")

    yAxis => plt%get_y_axis()
    call yAxis%set_title("x(t)")

    call d1%set_name("x(t)")
    call d1%set_line_width(2.0)
    call d1%set_line_color(CLR_BLUE)
    call d1%define_data(x(:,1), x(:,2))

    call plt%push(d1)
    call plt%draw()

contains
    ! Van Der Pol Equation
    subroutine vdp(t, x, dxdt)
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: dxdt

        real(real64), parameter :: mu = 20.0d0

        dxdt(1) = x(2)
        dxdt(2) = mu * (1.0d0 - x(1)**2) * x(2) - x(1)
    end subroutine
end program
