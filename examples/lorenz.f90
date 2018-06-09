! lorenz.f90

program example
    use iso_fortran_env
    use integral_core
    use fplot_core
    implicit none

    ! Local Variables
    type(ode_helper) :: fcn
    type(ode_auto) :: integrator
    procedure(ode_fcn), pointer :: ptr
    real(real64) :: ic(3), t(2)
    real(real64), allocatable, dimension(:,:) :: x
    type(plot_3d) :: plt
    type(plot_data_3d) :: d1
    class(plot_axis), pointer :: xAxis, yAxis, zAxis

    ! Set up the integrator
    ptr => lorenz
    call fcn%define_equations(3, ptr)
    ic = [1.0d0, 1.0d0, 1.0d0]
    t = [0.0d0, 1.0d2]

    ! Integrate
    x = integrator%integrate(fcn, t, ic)

    ! Plot
    call plt%initialize()
    call plt%set_font_size(14)

    xAxis => plt%get_x_axis()
    call xAxis%set_title("x(t)")

    yAxis => plt%get_y_axis()
    call yAxis%set_title("y(t)")

    zAxis => plt%get_z_axis()
    call zAxis%set_title("z(t)")

    call d1%set_line_color(CLR_BLUE)
    call d1%define_data(x(:,2), x(:,3), x(:,4))

    call plt%push(d1)
    call plt%draw()

contains
    ! The Lorenz system of equations:
    ! REF: https://en.wikipedia.org/wiki/Lorenz_system
    !
    ! x' = s * (y - x)
    ! y' = x * (r - z) - y
    ! z' = x * y - b * z
    subroutine lorenz(t, x, dxdt)
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: dxdt

        ! Parameters
        real(real64), parameter :: r = 28.0d0
        real(real64), parameter :: s = 10.0d0
        real(real64), parameter :: b = 8.0d0 / 3.0d0

        ! Equations
        dxdt(1) = s * (x(2) - x(1))
        dxdt(2) = x(1) * (r - x(3)) - x(2)
        dxdt(3) = x(1) * x(2) - b * x(3)
    end subroutine
end program
