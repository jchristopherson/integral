! integral_core.f90

module integral_core
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use ferror
    implicit none

    interface
        !> @brief Defines a function of one variable to be integrated.
        !!
        !! @param[in] x The value of the independent variable at which to
        !!  evaluate the function.
        !! @return The value of the function at @p x.
        function integral_fcn(x) result(f)
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in) :: x
            real(real64) :: f
        end function
    end interface

    !> @brief Defines control parameters for the integration process.
    type integration_controls
        !> @brief The tolerance on absolute error.
        real(real64) :: absolute_tolerance
        !> @brief The tolerance on relative error.
        real(real64) :: relative_tolerance
        !> @brief The maximum number of subintervals the integrator is allowed.
        integer(int32) :: max_subinterval_count
    end type

    !> @brief Provides information regarding the behavior of an integrator.
    type integration_behavior
        !> @brief An estimate of the absolute error of the integration process.
        real(real64) :: error_estimate
        !> @brief The number of integrand evaluations.
        integer(int32) :: evaluation_count
        !> @brief The number of subintervals into which the integration region
        !! was divided.
        integer(int32) :: subinterval_count
    end type

! ******************************************************************************
! INTEGRAL_QUADPACK.F90
! ------------------------------------------------------------------------------
    interface
        module function qags(f, a, b, cntrls, bhvr, err) result(rst)
            procedure(integral_fcn), pointer, intent(in) :: f
            real(real64), intent(in) :: a, b
            type(integration_controls), intent(in), optional :: cntrls
            type(integration_behavior), intent(out), optional :: bhvr
            class(errors), intent(inout), optional, target :: err
            real(real64) :: rst
        end function
    end interface
end module
