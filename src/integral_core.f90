! integral_core.f90

module integral_core
    use, intrinsic :: iso_fortran_env, only : int32, real64
    use ferror
    implicit none
    private
    public :: INT_OUT_OF_MEMORY_ERROR
    public :: INT_COUNT_EXCEEDED_ERROR
    public :: INT_ROUND_OFF_ERROR
    public :: INT_INTEGRAND_BEHAVIOR_ERROR
    public :: INT_CONVERGENCE_ERROR
    public :: INT_DIVERGENT_ERROR
    public :: INT_INVALID_INPUT_ERROR
    public :: integrand
    public :: integration_controls
    public :: integration_behavior
    public :: integrate_adaptive

! ------------------------------------------------------------------------------
    !> @brief An error flag indicating insufficient memory.
    integer(int32), parameter :: INT_OUT_OF_MEMORY_ERROR = 1
    !> @brief An error flag indicates a specific parameter count has been
    !!  exceeded.
    integer(int32), parameter :: INT_COUNT_EXCEEDED_ERROR = 2
    !> @brief An error flag indicating round off error has become an issue.
    integer(int32), parameter :: INT_ROUND_OFF_ERROR = 3
    !> @brief An error flag indicating difficult integrand behavior.
    integer(int32), parameter :: INT_INTEGRAND_BEHAVIOR_ERROR = 4
    !> @brief An error flag indicating convergence issues.
    integer(int32), parameter :: INT_CONVERGENCE_ERROR = 5
    !> @brief An error flag indicating divergent behavior.
    integer(int32), parameter :: INT_DIVERGENT_ERROR = 6
    !> @brief An error flag indicating an invalid input.
    integer(int32), parameter :: INT_INVALID_INPUT_ERROR = 7

! ------------------------------------------------------------------------------
    interface
        !> @brief Defines a function of one variable to be integrated.
        !!
        !! @param[in] x The value of the independent variable at which to
        !!  evaluate the function.
        !! @return The value of the function at @p x.
        function integrand(x) result(f)
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in) :: x
            real(real64) :: f
        end function
    end interface

! ------------------------------------------------------------------------------
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
        !> @brief Utilizes an adaptive 21-point Gauss-Kronrod integrator to
        !! evaluate the integral of a general integrand over a finite interval.
        !!
        !! @param[in] f A pointer to the routine containing the integrand.
        !! @param[in] a The lower limit of integration.
        !! @param[in] b The upper limit of integration.
        !! @param[in] cntrls An optional parameter providing controls over the
        !!  the integration.  If nothing is specified, default controls are
        !!  utilized.
        !! @param[out] bhvr An optional output the can be used to gain
        !!  additional information about the integration process.
        !! @param[in,out] err An optional output that can be used to provide
        !!  an error handling mechanism.  If not provided, a default error
        !!  handling mechanism will be utilized.  Possible errors that may
        !!  be encountered are as follows.
        !!  - INT_OUT_OF_MEMORY_ERROR: There is insufficient memory available to
        !!      complete this operation.
        !!
        !! @par Remarks
        !! This routine utilizes the QUADPACK routine QAGS.  For more
        !! information on this routine see http://www.netlib.org/quadpack/.
        module function integrate_adaptive(f, a, b, cntrls, bhvr, err) result(rst)
            procedure(integrand), pointer, intent(in) :: f
            real(real64), intent(in) :: a, b
            type(integration_controls), intent(in), optional :: cntrls
            type(integration_behavior), intent(out), optional :: bhvr
            class(errors), intent(inout), optional, target :: err
            real(real64) :: rst
        end function
    end interface
end module
