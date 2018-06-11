! integral_core.f90

!> @mainpage
!!
!! @section intro_sec Introduction
!! INTEGRAL is a Fortran library providing an easy-to-use, object-oriented
!! interface to several integration routines.  Integration of functions of one
!! variable are provided by the QUADPACK library.  Integration of systems of
!! ordinary differential equations are provided by the ODEPACK library.
!!
!! @image html double_pendulum_example_short_time_no_legend.png

!> @brief \b integral_core
!!
!! @par Purpose
!! Provides types and routines allowing for the integration of functions and
!! systems of ordinary differential equations.
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
    public :: INT_LACK_OF_DEFINITION_ERROR
    public :: INT_ARRAY_SIZE_MISMATCH_ERROR
    public :: INT_EXCESSIVE_WORK_ERROR
    public :: INT_IMPRACTICAL_TOLERANCE_ERROR
    public :: INT_REPEATED_ERROR_TEST_FAILURE
    public :: INT_15_POINT_RULE
    public :: INT_21_POINT_RULE
    public :: INT_31_POINT_RULE
    public :: INT_41_POINT_RULE
    public :: INT_51_POINT_RULE
    public :: INT_61_POINT_RULE
    public :: integrand
    public :: ode_fcn
    public :: ode_jacobian
    public :: ode_constraint
    public :: integration_behavior
    public :: integrator_base
    public :: finite_interval_integrator
    public :: finite_interval_fcn
    public :: adaptive_integrator
    public :: nonadaptive_integrator
    public :: ode_helper
    public :: ode_integrator
    public :: ode_integrator_interface
    public :: ode_integrator_reset
    public :: ode_auto

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
    !> @brief An error flag denoting a lack of definition.
    integer(int32), parameter :: INT_LACK_OF_DEFINITION_ERROR = 8
    !> @brief An error indicating an inappropriately sized array.
    integer(int32), parameter :: INT_ARRAY_SIZE_MISMATCH_ERROR = 9
    !> @brief An error indicating the too many iterations have been performed
    !! suggesting that excessive amounts of work are required to continue the
    !! solution process.
    integer(int32), parameter :: INT_EXCESSIVE_WORK_ERROR = 10
    !> @brief An error indicating that the user-defined tolerances are too
    !! stringent to be practical for the problem at hand.
    integer(int32), parameter :: INT_IMPRACTICAL_TOLERANCE_ERROR = 11
    !> @brief An error that occurs if integrator error tests fail repeatadly.
    integer(int32), parameter :: INT_REPEATED_ERROR_TEST_FAILURE = 12

! ------------------------------------------------------------------------------
    !> @brief Defines a 15-point Gauss-Kronrod integration rule.
    integer(int32), parameter :: INT_15_POINT_RULE = 15
    !> @brief Defines a 21-point Gauss-Kronrod integration rule.
    integer(int32), parameter :: INT_21_POINT_RULE = 21
    !> @brief Defines a 31-point Gauss-Kronrod integration rule.
    integer(int32), parameter :: INT_31_POINT_RULE = 31
    !> @brief Defines a 41-point Gauss-Kronrod integration rule.
    integer(int32), parameter :: INT_41_POINT_RULE = 41
    !> @brief Defines a 51-point Gauss-Kronrod integration rule.
    integer(int32), parameter :: INT_51_POINT_RULE = 51
    !> @brief Defines a 61-point Gauss-Kronrod integration rule.
    integer(int32), parameter :: INT_61_POINT_RULE = 61

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

        !> @brief Defines a routine containing a system of first order
        !! ODEs.
        !!
        !! @param[in] x The value of the independent variable at which to
        !!  evalaute the ODEs.
        !! @param[in] y An N-element array cotnaining the current estimates
        !!  of the dependent variables as evaluated at @p x.
        !! @param[out] dydx An N-element array where the values of the ODEs
        !!  as evaluated at @p x are to be written.
        subroutine ode_fcn(x, y, dydx)
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: dydx
        end subroutine

        !> @brief Defines a routine capable of computing the Jacobian matrix
        !! of a system of N first order ODEs of N variables.
        !!
        !! @param[in] x The value of the independent variable at which the
        !!  Jacobian is to be evaluated.
        !! @param[in] y An N-element array containing the values of the
        !!  dependent variables at @p x.
        !! @param[out] jac An N-by-N matrix where the Jacobian is to be written.
        subroutine ode_jacobian(x, y, jac)
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:,:) :: jac
        end subroutine

        !> @brief Defines a routine capable of describing constraints on a
        !! system of ODEs of N variables exposed to M constraints.
        !!
        !! @param[in] x The value of the independent variable at which the
        !!  Jacobian is to be evaluated.
        !! @param[in] y An N-element array containing the values of the
        !!  dependent variables at @p x.
        !! @param[out] f An M-element array where the values of the M
        !!  constraint equations are to be written.
        subroutine ode_constraint(x, y, f)
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: f
        end subroutine
    end interface

! ------------------------------------------------------------------------------
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
! INTEGRAL_INTEGRATOR_BASE.F90
! ------------------------------------------------------------------------------
    !> @brief Defines a base type for integrator types.
    type integrator_base
    private
        !> @brief The absolute tolerance.
        real(real64) :: m_absTol = 1.0d-8
        !> @brief The relative tolerance.
        real(real64) :: m_relTol = 1.0d-8
        !> @brief The maximum number of subintervals.
        integer(int32) :: m_maxInt = 10
    contains
        !> @brief Gets the absolute tolerance value.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure real(real64) function get_abs_tol(class(integrator_base) this)
        !! @endcode
        !!
        !! @param[in] this The integrator_base object.
        !! @return The tolerance value.
        procedure, public :: get_abs_tol => ib_get_abs_tol
        !> @brief Sets the absolute tolerance value.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_abs_tol(class(integrator_base) this, real(real64) x)
        !! @endcode
        !!
        !! @param[in,out] this The integrator_base object.
        !! @param[in] The tolerance value.
        procedure, public :: set_abs_tol => ib_set_abs_tol
        !> @brief Gets the relative tolerance value.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure real(real64) function get_rel_tol(class(integrator_base) this)
        !! @endcode
        !!
        !! @param[in] this The integrator_base object.
        !! @return The tolerance value.
        procedure, public :: get_rel_tol => ib_get_rel_tol
        !> @brief Sets the relative tolerance value.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_rel_tol(class(integrator_base) this, real(real64) x)
        !! @endcode
        !!
        !! @param[in,out] this The integrator_base object.
        !! @param[in] The tolerance value.
        procedure, public :: set_rel_tol => ib_set_rel_tol
        !> @brief Gets the maximum number of subintervals into which the
        !! integrator may divide the problem region.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) function get_max_subintervals(class(integrator_base) this)
        !! @endcode
        !!
        !! @param[in] this The integrator_base object.
        !! @return The number of intervals.
        procedure, public :: get_max_subintervals => ib_get_max_subintervals
        !> @brief Sets the maximum number of subintervals into which the
        !! integrator may divide the problem region.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_max_subintervals(class(integrator_base) this, integer(int32) x)
        !! @endcode
        !!
        !! @param[in,out] this The integrator_base object.
        !! @param[in] x The number of intervals.
        procedure, public :: ib_set_max_subintervals => ib_set_max_subintervals
    end type

! ------------------------------------------------------------------------------
    interface
        pure module function ib_get_abs_tol(this) result(x)
            class(integrator_base), intent(in) :: this
            real(real64) :: x
        end function

        module subroutine ib_set_abs_tol(this, x)
            class(integrator_base), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        pure module function ib_get_rel_tol(this) result(x)
            class(integrator_base), intent(in) :: this
            real(real64) :: x
        end function

        module subroutine ib_set_rel_tol(this, x)
            class(integrator_base), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        pure module function ib_get_max_subintervals(this) result(x)
            class(integrator_base), intent(in) :: this
            integer(int32) :: x
        end function

        module subroutine ib_set_max_subintervals(this, x)
            class(integrator_base), intent(inout) :: this
            integer(int32), intent(in) :: x
        end subroutine
    end interface

! ******************************************************************************
    !> @brief A type that defines an integrator meant to operate on integrands
    !! over a finite region.
    type, abstract, extends(integrator_base) :: finite_interval_integrator
    contains
        !> @brief Performs the actual integration.
        procedure(finite_interval_fcn), public, deferred, pass :: integrate
    end type

    interface
        !> @brief Defines the signature of a routine used to integrate
        !! a function of one variable over a finite interval.
        !!
        !! @param[in,out] this The finite_interval_integrator object.
        !! @param[in] fcn The integrand.
        !! @param[in] a The lower limit of integration.
        !! @param[in] b The upper limit of integration.
        !! @param[out] info An optional output providing information regarding
        !!  behavior of the integrator.
        !! @param[out] err An optional argument that may be used to provide
        !!  customized error handling behaviors.
        !! @return The value of the integral over the specified range.
        function finite_interval_fcn(this, fcn, a, b, info, err) result(rst)
            use, intrinsic :: iso_fortran_env, only : real64
            use ferror
            import finite_interval_integrator
            import integrand
            import integration_behavior
            class(finite_interval_integrator), intent(inout) :: this
            procedure(integrand), intent(in), pointer :: fcn
            real(real64), intent(in) :: a, b
            type(integration_behavior), intent(out), optional :: info
            class(errors), intent(inout), optional, target :: err
            real(real64) :: rst
        end function
    end interface

! ******************************************************************************
! INTEGRAL_ADAPTIVE_INTEGRATOR.F90
! ------------------------------------------------------------------------------
    !> @brief Defines an integrator that uses an adaptive Gauss-Kronrod method
    !! to compute the integral of a function of one variable over a finite
    !! interval.
    !!
    !! @par Example
    !! The following example illustrates how to use the adaptive_integrator to
    !! compute the integral of a function over a finite interval.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use integral_core
    !!     implicit none
    !!
    !!     ! Variables
    !!     real(real64) :: ans, y, pi, a, b
    !!     procedure(integrand), pointer :: fcn
    !!     type(adaptive_integrator) :: integrator
    !!
    !!     ! Define the integration limits
    !!     pi = 2.0d0 * acos(0.0d0)
    !!     a = pi / 6.0d0
    !!     b = pi / 4.0d0
    !!
    !!     ! Evaluate the integral
    !!     fcn => int_fcn
    !!     y = integrator%integrate(fcn, a, b)
    !!
    !!     ! Display the results
    !!     ans = 5.0d0 * pi / 12.0d0 - 2.0d0 * sqrt(2.0d0) + 4.0d0 / sqrt(3.0d0)
    !!     print '(AEN13.5AEN13.5A)', "The solution is: ", ans, &
    !!         ", the integrator computed: ", y, "."
    !!
    !! contains
    !!     ! This example is from http://tutorial.math.lamar.edu/Classes/CalcI/ComputingDefiniteIntegrals.aspx#Int_CompDef_Ex3a
    !!     ! The integrand is: f(x) = 5 - 2 sec(x) tan(x).
    !!     ! If the integral is considered over the range [pi/6, pi/4], the solution
    !!     ! is 5 pi / 12 - 2 sqrt(2) + 4 / sqrt(3).
    !!     function int_fcn(x) result(f)
    !!         real(real64), intent(in) :: x
    !!         real(real64) :: f
    !!         f = 5.0d0 - 2.0d0 * tan(x) / cos(x) ! Remember, sec(x) = 1 / cos(x)
    !!     end function
    !! end program
    !! @endcode
    !! The above program produces the following output.
    !! @code{.txt}
    !! The solution is: 789.97089E-03, the integrator computed: 789.97089E-03.
    !! @endcode
    type, extends(finite_interval_integrator) :: adaptive_integrator
    private
        !> @brief A workspace array.
        real(real64), allocatable, dimension(:) :: m_work
        !> @brief A workspace array.
        integer(int32), allocatable, dimension(:) :: m_iwork
        !> @brief True for user defined breakpoints; else, false.
        logical :: m_userDefinedBreaks = .false.
        !> @brief A list of user-defined breakpoints.
        real(real64), allocatable, dimension(:) :: m_breakpoints
    contains
        !> @brief Initializes the adaptive_integrator object.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine initialize(class(adaptive_integrator) this, optional class(errors) err)
        !! @endcode
        !!
        !! @param[in] this The adaptive_integrator object.
        !! @param[in,out] err An optional output that can be used to provide
        !!  an error handling mechanism.  If not provided, a default error
        !!  handling mechanism will be utilized.  Possible errors that may
        !!  be encountered are as follows.
        !!  - INT_OUT_OF_MEMORY_ERROR: There is insufficient memory available to
        !!      complete this operation.
        procedure, public :: initialize => ai_init
        !> @brief Performs the actual integration.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function integrate(class(adaptive_integrator) this, procedure(integrand) pointer fcn, real(real64) a, real(real64) b, optional type(integration_behavior) info, optional class(errors) err)
        !! @endcode
        !!
        !! @param[in,out] this The adaptive_integrator object.
        !! @param[in] fcn The integrand.
        !! @param[in] a The lower limit of integration.
        !! @param[in] b The upper limit of integration.
        !! @param[out] info An optional output providing information regarding
        !!  behavior of the integrator.
        !! @param[in,out] err An optional output that can be used to provide
        !!  an error handling mechanism.  If not provided, a default error
        !!  handling mechanism will be utilized.  Possible errors that may
        !!  be encountered are as follows.
        !!  - INT_OUT_OF_MEMORY_ERROR: There is insufficient memory available to
        !!      complete this operation.
        !!  - INT_COUNT_EXCEEDED_ERROR: The maximum number of subdivisions has
        !!      been reached.
        !!  - INT_ROUND_OFF_ERROR: The occurence of roundoff error is preventing
        !!      the integrator from reaching the requested tolerances.
        !!  - INT_INTEGRAND_BEHAVIOR_ERROR: The integrand appears to be
        !!      behaving too poorly to proceed.
        !!  - INT_CONVERGENCE_ERROR: The algorithm could not converge with the
        !!      given constraints.
        !!  - INT_DIVERGENT_ERROR: The integral is likely divergent.
        !!  - INT_INVALID_INPUT_ERROR: An invalid input was supplied.
        !!
        !! @return The value of the integral over the specified range.
        !!
        !! @par Remarks
        !! This routine utilizes the QUADPACK routine QAGS.  For more
        !! information on this routine see http://www.netlib.org/quadpack/.
        procedure, public :: integrate => ai_integrate
        !> @brief Gets a flag determining if user-defined breakpoints should
        !! be used to inform the integrator of points of difficulty within
        !! the integration region.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical function get_use_breakpoints(class(adaptive_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The adaptive_integrator object.
        !! @return Returns true if the user-defined breakpoints should be
        !!  used; else, false.
        procedure, public :: get_use_breakpoints => ai_get_use_brkpnts
        !> @brief Sets a flag determining if user-defined breakpoints should
        !! be used to inform the integrator of points of difficulty within
        !! the integration region.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_use_breakpoints(class(adaptive_integrator) this, logical x)
        !! @endcode
        !!
        !! @param[in,out] this The adaptive_integrator object.
        !! @param[in] x Set to true if the user-defined breakpoints should be
        !!  used; else, set to false.
        procedure, public :: set_use_breakpoints => ai_set_use_brkpnts
        !> @brief Gets an array of user-defined breakpoints.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64)(:) allocatable function get_breakpoints(class(adaptive_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The adaptive_integrator object.
        !! @return The array of points.  This array will not be allocated if
        !!  no points have been defined.
        procedure, public :: get_breakpoints => ai_get_breakpoints
        !> @brief Sets an array of user-defined breakpoints.  Additionally,
        !! this array call set_use_breakpoints with an argument of true such
        !! that the integrator will use the breakpoints defined by this routine.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_breakpoints(class(adaptive_integrator) this, real(real64) x(:))
        !! @endcode
        !!
        !! @param[in,out] this The adaptive_integrator object.
        !! @param[in] x The array of breakpoints.  Notice, each point must lie
        !!  within the integration interval.
        procedure, public :: set_breakpoints => ai_set_breakpoints
    end type

! ------------------------------------------------------------------------------
    interface
        module subroutine ai_init(this, err)
            class(adaptive_integrator), intent(inout) :: this
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module function ai_integrate(this, fcn, a, b, info, err) result(rst)
            class(adaptive_integrator), intent(inout) :: this
            procedure(integrand), intent(in), pointer :: fcn
            real(real64), intent(in) :: a, b
            type(integration_behavior), intent(out), optional :: info
            class(errors), intent(inout), optional, target :: err
            real(real64) :: rst
        end function

        pure module function ai_get_use_brkpnts(this) result(x)
            class(adaptive_integrator), intent(in) :: this
            logical :: x
        end function

        module subroutine ai_set_use_brkpnts(this, x)
            class(adaptive_integrator), intent(inout) :: this
            logical, intent(in) :: x
        end subroutine

        module function ai_get_breakpoints(this) result(x)
            class(adaptive_integrator), intent(in) :: this
            real(real64), allocatable, dimension(:) :: x
        end function

        module subroutine ai_set_breakpoints(this, x)
            class(adaptive_integrator), intent(inout) :: this
            real(real64), intent(in), dimension(:) :: x
        end subroutine
    end interface

! ******************************************************************************
! INTEGRAL_NONADAPTIVE_INTEGRATOR.F90
! ------------------------------------------------------------------------------
    !> @brief Defines an integrator that uses a non-adaptive Gauss-Kronrod
    !! method to compute the integral of a function of one variable over a
    !! finite interval.
    type, extends(finite_interval_integrator) :: nonadaptive_integrator
    private
        !> @brief The integration rule to use.
        integer(int32) :: m_rule = INT_15_POINT_RULE
    contains
        !> @brief Performs the actual integration.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) function integrate(class(nonadaptive_integrator) this, procedure(integrand) pointer fcn, real(real64) a, real(real64) b, optional type(integration_behavior) info, optional class(errors) err)
        !! @endcode
        !!
        !! @param[in,out] this The nonadaptive_integrator object.
        !! @param[in] fcn The integrand.
        !! @param[in] a The lower limit of integration.
        !! @param[in] b The upper limit of integration.
        !! @param[out] info An optional output providing information regarding
        !!  behavior of the integrator.
        !! @param[in,out] err An optional output that can be used to provide
        !!  an error handling mechanism.  If not provided, a default error
        !!  handling mechanism will be utilized.  Possible errors that may
        !!  be encountered are as follows.
        !!  - INT_INVALID_INPUT_ERROR: An invalid input was supplied.
        !!
        !! @return The value of the integral over the specified range.
        !!
        !! @par Remarks
        !! This routine utilizes the QUADPACK routine QAGS.  For more
        !! information on this routine see http://www.netlib.org/quadpack/.
        procedure, public :: integrate => ni_integrate
        !> @brief Gets the integration rule being utilized.  The default
        !! integration rule is given by INT_15_POINT_RULE.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) function get_rule(class(nonadaptive_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The nonadaptive_integrator object.
        !! @return The integration rule identifier.  The available rules are as
        !!  follows.
        !!  - INT_15_POINT_RULE
        !!  - INT_21_POINT_RULE
        !!  - INT_31_POINT_RULE
        !!  - INT_41_POINT_RULE
        !!  - INT_51_POINT_RULE
        !!  - INT_61_POINT_RULE
        procedure, public :: get_rule => ni_get_rule
        !> @brief Sets the integration rule being utilized.  The default
        !! integration rule is given by INT_15_POINT_RULE.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_rule(class(nonadaptive_integrator) this, integer(int32) x, optional class(errors) err)
        !! @endcode
        !!
        !! @param[in,out] this The nonadaptive_integrator object.
        !! @param[in] x The integration rule identifier.  The available rules
        !!  are as follows.
        !!  - INT_15_POINT_RULE
        !!  - INT_21_POINT_RULE
        !!  - INT_31_POINT_RULE
        !!  - INT_41_POINT_RULE
        !!  - INT_51_POINT_RULE
        !!  - INT_61_POINT_RULE
        !! @param[in,out] err An optional output that can be used to provide
        !!  an error handling mechanism.  If not provided, a default error
        !!  handling mechanism will be utilized.  Possible errors that may
        !!  be encountered are as follows.
        !!  - INT_INVALID_INPUT_ERROR: An invalid input was supplied.
        procedure, public :: set_rule => ni_set_rule
    end type

! ------------------------------------------------------------------------------
    interface
        module function ni_integrate(this, fcn, a, b, info, err) result(rst)
            class(nonadaptive_integrator), intent(inout) :: this
            procedure(integrand), intent(in), pointer :: fcn
            real(real64), intent(in) :: a, b
            type(integration_behavior), intent(out), optional :: info
            class(errors), intent(inout), optional, target :: err
            real(real64) :: rst
        end function

        pure module function ni_get_rule(this) result(x)
            class(nonadaptive_integrator), intent(in) :: this
            integer(int32) :: x
        end function

        module subroutine ni_set_rule(this, x, err)
            class(nonadaptive_integrator), intent(inout) :: this
            integer(int32), intent(in) :: x
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

! ******************************************************************************
! INTEGRAL_ODE_HELPER.F90
! ------------------------------------------------------------------------------
    !> @brief Defines a type used to pass information regarding the ODEs to
    !! the solver.
    type :: ode_helper
    private
        !> @brief A pointer to the routine containing the ODEs to integrate.
        procedure(ode_fcn), pointer, nopass :: m_fcn => null()
        !> @brief A pointer to the routine containing the Jacobian.
        procedure(ode_jacobian), pointer, nopass :: m_jac => null()
        !> @brief A pointer to a routine containing any constraint equations.
        procedure(ode_constraint), pointer, nopass :: m_rts => null()
        !> @brief The number of first order ODEs to integrate.
        integer(int32) :: m_count = 0
        !> @brief The number of constraint equations.
        integer(int32) :: m_constraints = 0
    contains
        !> @brief Defines the system of ODEs to solve.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine define_equations(class(ode_helper) this, integer(int32) neqn, procedure(ode_fcn) pointer fcn, optional procedure(ode_jacobian) pointer jac)
        !! @endcode
        !!
        !! @param[in,out] this The ode_helper object.
        !! @param[in] neqn The number of first order ODEs.
        !! @param[in] fcn A pointer to the routine containing the ODEs.
        !! @param[in] jac An optional pointer to the routine used to compute
        !!  the Jacobian of the system of equations given in @p fcn.  If the
        !!  solver requires a Jacobian, and no routine is provided, the solver
        !!  will generate a numerical estimate of the Jacobian matrix to use
        !!  in its place.
        !!
        !! @par Example
        !! The following example illustrates how to define and solve a system of
        !! differential equations.  The Lorenz system is illustrated.
        !! @code{.f90}
        !! program example
        !!     use iso_fortran_env
        !!     use integral_core
        !!     use fplot_core
        !!     implicit none
        !!
        !!     ! Local Variables
        !!     type(ode_helper) :: fcn
        !!     type(ode_auto) :: integrator
        !!     procedure(ode_fcn), pointer :: ptr
        !!     real(real64) :: ic(3), t(2)
        !!     real(real64), allocatable, dimension(:,:) :: x
        !!     type(plot_3d) :: plt
        !!     type(plot_data_3d) :: d1
        !!     class(plot_axis), pointer :: xAxis, yAxis, zAxis
        !!
        !!     ! Set up the integrator
        !!     ptr => lorenz
        !!     call fcn%define_equations(3, ptr)
        !!     ic = [1.0d0, 1.0d0, 1.0d0]
        !!     t = [0.0d0, 1.0d2]
        !!
        !!     ! Integrate
        !!     x = integrator%integrate(fcn, t, ic)
        !!
        !!     ! Plot
        !!     call plt%initialize()
        !!     call plt%set_font_size(14)
        !!
        !!     xAxis => plt%get_x_axis()
        !!     call xAxis%set_title("x(t)")
        !!
        !!     yAxis => plt%get_y_axis()
        !!     call yAxis%set_title("y(t)")
        !!
        !!     zAxis => plt%get_z_axis()
        !!     call zAxis%set_title("z(t)")
        !!
        !!     call d1%set_line_color(CLR_BLUE)
        !!     call d1%define_data(x(:,2), x(:,3), x(:,4))
        !!
        !!     call plt%push(d1)
        !!     call plt%draw()
        !!
        !! contains
        !!     ! The Lorenz system of equations:
        !!     ! REF: https://en.wikipedia.org/wiki/Lorenz_system
        !!     !
        !!     ! x' = s * (y - x)
        !!     ! y' = x * (r - z) - y
        !!     ! z' = x * y - b * z
        !!     subroutine lorenz(t, x, dxdt)
        !!         real(real64), intent(in) :: t
        !!         real(real64), intent(in), dimension(:) :: x
        !!         real(real64), intent(out), dimension(:) :: dxdt
        !!
        !!         ! Parameters
        !!         real(real64), parameter :: r = 28.0d0
        !!         real(real64), parameter :: s = 10.0d0
        !!         real(real64), parameter :: b = 8.0d0 / 3.0d0
        !!
        !!         ! Equations
        !!         dxdt(1) = s * (x(2) - x(1))
        !!         dxdt(2) = x(1) * (r - x(3)) - x(2)
        !!         dxdt(3) = x(1) * x(2) - b * x(3)
        !!     end subroutine
        !! end program
        !! @endcode
        !! The above program produces the following output.
        !! @image html lorenz_example.png
        procedure, public :: define_equations => oh_init
        !> @brief Gets the number of ODEs.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) function get_equation_count(class(ode_helper) this)
        !! @endcode
        !!
        !! @param[in] this The ode_helper object.
        !! @return The number of ODEs.
        procedure, public :: get_equation_count => oh_get_count
        !> @brief Gets a flag denoting if the user has defined the system
        !! of equations.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical function get_equations_defined(class(ode_helper) this)
        !! @endcode
        !!
        !! @param[in] this The ode_helper object.
        !! @return Returns true if the equations have been defined (by calling
        !!  @p define_equations); else, false if they have not yet been
        !!  defined.
        procedure, public :: get_equations_defined => oh_is_fcn_defined
        !> @brief Gets a flag denoting if the user has defined the routine
        !! used to compute the Jacobian matrix.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical function get_jacobian_defined(class(ode_helper) this)
        !! @endcode
        !!
        !! @param[in] this The ode_helper object.
        !! @return Returns true if the routine for computing the Jacobian
        !!  has been defined (by calling @p define_equations); else, false if
        !!  it has not.
        procedure, public :: get_jacobian_defined => oh_is_jac_defined
        !> @brief Defines a routine for applying constraints to the ODEs.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine define_constraints(class(ode_helper) this, integer(int32) n, procedure(ode_constraint) pointer fcn)
        !! @endcode
        !!
        !! @param[in,out] this The ode_helper object.
        !! @param[in] n The number of constraint equations.
        !! @param[in] fcn A pointer to the routine containing the constraint
        !!  equations.
        !!
        !! @par Example
        !! The following example illustrates how to employ constraints by
        !! modeling a bouncing ball.
        !! @code{.f90}
        !! program example
        !!     use iso_fortran_env
        !!     use integral_core
        !!     use fplot_core
        !!     implicit none
        !!
        !!     ! Parameters
        !!     real(real64), parameter :: g = 9.81d0 ! Gravitational acceleration
        !!     real(real64), parameter :: k = -0.8d0 ! Coefficient of restitution
        !!
        !!     ! Local Variables
        !!     procedure(ode_fcn), pointer :: ptr
        !!     procedure(ode_constraint), pointer :: cptr
        !!     type(ode_helper) :: fcn
        !!     type(ode_auto) :: integrator
        !!     integer(int32) :: n
        !!     real(real64) :: ic(2), t(2)
        !!     real(real64), allocatable, dimension(:,:) :: x1, x2, x3, x4
        !!     type(plot_2d) :: plt
        !!     type(plot_data_2d) :: d1, d2, d3, d4
        !!     class(plot_axis), pointer :: xAxis, yAxis
        !!     type(legend), pointer :: lgnd
        !!
        !!     ! Set up the integrator
        !!     ptr => ball
        !!     cptr => ground_constraint
        !!     call fcn%define_equations(2, ptr)
        !!     call fcn%define_constraints(1, cptr)
        !!     call integrator%set_max_step_size(1.0d-3)
        !!     call integrator%set_limit_step_size(.true.)
        !!
        !!     ! Compute the solution
        !!     t = [0.0d0, 1.0d1]
        !!     ic = [1.0d1, 5.0d0]
        !!     x1 = integrator%integrate(fcn, t, ic)
        !!
        !!     ! The integrator stops when the ball first makes contact.  As a result, lets
        !!     ! reset the time limits and initial conditions to continue the integration
        !!     n = size(x1, 1)
        !!     t(1) = x1(n,1)
        !!     ic = [abs(x1(n,2)), k * x1(n,3)]
        !!     call integrator%reset()
        !!     x2 = integrator%integrate(fcn, t, ic)
        !!
        !!     ! Again
        !!     n = size(x2, 1)
        !!     t(1) = x2(n,1)
        !!     ic = [abs(x2(n,2)), k * x2(n,3)]
        !!     call integrator%reset()
        !!     x3 = integrator%integrate(fcn, t, ic)
        !!
        !!     ! Again
        !!     n = size(x3, 1)
        !!     t(1) = x3(n,1)
        !!     ic = [abs(x3(n,2)), k * x3(n,3)]
        !!     call integrator%reset()
        !!     x4 = integrator%integrate(fcn, t, ic)
        !!
        !!
        !!     ! Plot the solution
        !!     call plt%initialize()
        !!     call plt%set_font_size(14)
        !!
        !!     lgnd => plt%get_legend()
        !!     call lgnd%set_is_visible(.false.)
        !!
        !!     xAxis => plt%get_x_axis()
        !!     call xAxis%set_title("t")
        !!
        !!     yAxis => plt%get_y_axis()
        !!     call yAxis%set_title("x(t)")
        !!
        !!     call d1%set_line_color(CLR_BLUE)
        !!     call d1%set_line_width(2.0)
        !!     call d1%define_data(x1(:,1), x1(:,2))
        !!
        !!     call d2%set_line_color(CLR_BLUE)
        !!     call d2%set_line_width(2.0)
        !!     call d2%define_data(x2(:,1), x2(:,2))
        !!
        !!     call d3%set_line_color(CLR_BLUE)
        !!     call d3%set_line_width(2.0)
        !!     call d3%define_data(x3(:,1), x3(:,2))
        !!
        !!     call d4%set_line_color(CLR_BLUE)
        !!     call d4%set_line_width(2.0)
        !!     call d4%define_data(x4(:,1), x4(:,2))
        !!
        !!     call plt%push(d1)
        !!     call plt%push(d2)
        !!     call plt%push(d3)
        !!     call plt%push(d4)
        !!     call plt%draw()
        !!
        !! contains
        !!     ! A bouncing ball can be described by the following equation:
        !!     ! x" = -g
        !!     !
        !!     ! Where g = gravitational acceleration
        !!     subroutine ball(t, x, dxdt)
        !!         real(real64), intent(in) :: t
        !!         real(real64), intent(in), dimension(:) :: x
        !!         real(real64), intent(out), dimension(:) :: dxdt
        !!         dxdt(1) = x(2)
        !!         dxdt(2) = -g
        !!     end subroutine
        !!
        !!     ! The constraint function
        !!     subroutine ground_constraint(t, x, f)
        !!         real(real64), intent(in) :: t
        !!         real(real64), intent(in), dimension(:) :: x
        !!         real(real64), intent(out), dimension(:) :: f
        !!         f(1) = x(1)     ! Find when x == 0
        !!     end subroutine
        !! end program
        !! @endcode
        !! The above program produces the following output.
        !! @image html bouncing_ball_example.png
        procedure, public :: define_constraints => oh_define_constraints
        !> @brief Gets the number of constraint equations.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) function get_constraint_count(class(ode_helper) this)
        !! @endcode
        !!
        !! @param[in] this The ode_helper object.
        !! @return The number of constraint equations.
        procedure, public :: get_constraint_count => oh_get_constraint_count
        !> @brief Determines if any constraints have been defined.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical function get_constraints_defined(class(ode_helper) this)
        !! @endcode
        !!
        !! @param[in] this The ode_helper object.
        !! @return True if the constraint equation routine is defined; else,
        !!  false.
        procedure, public :: get_constraints_defined => oh_get_constraints_defined
        !! @brief Evaluates the routine containing the system of ODEs.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine evaluate_ode(class(ode_helper) this, real(real64) x, real(real64) y(:), real(real64) dydx(:))
        !! @endcode
        !!
        !! @param[in,out] this The ode_helper object.
        !! @param[in] x The value of the independent variable at which to
        !!  evalaute the ODEs.
        !! @param[in] y An N-element array cotnaining the current estimates
        !!  of the dependent variables as evaluated at @p x.
        !! @param[out] dydx An N-element array where the values of the ODEs
        !!  as evaluated at @p x are to be written.
        procedure, public :: evaluate_ode => oh_eval
        !> @brief Evaluates the Jacobian matrix.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine evalaute_jacobian(class(ode_helper) this, real(real64) x, real(real64) y(:), real(real64) jac(:,:))
        !! @endcode
        !!
        !! @param[in,out] this The ode_helper object.
        !! @param[in] x The value of the independent variable at which the
        !!  Jacobian is to be evaluated.
        !! @param[in] y An N-element array containing the values of the
        !!  dependent variables at @p x.
        !! @param[out] jac An N-by-N matrix where the Jacobian is to be written.
        procedure, public :: evaluate_jacobian => oh_eval_jac
        !> @brief Evaluates the constraint equations.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine evaluate_constraints(class(ode_helper) this, real(real64) x, real(real64) y(:), real(real64) f(:))
        !! @endcode
        !!
        !! @param[in,out] this The ode_helper object.
        !! @param[in] x The value of the independent variable at which the
        !!  Jacobian is to be evaluated.
        !! @param[in] y An N-element array containing the values of the
        !!  dependent variables at @p x where N is the number of ODEs.
        !! @param[out] f An M-element array where the values of the M
        !!  constraint equations are to be written.
        procedure, public :: evaluate_constraints => oh_eval_constraints
    end type

    interface
        module subroutine oh_init(this, neqn, fcn, jac)
            class(ode_helper), intent(inout) :: this
            integer(int32), intent(in) :: neqn
            procedure(ode_fcn), pointer, intent(in) :: fcn
            procedure(ode_jacobian), pointer, intent(in), optional :: jac
        end subroutine

        pure module function oh_get_count(this) result(x)
            class(ode_helper), intent(in) :: this
            integer(int32) :: x
        end function

        pure module function oh_is_fcn_defined(this) result(x)
            class(ode_helper), intent(in) :: this
            logical :: x
        end function

        pure module function oh_is_jac_defined(this) result(x)
            class(ode_helper), intent(in) :: this
            logical :: x
        end function

        module subroutine oh_define_constraints(this, n, fcn)
            class(ode_helper), intent(inout) :: this
            integer(int32), intent(in) :: n
            procedure(ode_constraint), intent(in), pointer :: fcn
        end subroutine

        pure module function oh_get_constraint_count(this) result(x)
            class(ode_helper), intent(in) :: this
            integer(int32) :: x
        end function

        pure module function oh_get_constraints_defined(this) result(x)
            class(ode_helper), intent(in) :: this
            logical :: x
        end function

        module subroutine oh_eval(this, x, y, dydx)
            class(ode_helper), intent(inout) :: this
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: dydx
        end subroutine

        module subroutine oh_eval_jac(this, x, y, jac)
            class(ode_helper), intent(inout) :: this
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:,:) :: jac
        end subroutine

        module subroutine oh_eval_constraints(this, x, y, f)
            class(ode_helper), intent(inout) :: this
            real(real64), intent(in) :: x
            real(real64), intent(in), dimension(:) :: y
            real(real64), intent(out), dimension(:) :: f
        end subroutine
    end interface

! ******************************************************************************
! INTEGRAL_ODE_INTEGRATOR.F90
! ------------------------------------------------------------------------------
    type, abstract :: ode_integrator
        !> @brief An array of relative error tolerance values for each ODE.
        real(real64), allocatable, dimension(:) :: m_rtol
        !> @brief An array of absolute error tolerance values for each ODE.
        real(real64), allocatable, dimension(:) :: m_atol
        !> @brief Determines if output should be provided wherever the
        !! integrator computes the solution as well as at user-specified points.
        logical :: m_allOutput = .true.
        !> @brief Determines if the integrator is allowed to overshoot a
        !! critical point and interpolate back to achieve output at
        !! the correct location.
        logical :: m_canOvershoot = .true.
        !> @brief The critical point that cannot be overshot.  This is only
        !! utilized in the event that m_canOvershoot is set to true.
        real(real64) :: m_criticalPoint = 0.0d0
        !> @brief Defines the minimum buffer size
        integer(int32) :: m_minBufferSize = 500
        !> @brief The maximum allowable step size.  This is only used if
        !! m_limitStepSize is set to true.
        real(real64) :: m_maxStepSize = 1.0d0
        !> @brief Determines if limits should be imposed upon maximum step
        !! size.
        logical :: m_limitStepSize = .false.
        !> @brief Determines the iteration step limit per integration step
        integer(int32) :: m_maxStepCount = 500
    contains
        !> @brief Performs the integration.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! logical function integrate(class(ode_integrator) this, class(ode_helper) fcnobj, real(real64) x(:), real(real64) y(:), optional class(errors) err)
        !! @endcode
        !!
        !! @param[in,out] this The ode_integrator object.
        !! @param[in,out] fcnobj The ode_helper object containing the equations
        !!  to integrate.
        !! @param[in] x An array containing the values of the independent
        !!  variable at which the solution is desired.  There must be at least
        !!  two values in this array.
        !! @param[in] y An N-element array containing the initial conditions for
        !!  each of the N ODEs.
        !! @param[in,out] err An optional output that can be used to provide
        !!  an error handling mechanism.  If not provided, a default error
        !!  handling mechanism will be utilized.  Possible errors that may
        !!  be encountered are as follows.
        !!  - INT_INVALID_INPUT_ERROR: An invalid input was supplied.
        !!  - INT_OUT_OF_MEMORY_ERROR: There is insufficient memory available.
        !!  - INT_LACK_OF_DEFINITION_ERROR: Occurs if no equations have been
        !!      defined.
        !!  - INT_ARRAY_SIZE_MISMATCH_ERROR: Occurs if @p y is not sized to
        !!      match the problem as defined in @p fcnobj, or if the tolerance
        !!      arrays are not sized to match the problem as defined in
        !!      @p fcnobj.
        !!  Notice, specific integrators may have additional errors.  See the
        !!  @p step routine of the appropriate integrator for more information.
        !!
        !! @return Returns the solution in a matrix of N+1 columns.  The
        !!  first column contains the values of the independent variable at
        !!  which the solution was computed.  The remaining columns contain the
        !!  solution points for each ODE.
        procedure, public :: integrate => oi_integrate
        !> @brief Gets a value determining if all integrator output should be
        !! reported.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical function get_provide_all_output(class(ode_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The ode_integrator object.
        !! @return Returns true if all integrator output should be reported;
        !!  else, returns false.
        procedure, public :: get_provide_all_output => oi_get_use_all_output
        !> @brief Sets a value determining if all integrator output should be
        !! reported.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_provide_all_output(class(ode_integrator) this, logical x)
        !! @endcode
        !!
        !! @param[in,out] this The ode_integrator object.
        !! @param[in] x Set to true if all integrator output should be
        !!  reported; else, set to false to only allow output at user-defined
        !!  points.
        procedure, public :: set_provide_all_output => oi_set_use_all_output
        !> @brief Gets a value determining if the integrator can overshoot the
        !! integration limits and interpolate back to achieve the desired
        !! solution point.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical function get_allow_overshoot(class(ode_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The ode_integrator object.
        !! @return Returns true if the integrator can overshoot; else, false.
        procedure, public :: get_allow_overshoot => oi_get_allow_overshoot
        !> @brief Sets a value determining if the integrator can overshoot the
        !! integration limits and interpolate back to achieve the desired
        !! solution point.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_allow_overshoot(class(ode_integrator) this, logical x)
        !! @endcode
        !!
        !! @param[in,out] this The ode_integrator object.
        !! @param[in] x Set to true to allow the integrator to overshoot the
        !!  integration limits; else, false to prevent overshoot.
        procedure, public :: set_allow_overshoot => oi_set_allow_overshoot
        !> @brief Gets a critical integration limit that the integrator is not
        !! allowed to step beyond.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure real(real64) function get_integration_limit(class(ode_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The ode_integrator object.
        !! @return The integration limit.
        procedure, public :: get_integration_limit => oi_get_critical_point
        !> @brief Sets a critical integration limit that the integrator is not
        !! allowed to step beyond.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_integration_limit(class(ode_integrator) this, real(real64) x)
        !! @endcode
        !!
        !! @param[in,out] this The ode_integrator object.
        !! @param[in] x The integration limit.
        procedure, public :: set_integration_limit => oi_set_critical_point
        !> @brief Gets the minimum internal storage buffer size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) get_min_buffer_size(class(ode_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The ode_integrator object.
        !! @return The minimum buffer size.
        procedure, public :: get_min_buffer_size => oi_get_min_buffer_size
        !> @brief Sets the minimum internal storage buffer size.  Properly
        !! sizing the buffer has the potential of improving the performance of
        !! the integrator.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_min_buffer_size(class(ode_integrator) this, integer(int32) x)
        !! @endcode
        !!
        !! @param[in,out] this The ode_integrator object.
        !! @param[in] x The minimum buffer size.
        procedure, public :: set_min_buffer_size => oi_set_min_buffer_size
        !> @brief Gets the maximum allowed step size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure real(real64) function get_max_step_size(class(ode_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The ode_integrator object.
        !! @return The maximum step size.
        procedure, public :: get_max_step_size => oi_get_max_step_size
        !> @brief Sets the maximum allowed step size.  This value is only
        !! honored if @p get_limit_step_size is set to true.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_max_step_size(class(ode_integrator) this, real(real64) x)
        !! @endcode
        !!
        !! @param[in,out] this The ode_integrator object.
        !! @param[in] x The maximum step size.
        procedure, public :: set_max_step_size => oi_set_max_step_size
        !> @brief Gets a value determining if the step size should be limited.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure logical get_limit_step_size(class(ode_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The ode_integrator object.
        !! @return Returns true if the integrator's step size should be limited;
        !!  else, false.
        procedure, public :: get_limit_step_size => oi_get_limit_step_size
        !> @brief Sets a value determining if the step size should be limited.
        !! This value must be set to true to use the value defined by
        !! @p set_max_step_size.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_limit_step_size(class(ode_integrator) this, logical x)
        !! @endcode
        !!
        !! @param[in,out] this The ode_integrator object.
        !! @param[in] x Set to true if the integrator's step size should be
        !!  limited; else, set to false.
        procedure, public :: set_limit_step_size => oi_set_limit_step_size
        !> @brief Gets the limit on the number of iterations allowed per
        !! integration step.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! pure integer(int32) function get_iteration_limit(class(ode_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The ode_integrator object.
        !! @return The iteration limit.
        procedure, public :: get_iteration_limit => oi_get_iteration_limit
        !> @brief Sets the limit on the number of iterations allowed per
        !! integration step.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_iteration_limit(class(ode_integrator) this, integer(int32) x)
        !! @endcode
        !!
        !! @param[in,out] this The ode_integrator object.
        !! @param[in] x The iteration limit.
        procedure, public :: set_iteration_limit => oi_set_iteration_limit
        !> @brief Gets the relative tolerances for each ODE.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64)(:) function get_relative_tolerances(class(ode_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The ode_integrator object.
        !! @return An array containing the tolerance value for each ODE.
        procedure, public :: get_relative_tolerances => oi_get_rtol
        !> @brief Sets the relative tolerances for each ODE.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_relative_tolerances(class(ode_integrator) this, real(real64) x(:))
        !! @endcode
        !!
        !! @param[in,out] this The ode_integrator object.
        !! @param[in] x An array containing the tolerance value for each ODE.
        procedure, public :: set_relative_tolerances => oi_set_rtol
        !> @brief Gets the absolute tolerances for each ODE.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64)(:) function get_absolute_tolerances(class(ode_integrator) this)
        !! @endcode
        !!
        !! @param[in] this The ode_integrator object.
        !! @return An array containing the tolerance value for each ODE.
        procedure, public :: get_absolute_tolerances => oi_get_atol
        !> @brief Sets the absolute tolerances for each ODE.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine set_absolute_tolerances(class(ode_integrator) this, real(real64) x(:))
        !! @endcode
        !!
        !! @param[in,out] this The ode_integrator object.
        !! @param[in] x An array containing the tolerance value for each ODE.
        procedure, public :: set_absolute_tolerances => oi_set_atol
        !> @brief Takes a single integration step towards the desired point.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! logical function step(class(ode_integrator) this, class(ode_helper) fcn, real(real64) x, real(real64) y(:), real(real64) xout, real(real64) rtol(:), real(real64) atol(:), optional class(errors) err)
        !! @endcode
        !!
        !! @param[in,out] this The ode_integrator object.
        !! @param[in,out] fcn An ode_helper object containing the ODEs to
        !!  integrate.
        !! @param[in,out] x On input, the value of the independent variable at
        !!  which to start.  On output, the value of the independent variable at
        !!  which integration terminated.
        !! @param[in,out] y On input, the value(s) of the dependent variable(s)
        !!  at the initial value given in @p x.  On output, the value(s) of the
        !!  dependent variable(s) as evaluated at the output given in @p x.
        !! @param[in] xout The value of the independent variable at which the
        !!  solution is desired.
        !! @param[in] rtol An array containing relative tolerance information
        !!  for each ODE.
        !! @param[in] atol An array containing absolute tolerance information
        !!  for each ODE.
        !! @param[in,out] err An optional argument that can be used to control
        !!  the error handling behavior of the integrator.
        !! @return Returns true if the integrator requests a stop; else, false,
        !!  to continue as normal.  A possible sitaution resulting in a true
        !!  value is in the event a constraint has been satisfied.
        procedure(ode_integrator_interface), public, deferred, pass :: step
        !> @brief Resets the state of the integrator.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine reset(class(ode_integrator) this)
        !! @endcode
        !!
        !! @param[in,out] this The ode_integrator object.
        procedure(ode_integrator_reset), public, deferred, pass :: reset
    end type

! ------------------------------------------------------------------------------
    interface
        !> @brief Defines a routine for computing a single integration step in
        !! the direction of @p xout.
        !!
        !! @param[in,out] this The ode_integrator object.
        !! @param[in,out] fcn An ode_helper object containing the ODEs to
        !!  integrate.
        !! @param[in,out] x On input, the value of the independent variable at
        !!  which to start.  On output, the value of the independent variable at
        !!  which integration terminated.
        !! @param[in,out] y On input, the value(s) of the dependent variable(s)
        !!  at the initial value given in @p x.  On output, the value(s) of the
        !!  dependent variable(s) as evaluated at the output given in @p x.
        !! @param[in] xout The value of the independent variable at which the
        !!  solution is desired.
        !! @param[in] rtol An array containing relative tolerance information
        !!  for each ODE.
        !! @param[in] atol An array containing absolute tolerance information
        !!  for each ODE.
        !! @param[in,out] err An optional argument that can be used to control
        !!  the error handling behavior of the integrator.
        !! @return Returns true if the integrator requests a stop; else, false,
        !!  to continue as normal.
        function ode_integrator_interface(this, fcn, x, y, xout, rtol, atol, err) result(brk)
            use, intrinsic :: iso_fortran_env, only : int32, real64
            use ferror
            import ode_integrator
            import ode_helper
            class(ode_integrator), intent(inout) :: this
            class(ode_helper), intent(inout) :: fcn
            real(real64), intent(inout) :: x
            real(real64), intent(inout), dimension(:) :: y
            real(real64), intent(in) :: xout
            real(real64), intent(in), dimension(:) :: rtol, atol
            class(errors), intent(inout), optional, target :: err
            logical :: brk
        end function

        !> @brief Resets the state of the integrator.
        !!
        !! @param[in,out] this The ode_integrator object.
        subroutine ode_integrator_reset(this)
            import ode_integrator
            class(ode_integrator), intent(inout) :: this
        end subroutine

        module function oi_integrate(this, fcnobj, x, y, err) result(rst)
            class(ode_integrator), intent(inout) :: this
            class(ode_helper), intent(inout) :: fcnobj
            real(real64), intent(in), dimension(:) :: x, y
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable, dimension(:,:) :: rst
        end function

        pure module function oi_get_use_all_output(this) result(x)
            class(ode_integrator), intent(in) :: this
            logical :: x
        end function

        module subroutine oi_set_use_all_output(this, x)
            class(ode_integrator), intent(inout) :: this
            logical, intent(in) :: x
        end subroutine

        pure module function oi_get_allow_overshoot(this) result(x)
            class(ode_integrator), intent(in) :: this
            logical :: x
        end function

        module subroutine oi_set_allow_overshoot(this, x)
            class(ode_integrator), intent(inout) :: this
            logical, intent(in) :: x
        end subroutine

        pure module function oi_get_critical_point(this) result(x)
            class(ode_integrator), intent(in) :: this
            real(real64) :: x
        end function

        module subroutine oi_set_critical_point(this, x)
            class(ode_integrator), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        pure module function oi_get_min_buffer_size(this) result(x)
            class(ode_integrator), intent(in) :: this
            integer(int32) :: x
        end function

        module subroutine oi_set_min_buffer_size(this, x)
            class(ode_integrator), intent(inout) :: this
            integer(int32), intent(in) :: x
        end subroutine

        pure module function oi_get_max_step_size(this) result(x)
            class(ode_integrator), intent(in) :: this
            real(real64) :: x
        end function

        module subroutine oi_set_max_step_size(this, x)
            class(ode_integrator), intent(inout) :: this
            real(real64), intent(in) :: x
        end subroutine

        pure module function oi_get_limit_step_size(this) result(x)
            class(ode_integrator), intent(in) :: this
            logical :: x
        end function

        module subroutine oi_set_limit_step_size(this, x)
            class(ode_integrator), intent(inout) :: this
            logical, intent(in) :: x
        end subroutine

        pure module function oi_get_iteration_limit(this) result(x)
            class(ode_integrator), intent(in) :: this
            integer(int32) :: x
        end function

        module subroutine oi_set_iteration_limit(this, x)
            class(ode_integrator), intent(inout) :: this
            integer(int32), intent(in) :: x
        end subroutine

        module function oi_get_rtol(this) result(x)
            class(ode_integrator), intent(in) :: this
            real(real64), allocatable, dimension(:) :: x
        end function

        module subroutine oi_set_rtol(this, x)
            class(ode_integrator), intent(inout) :: this
            real(real64), intent(in), dimension(:) :: x
        end subroutine

        module function oi_get_atol(this) result(x)
            class(ode_integrator), intent(in) :: this
            real(real64), allocatable, dimension(:) :: x
        end function

        module subroutine oi_set_atol(this, x)
            class(ode_integrator), intent(inout) :: this
            real(real64), intent(in), dimension(:) :: x
        end subroutine
    end interface

! ******************************************************************************
! INTEGRAL_ODE_AUTO.F90
! ------------------------------------------------------------------------------
    !> @brief Defines an integrator for systems of first order ODEs that is
    !! capable of switching between an Adams method and a BDF method
    !! automatically.  This integrator is able to handle both stiff and
    !! non-stiff systems of equations.
    !!
    !! @par Example
    !! The following example illustrates the use of this integrator to solve
    !! the Van Der Pol equation, which is a second-order ODE that does exhibit
    !! difficult and stiff behavior.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use integral_core
    !!     use fplot_core
    !!     implicit none
    !!
    !!     ! Local Variables
    !!     type(ode_helper) :: fcn
    !!     type(ode_auto) :: integrator
    !!     procedure(ode_fcn), pointer :: ptr
    !!     real(real64) :: ic(2), t(2)
    !!     real(real64), allocatable, dimension(:,:) :: x
    !!     type(plot_2d) :: plt
    !!     type(plot_data_2d) :: d1
    !!     class(plot_axis), pointer :: xAxis, yAxis
    !!     class(legend), pointer :: lgnd
    !!
    !!     ! Set up the integrator
    !!     ptr => vdp
    !!     call fcn%define_equations(2, ptr)
    !!
    !!     ! Define the initial conditions
    !!     t = [0.0d0, 8.0d1]
    !!     ic = [2.0d0, 0.0d0]
    !!
    !!     ! Compute the solution
    !!     x = integrator%integrate(fcn, t, ic)
    !!
    !!     ! Plot the solution
    !!     call plt%initialize()
    !!     call plt%set_font_size(14)
    !!
    !!     lgnd => plt%get_legend()
    !!     call lgnd%set_is_visible(.false.)
    !!
    !!     xAxis => plt%get_x_axis()
    !!     call xAxis%set_title("t")
    !!
    !!     yAxis => plt%get_y_axis()
    !!     call yAxis%set_title("x(t)")
    !!
    !!     call d1%set_name("x(t)")
    !!     call d1%set_line_width(2.0)
    !!     call d1%set_line_color(CLR_BLUE)
    !!     call d1%define_data(x(:,1), x(:,2))
    !!
    !!     call plt%push(d1)
    !!     call plt%draw()
    !!
    !! contains
    !!     ! Van Der Pol Equation
    !!     ! x" + x - mu * (1 - x**2) * x' = 0
    !!     subroutine vdp(t, x, dxdt)
    !!         real(real64), intent(in) :: t
    !!         real(real64), intent(in), dimension(:) :: x
    !!         real(real64), intent(out), dimension(:) :: dxdt
    !!
    !!         real(real64), parameter :: mu = 20.0d0
    !!
    !!         dxdt(1) = x(2)
    !!         dxdt(2) = mu * (1.0d0 - x(1)**2) * x(2) - x(1)
    !!     end subroutine
    !! end program
    !! @endcode
    !! The above program produces the following output.
    !! @image html vanderpol_example.png
    type, extends(ode_integrator) :: ode_auto
        !> A workspace array.
        real(real64), allocatable, dimension(:) :: m_rwork
        !> An integer workspace array.
        integer(int32), allocatable, dimension(:) :: m_iwork
        !> An array used to contain information regarding constraint status.
        !! A value of 1 indicates the constraint was satisfied; else, a value
        !! of 0 indicates the constraint wasn't satisfied.  There is one entry
        !! in the array for each constraint.
        integer(int32), allocatable, dimension(:) :: m_rststats
        !> This flag is used directly by ODEPACK.  Set to 1 for initial call.
        integer(int32) :: m_istate = 1
    contains
        !> @brief Takes a single integration step towards the desired point.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! logical function step(class(ode_auto) this, class(ode_helper) fcn, real(real64) x, real(real64) y(:), real(real64) xout, real(real64) rtol(:), real(real64) atol(:), optional class(errors) err)
        !! @endcode
        !!
        !! @param[in,out] this The ode_auto object.
        !! @param[in,out] fcn An ode_helper object containing the ODEs to
        !!  integrate.
        !! @param[in,out] x On input, the value of the independent variable at
        !!  which to start.  On output, the value of the independent variable at
        !!  which integration terminated.
        !! @param[in,out] y On input, the value(s) of the dependent variable(s)
        !!  at the initial value given in @p x.  On output, the value(s) of the
        !!  dependent variable(s) as evaluated at the output given in @p x.
        !! @param[in] xout The value of the independent variable at which the
        !!  solution is desired.
        !! @param[in] rtol An array containing relative tolerance information
        !!  for each ODE.
        !! @param[in] atol An array containing absolute tolerance information
        !!  for each ODE.
        !! @param[in,out] err An optional output that can be used to provide
        !!  an error handling mechanism.  If not provided, a default error
        !!  handling mechanism will be utilized.  Possible errors that may
        !!  be encountered are as follows.
        !!  - INT_EXCESSIVE_WORK_ERROR: Occurs if excessive work has been done.
        !!  - INT_IMPRACTICAL_TOLERANCE_ERROR: Occurs if the user-defined
        !!      tolerances are not practically achievable.
        !!  - INT_INVALID_INPUT_ERROR: Occurs if an invalid input was supplied.
        !!  - INT_REPEATED_ERROR_TEST_FAILURE: Occurs if error testing
        !!      repeatadly fails.
        !!  - INT_CONVERGENCE_ERROR: Occurs if the iteration process cannot
        !!      converge.
        !!  - INT_INTEGRAND_BEHAVIOR_ERROR: Occurs if a component of the
        !!      solution appears to vanish.
        !!
        !! @return Returns true if the integrator requests a stop; else, false,
        !!  to continue as normal.  A possible sitaution resulting in a true
        !!  value is in the event a constraint has been satisfied.
        procedure, public :: step => oa_step
        !> @brief Resets the state of the integrator.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! subroutine reset(class(ode_auto) this)
        !! @endcode
        !!
        !! @param[in,out] this The ode_auto object.
        procedure, public :: reset => oa_reset_integrator
        !> @brief Gets the status of each constraint equation.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! logical(:) function get_constraint_status(class(ode_auto) this)
        !! @endcode
        !!
        !! @param[in] this The ode_auto object.
        !! @return An array containing the status of each constraint equation.
        !!  A value of true indicates that the constraint was satisfied; else,
        !!  a value of false indicates the constraint was not satisfied.
        procedure, public :: get_constraint_status => oa_get_constraint_info

        procedure, private :: init_workspace => oa_init_workspace
    end type

    interface
        module function oa_step(this, fcn, x, y, xout, rtol, atol, err) result(brk)
            class(ode_auto), intent(inout) :: this
            class(ode_helper), intent(inout) :: fcn
            real(real64), intent(inout) :: x
            real(real64), intent(inout), dimension(:) :: y
            real(real64), intent(in) :: xout
            real(real64), intent(in), dimension(:) :: rtol, atol
            class(errors), intent(inout), optional, target :: err
            logical :: brk
        end function

        module subroutine oa_init_workspace(this, liw, lrw, ncnsts, err)
            class(ode_auto), intent(inout) :: this
            integer(int32), intent(in) :: liw, lrw, ncnsts
            class(errors), intent(inout) :: err
        end subroutine

        module subroutine oa_reset_integrator(this)
            class(ode_auto), intent(inout) :: this
        end subroutine

        module function oa_get_constraint_info(this) result(x)
            class(ode_auto), intent(in) :: this
            logical, allocatable, dimension(:) :: x
        end function
    end interface
end module
