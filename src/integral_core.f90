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
    integer(int32), parameter :: INT_EXCESSIVE_WORK_ERROR = 10
    integer(int32), parameter :: INT_IMPRACTICAL_TOLERANCE_ERROR = 11
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
    contains
        procedure, public :: integrate => oi_integrate
        procedure, public :: get_provide_all_output => oi_get_use_all_output
        procedure, public :: set_provide_all_output => oi_set_use_all_output
        procedure, public :: get_allow_overshoot => oi_get_allow_overshoot
        procedure, public :: set_allow_overshoot => oi_set_allow_overshoot
        procedure, public :: get_integration_limit => oi_get_critical_point
        procedure, public :: set_integration_limit => oi_set_critical_point
        procedure, public :: get_min_buffer_size => oi_get_min_buffer_size
        procedure, public :: set_min_buffer_size => oi_set_min_buffer_size
        !> @brief Takes a single integration step towards the desired point.
        procedure(ode_integrator_interface), public, deferred, pass :: step
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
    end interface

! ******************************************************************************
! INTEGRAL_ODE_AUTO.F90
! ------------------------------------------------------------------------------
    type, extends(ode_integrator) :: ode_auto
        real(real64), allocatable, dimension(:) :: m_rwork
        integer(int32), allocatable, dimension(:) :: m_iwork
        integer(int32), allocatable, dimension(:) :: m_rststats
        !> This flag is used directly by ODEPACK.  Set to 1 for initial call.
        integer(int32) :: m_istate = 1
    contains
        procedure, public :: step => oa_step
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
    end interface
end module
