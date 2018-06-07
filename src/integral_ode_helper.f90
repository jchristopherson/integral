! integral_ode_helper.f90

submodule (integral_core) integral_ode_helper
contains
! ------------------------------------------------------------------------------
    module subroutine oh_init(this, neqn, fcn, jac)
        class(ode_helper), intent(inout) :: this
        integer(int32), intent(in) :: neqn
        procedure(ode_fcn), pointer, intent(in) :: fcn
        procedure(ode_jacobian), pointer, intent(in), optional :: jac
        nullify(this%m_fcn)
        nullify(this%m_jac)
        this%m_count = neqn
        this%m_fcn => fcn
        if (present(jac)) this%m_jac => jac
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oh_get_count(this) result(x)
        class(ode_helper), intent(in) :: this
        integer(int32) :: x
        x = this%m_count
    end function

! ------------------------------------------------------------------------------
    pure module function oh_is_fcn_defined(this) result(x)
        class(ode_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_fcn)
    end function

! ------------------------------------------------------------------------------
    pure module function oh_is_jac_defined(this) result(x)
        class(ode_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_jac)
    end function

! ------------------------------------------------------------------------------
    module subroutine oh_define_constraints(this, n, fcn)
        class(ode_helper), intent(inout) :: this
        integer(int32), intent(in) :: n
        procedure(ode_constraint), intent(in), pointer :: fcn
        this%m_constraints = n
        this%m_rts => fcn
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oh_get_constraint_count(this) result(x)
        class(ode_helper), intent(in) :: this
        integer(int32) :: x
        x = this%m_constraints
    end function

! ------------------------------------------------------------------------------
    pure module function oh_get_constraints_defined(this) result(x)
        class(ode_helper), intent(in) :: this
        logical :: x
        x = associated(this%m_rts)
    end function

! ------------------------------------------------------------------------------
    module subroutine oh_eval(this, x, y, dydx)
        class(ode_helper), intent(inout) :: this
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: y
        real(real64), intent(out), dimension(:) :: dydx
        call this%m_fcn(x, y, dydx)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine oh_eval_jac(this, x, y, jac)
        class(ode_helper), intent(inout) :: this
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: y
        real(real64), intent(out), dimension(:,:) :: jac
        call this%m_jac(x, y, jac)
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine oh_eval_constraints(this, x, y, f)
        class(ode_helper), intent(inout) :: this
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: y
        real(real64), intent(out), dimension(:) :: f
        call this%m_rts(x, y, f)
    end subroutine

! ------------------------------------------------------------------------------
end submodule
