! integral_ode_helper.f90

submodule (integral_core) integral_ode_helper
contains
! ------------------------------------------------------------------------------
    module subroutine oh_init(this, neqn, fcn, jac)
        class(ode_helper), intent(inout) :: this
        integer(int32), intent(in) :: neqn
        procedure(ode_fcn), pointer, intent(in) :: fcn
        procedure(ode_jacobian), pointer, intent(in), optional :: jac
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oh_get_count(this) result(x)
        class(ode_helper), intent(in) :: this
        integer(int32) :: x
    end function

! ------------------------------------------------------------------------------
    pure module function oh_is_fcn_defined(this) result(x)
        class(ode_helper), intent(in) :: this
        logical :: x
    end function

! ------------------------------------------------------------------------------
    pure module function oh_is_jac_defined(this) result(x)
        class(ode_helper), intent(in) :: this
        logical :: x
    end function

! ------------------------------------------------------------------------------
end submodule
