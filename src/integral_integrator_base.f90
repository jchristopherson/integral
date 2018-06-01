! integral_integrator_base.f90

submodule (integral_core) integral_integrator_base
contains
! ------------------------------------------------------------------------------
    pure module function ib_get_abs_tol(this) result(x)
        class(integrator_base), intent(in) :: this
        real(real64) :: x
        x = this%m_absTol
    end function

    module subroutine ib_set_abs_tol(this, x)
        class(integrator_base), intent(inout) :: this
        real(real64), intent(in) :: x
        this%m_absTol = x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function ib_get_rel_tol(this) result(x)
        class(integrator_base), intent(in) :: this
        real(real64) :: x
        x = this%m_relTol
    end function

    module subroutine ib_set_rel_tol(this, x)
        class(integrator_base), intent(inout) :: this
        real(real64), intent(in) :: x
        this%m_relTol = x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function ib_get_max_subintervals(this) result(x)
        class(integrator_base), intent(in) :: this
        integer(int32) :: x
        x = this%m_maxInt
    end function

    module subroutine ib_set_max_subintervals(this, x)
        class(integrator_base), intent(inout) :: this
        integer(int32), intent(in) :: x
        this%m_maxInt = x
    end subroutine

! ------------------------------------------------------------------------------
end submodule
