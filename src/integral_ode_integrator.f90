! integral_ode_integrator.f90

submodule (integral_core) integral_ode_integrator
contains
! ------------------------------------------------------------------------------
    module function oi_integrate(this, fcnobj, x, y, err) result(rst)
        ! Arguments
        class(ode_integrator), intent(inout) :: this
        class(ode_helper), intent(in) :: fcnobj
        real(real64), intent(in), dimension(:) :: x, y
        class(errors), intent(inout), optional, target :: err
        real(real64), allocatable, dimension(:,:) :: rst

        ! Local Variables
        integer(int32) :: neqn, ncols, flag, i
        real(real64), allocatable, dimension(:,:) :: buffer
        real(real64), allocatable, dimension(:) :: ytemp
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        neqn = fcnobj%get_equation_count()
        ncols = neqn + 1
        allocate(ytemp(neqn), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("oi_integrate", &
                "Insufficient memory available.", INT_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Quick Return
        if (neqn <= 0) return

        ! Input Checking
        if (fcnobj%get_equations_defined()) then
            ! ERROR: Equations undefined
        end if
        if (neqn /= size(y)) then
            ! ERROR: # of equations does not match initial condition vector size
        end if
        if (size(x) < 2) then
            ! ERROR: There must be at least a starting and ending point for the
            ! integration routine
        end if

        ! Ensure tolerances are defined.  If not, utilize defaults
        if (allocated(this%m_rtol)) then
            if (size(this%m_rtol) /= neqn) then
                ! WARNING: Size mismatch in tolerance array - using defaults
                deallocate(this%m_rtol)
                this%m_rtol = alloc_default_rtol(neqn)
            end if
        else
            this%m_rtol = alloc_default_rtol(neqn)
        end if

        if (allocated(this%m_atol)) then
            if (size(this%m_rtol) /= neqn) then
                ! WARNING: Size mismatch in tolerance array - using defaults
                deallocate(this%m_atol)
                this%m_atol = alloc_default_atol(neqn)
            end if
        else
            this%m_atol = alloc_default_atol(neqn)
        end if

        !
        if (this%get_provide_all_output()) then
            ! The integrator is expected to provide all outputs, regardless if
            ! the actual output lies on a requested output point
            allocate(buffer(this%get_min_buffer_size(), ncols), stat = flag)
            if (flag /= 0) then
                call errmgr%report_error("oi_integrate", &
                    "Insufficient memory available.", INT_OUT_OF_MEMORY_ERROR)
                return
            end if

            ! Store the initial conditions
            buffer(1,1) = x(1)
            buffer(1,2:ncols) = y

            ! Integrate
            i = 1
            do
                ! Copy the previous step's output
                ytemp = buffer(i,2:ncols)

                ! Increment i
                i = i + 1

                ! Ensure the buffer is sufficiently sized
                if (i > size(buffer, 1)) then
                    call realloc_buffer(buffer, &
                        size(buffer, 1) + this%get_min_buffer_size(), errmgr)
                    if (errmgr%has_error_occurred()) return
                end if

                ! Take the next integration step

                ! Store the output

                ! Break if X exceeds XOUT
            end do
        else
        end if

        ! Resize the buffer to tightly fit the data, if necessary
    end function

! ------------------------------------------------------------------------------
    pure module function oi_get_use_all_output(this) result(x)
        class(ode_integrator), intent(in) :: this
        logical :: x
    end function

! --------------------
    module subroutine oi_set_use_all_output(this, x)
        class(ode_integrator), intent(inout) :: this
        logical, intent(in) :: x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oi_get_allow_overshoot(this) result(x)
        class(ode_integrator), intent(in) :: this
        logical :: x
    end function

! --------------------
    module subroutine oi_set_allow_overshoot(this, x)
        class(ode_integrator), intent(inout) :: this
        logical, intent(in) :: x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oi_get_critical_point(this) result(x)
        class(ode_integrator), intent(in) :: this
        real(real64) :: x
    end function

! --------------------
    module subroutine oi_set_critical_point(this, x)
        class(ode_integrator), intent(inout) :: this
        real(real64), intent(in) :: x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oi_get_min_buffer_size(this) result(x)
        class(ode_integrator), intent(in) :: this
        integer(int32) :: x
    end function

! --------------------
    module subroutine oi_set_min_buffer_size(this, x)
        class(ode_integrator), intent(inout) :: this
        integer(int32), intent(in) :: x
    end subroutine

! ******************************************************************************
! PRIVATE SUBMODULE ROUTINES
! ------------------------------------------------------------------------------
    ! Establishes a default relative tolerance array.
    function alloc_default_rtol(n) result(x)
        integer(int32), intent(in) :: n
        real(real64), dimension(n) :: x
        x = 1.0d-6
    end function

! ------------------------------------------------------------------------------
    ! Establishes a default absolute tolerance array.
    function alloc_default_atol(n) result(x)
        integer(int32), intent(in) :: n
        real(real64), dimension(n) :: x
        x = 1.0d-8
    end function

! ------------------------------------------------------------------------------
    subroutine realloc_buffer(x, new_row_count, err)
        ! Arguments
        real(real64), intent(inout), allocatable, dimension(:,:) :: xold
        integer(int32), intent(in) :: new_row_count
        class(errors), intent(inout) :: err

        ! Local Variables
        integer(int32) :: oldm, n
        real(real64), allocatable, dimension(:,:) :: temp

        ! Initialization
        oldm = size(x, 1)
        n = size(x, 2)

        ! If oldm == new_row_count, do nothing
        if (oldm == new_row_count) return

        ! Allocate space for a temporary array that allows for resizing of X
        allocate(temp(oldm, n), stat = flag)
        if (flag /= 0) then
            call err%report_error("realloc_buffer", &
                "Insufficient memory available.", INT_OUT_OF_MEMORY_ERROR)
            return
        end if
        temp = x

        ! Deallocate X, and resize
        deallocate(x)
        allocate(x(new_row_count, n), stat = flag)
        if (flag /= 0) then
            call err%report_error("realloc_buffer", &
                "Insufficient memory available.", INT_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Copy temp back into X
        if (new_row_count > oldm) then
            x(1:oldm,:) = temp
        else
            x = temp(1:new_row_count,:)
        end if
    end subroutine

! ------------------------------------------------------------------------------
end submodule
