! integral_ode_integrator.f90

submodule (integral_core) integral_ode_integrator
contains
! ------------------------------------------------------------------------------
    module function oi_integrate(this, fcnobj, x, y, err) result(rst)
        ! Arguments
        class(ode_integrator), intent(inout) :: this
        class(ode_helper), intent(inout) :: fcnobj
        real(real64), intent(in), dimension(:) :: x, y
        class(errors), intent(inout), optional, target :: err
        real(real64), allocatable, dimension(:,:) :: rst

        ! Local Variables
        integer(int32) :: neqn, ncols, flag, i, n, nbuffer
        real(real64), allocatable, dimension(:,:) :: buffer
        real(real64), allocatable, dimension(:) :: ytemp
        real(real64) :: xout, xi
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = 256) :: errmsg
        logical :: brk, ascending

        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        neqn = fcnobj%get_equation_count()
        ncols = neqn + 1
        n = size(x)
        if (this%get_provide_all_output()) then
            nbuffer = this%get_min_buffer_size()
        else
            nbuffer = n
        end if
        allocate(ytemp(neqn), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("oi_integrate", &
                "Insufficient memory available.", INT_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Quick Return
        if (neqn <= 0) return

        ! Input Checking
        if (.not.fcnobj%get_equations_defined()) then
            ! ERROR: Equations undefined
            call errmgr%report_error("oi_integrate", &
                "The routine expected to hold the ODEs is undefined.", &
                INT_LACK_OF_DEFINITION_ERROR)
            return
        end if
        if (neqn /= size(y)) then
            ! ERROR: # of equations does not match initial condition vector size
            write(errmsg, '(AI0AI0A)') "Expected an array of size ", neqn, &
                ", but found an array of size ", size(y), "."
            call errmgr%report_error("oi_integrate", trim(errmsg), &
                INT_ARRAY_SIZE_MISMATCH_ERROR)
            return
        end if
        if (n < 2) then
            ! ERROR: There must be at least a starting and ending point for the
            ! integration routine
        end if
        if (x(n) == x(1)) then
            ! ERROR: No integration range
            call errmgr%report_error("oi_integrate", &
                "The starting and ending integration points are the same.", &
                INT_INVALID_INPUT_ERROR)
            return
        end if

        ! Ensure tolerances are defined.  If not, utilize defaults
        if (allocated(this%m_rtol)) then
            if (size(this%m_rtol) /= neqn) then
                ! WARNING: Size mismatch in tolerance array - using defaults
                call errmgr%report_warning("oi_integrate", &
                    "The relative tolerance array is not correctly " // &
                    "sized for the problem.  Using a default tolerance " // &
                    "instead.", INT_ARRAY_SIZE_MISMATCH_ERROR)
                deallocate(this%m_rtol)
                this%m_rtol = alloc_default_rtol(neqn)
            end if
        else
            this%m_rtol = alloc_default_rtol(neqn)
        end if

        if (allocated(this%m_atol)) then
            if (size(this%m_rtol) /= neqn) then
                ! WARNING: Size mismatch in tolerance array - using defaults
                call errmgr%report_warning("oi_integrate", &
                    "The absolute tolerance array is not correctly " // &
                    "sized for the problem.  Using a default tolerance " // &
                    "instead.", INT_ARRAY_SIZE_MISMATCH_ERROR)
                deallocate(this%m_atol)
                this%m_atol = alloc_default_atol(neqn)
            end if
        else
            this%m_atol = alloc_default_atol(neqn)
        end if

        ! Additional Initialization
        allocate(buffer(nbuffer, ncols), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("oi_integrate", &
                "Insufficient memory available.", INT_OUT_OF_MEMORY_ERROR)
            return
        end if
        ascending = x(n) > x(1)

        ! Store the initial conditions
        buffer(1,1) = x(1)
        buffer(1,2:ncols) = y

        ! Integrate
        i = 1
        if (this%get_provide_all_output()) then
            xout = x(n)
            call this%set_integration_limit(x(n))
            call this%set_allow_overshoot(.false.)
        end if
        do
            ! Copy the previous step's output
            xi = buffer(i,1)
            ytemp = buffer(i,2:ncols)

            ! Increment i
            i = i + 1

            ! Ensure the buffer is sufficiently sized
            if (i > size(buffer, 1)) then
                call realloc_buffer(buffer, &
                    size(buffer, 1) + nbuffer, errmgr)
                if (errmgr%has_error_occurred()) return
            end if

            ! Take the next integration step
            if (.not.this%get_provide_all_output()) xout = x(i)
            brk = this%step(fcnobj, xi, ytemp, xout, this%m_rtol, this%m_atol, &
                errmgr)
            if (brk) then
                ! Store the output, and then exit
                buffer(i,1) = xi
                buffer(i,2:ncols) = ytemp
                exit
            end if
            if (errmgr%has_error_occurred()) return

            ! Store the output
            buffer(i,1) = xi
            buffer(i,2:ncols) = ytemp

            ! Break if X exceeds XOUT
            if (ascending .and. xi >= x(n)) then
                exit
            else if (.not.ascending .and. xi <= x(n)) then
                exit
            end if
        end do

        ! Resize the buffer to tightly fit the data, if necessary
        rst = buffer(1:i,:)
    end function

! ------------------------------------------------------------------------------
    pure module function oi_get_use_all_output(this) result(x)
        class(ode_integrator), intent(in) :: this
        logical :: x
        x = this%m_allOutput
    end function

! --------------------
    module subroutine oi_set_use_all_output(this, x)
        class(ode_integrator), intent(inout) :: this
        logical, intent(in) :: x
        this%m_allOutput = x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oi_get_allow_overshoot(this) result(x)
        class(ode_integrator), intent(in) :: this
        logical :: x
        x = this%m_canOvershoot
    end function

! --------------------
    module subroutine oi_set_allow_overshoot(this, x)
        class(ode_integrator), intent(inout) :: this
        logical, intent(in) :: x
        this%m_canOvershoot = x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oi_get_critical_point(this) result(x)
        class(ode_integrator), intent(in) :: this
        real(real64) :: x
        x = this%m_criticalPoint
    end function

! --------------------
    module subroutine oi_set_critical_point(this, x)
        class(ode_integrator), intent(inout) :: this
        real(real64), intent(in) :: x
        this%m_criticalPoint = x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oi_get_min_buffer_size(this) result(x)
        class(ode_integrator), intent(in) :: this
        integer(int32) :: x
        x = this%m_minBufferSize
    end function

! --------------------
    module subroutine oi_set_min_buffer_size(this, x)
        class(ode_integrator), intent(inout) :: this
        integer(int32), intent(in) :: x
        this%m_minBufferSize = x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oi_get_max_step_size(this) result(x)
        class(ode_integrator), intent(in) :: this
        real(real64) :: x
        x = this%m_maxStepSize
    end function

! --------------------
    module subroutine oi_set_max_step_size(this, x)
        class(ode_integrator), intent(inout) :: this
        real(real64), intent(in) :: x
        this%m_maxStepSize = x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oi_get_limit_step_size(this) result(x)
        class(ode_integrator), intent(in) :: this
        logical :: x
        x = this%m_limitStepSize
    end function

! --------------------
    module subroutine oi_set_limit_step_size(this, x)
        class(ode_integrator), intent(inout) :: this
        logical, intent(in) :: x
        this%m_limitStepSize = x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oi_get_iteration_limit(this) result(x)
        class(ode_integrator), intent(in) :: this
        integer(int32) :: x
        x = this%m_maxStepCount
    end function

! --------------------
    module subroutine oi_set_iteration_limit(this, x)
        class(ode_integrator), intent(inout) :: this
        integer(int32), intent(in) :: x
        this%m_maxStepCount = x
    end subroutine

! ------------------------------------------------------------------------------
    module function oi_get_rtol(this) result(x)
        class(ode_integrator), intent(in) :: this
        real(real64), allocatable, dimension(:) :: x
        x = this%m_rtol
    end function

! --------------------
    module subroutine oi_set_rtol(this, x)
        class(ode_integrator), intent(inout) :: this
        real(real64), intent(in), dimension(:) :: x
        this%m_rtol = x
    end subroutine

! ------------------------------------------------------------------------------
    module function oi_get_atol(this) result(x)
        class(ode_integrator), intent(in) :: this
        real(real64), allocatable, dimension(:) :: x
        x = this%m_atol
    end function

! --------------------
    module subroutine oi_set_atol(this, x)
        class(ode_integrator), intent(inout) :: this
        real(real64), intent(in), dimension(:) :: x
        this%m_atol = x
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
        real(real64), intent(inout), allocatable, dimension(:,:) :: x
        integer(int32), intent(in) :: new_row_count
        class(errors), intent(inout) :: err

        ! Local Variables
        integer(int32) :: oldm, n, flag
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
