! integral_ode_irk.f90

submodule (integral_core) integral_ode_irk
contains
! ------------------------------------------------------------------------------
    module function oirk_step(this, fcn, x, y, xout, rtol, atol, err) result(brk)
        ! Arguments
        class(ode_irk), intent(inout) :: this
        class(ode_helper), intent(inout) :: fcn
        real(real64), intent(inout) :: x
        real(real64), intent(inout), dimension(:) :: y
        real(real64), intent(in) :: xout
        real(real64), intent(in), dimension(:) :: rtol, atol
        class(errors), intent(inout), optional, target :: err
        logical :: brk

        ! Local Variables
        integer(int32) :: neq, lrw, liw, itol, ijac, mljac, mujac, imas, &
            mlmas, mumas, iout, ipar, ljac, le, lmas, idid, ipar(1)
        real(real64) :: h, rpar(1)
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Set up error handling
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Initialization
        brk = .false.
        neq = fcn%get_equation_count()
        itol = 1
        mljac = neq
        mujac = neq
        if (fcn%get_jacobian_defined()) then
            ijac = 1
        else
            ijac = 0
        end if
        mlmas = neq
        mumas = neq
        if (fcn%get_mass_matrix_defined()) then
            imas = 1
        else
            imas = 0
        end if

        ! Establish workspace arrays
        if (mljac == neq) then
            ljac = neq
            le = neq
        else
            ljac = mljac + mujac + 1
            le = 2 * mljac + mujac + 1
        end if
        lmas = 0
        if (imas == 1 .and. mumas == neq) lmas = neq
        if (lmas < neq) lmas = mlmas + mumas + 1
        liw = 3 * neq + 20
        lrw = neq * (ljac + lmas + 3 * le + 12) + 20

        call this%init_workspace(liw, lrw, errmgr)
        if (errmgr%has_error_occurred()) return

        ! Set solution parameters
        this%m_iwork(2) = this%get_iteration_limit()
        this%m_iwork(3) = this%get_newton_iteration_count()
        if (this%get_index_1_dimension() == 0) then
            this%m_iwork(5) = neq
        else
            this%m_iwork(5) = this%get_index_1_dimension()
        end if
        this%m_iwork(6) = this%get_index_2_dimension()
        this%m_iwork(7) = this%get_index_3_dimension()
        this%m_iwork(8) = this%get_controller()
        this%m_iwork(9) = 0
        this%m_iwork(10) = 0
        this%m_rwork(4) = this%get_newton_iteration_tolerance()
        if (this%get_limit_step_size()) then
            this%m_rwork(7) = min(this%get_max_step_size(), xout - x)
        else
            this%m_rwork(7) = xout - x
        end if

        ! Determine if all output should be provided
        if (this%get_provide_all_output()) then
            iout = 1    ! Call radausolout
        else
            iout = 0    ! Don't call radausolout
        end if

        ! Initial step size estimate
        if (this%m_stepSize == 0.0d0) h = this%m_rwork(7)

        ! Compute the solution
        call RADAU5(neq, radaufcn, x, y, xout, h, rtol, atol, itol, radaujac, &
            ijac, mljac, mujac, radaumass, imas, mlmas, mumas, radausolout, &
            iout, this%m_rwork, lrw, this%m_iwork, liw, rpar, ipar, idid)

        ! Update the step-size tracking variable
        this%m_stepSize = h

        ! Check for errors
        select case (idid)
        case (-1)
            call errmgr%report_error("oirk_step", &
                "Invalid input encountered.", INT_INVALID_INPUT_ERROR)
            return
        case (-2)
            call errmgr%report_error("oirk_step", &
                "The step could not converge in the allowed number of " // &
                "iterations.  Try increasing the number of allowed " // &
                "iterations.", INT_CONVERGENCE_ERROR)
            return
        case (-3)
            call errmgr%report_error("oirk_step", &
                "The step size has become too small.", &
                INT_STEP_SIZE_TOO_SMALL_ERROR)
            return
        case (-4)
            call errmgr%report_error("oirk_step", &
                "The matrix has become singular.", &
                INT_SINGULAR_MATRIX_ERROR)
            return
        end select
    
    contains
        subroutine radaufcn(n_, x_, y_, f_, rpar_, ipar_)
            integer(int32), intent(in) :: n_
            real(real64), intent(in) :: x_, y_(n_)
            real(real64), intent(out) :: f_(n_)
            real(real64), intent(inout) :: rpar_(*)
            integer(int32), intent(inout) :: ipar_(*)
            call fcn%evaluate_ode(x_, y_, f_)
        end subroutine

        subroutine radaujac(n_, x_, y_, dfy_, ldfy_, rpar_, ipar_)
            integer(int32), intent(in) :: n_, ldfy_
            real(real64), intent(in) :: x_, y_(n_)
            real(real64), intent(out) :: dfy_(ldfy_, n_)
            real(real64), intent(inout) :: rpar_(*)
            integer(int32), intent(inout) :: ipar_(*)
            call fcn%evaluate_jacobian(x_, y_, dfy_)
        end subroutine

        subroutine radaumass(n_, am_, lmas_, rpar_, ipar_)
            integer(int32), intent(in) :: n_, lmas_
            real(real64), intent(out) :: am_(lmas_, n_)
            real(real64), intent(inout) :: rpar_(*)
            integer(int32), intent(inout) :: ipar_(*)
            call fcn%evaluate_mass_matrix(am_)
        end subroutine

        subroutine radausolout(nr_, xold_, x_, y_, cont_, lrc_, n_, rpar_, ipar_, irtrn_)
            integer(int32), intent(in) :: nr_, lrc_, n_
            real(real64), intent(in) :: x_, y_(n_), cont_(lrc_)
            real(real64), intent(inout) :: rpar_(*)
            integer(int32), intent(inout) :: ipar_(*)
            integer(int32), intent(out) :: irtrn_ ! set to < 0 to force exit

            ! This routine is only called when dense output is requested.  For
            ! dense output, call a routine from here to store the results.  If
            ! for some reason, the user wishes to break from the program set
            ! irtrn_ to a value < 0.
            
            if (associated(this%m_collector)) call this%m_collector(x_, y_)
        end subroutine
    end function

! ------------------------------------------------------------------------------
    module subroutine oirk_reset_integrator(this)
        ! Arguments
        class(ode_irk), intent(inout) :: this

        ! Local Variables
        integer(int32) :: liw, lrw

        ! Initialization
        if (allocated(this%m_iwork)) then
            liw = size(this%m_iwork)
        else
            liw = 0
        end if
        if (allocated(this%m_rwork)) then
            lrw = size(this%m_rwork)
        else
            lrw = 0
        end if

        ! Quick Return
        if (liw == 0 .or. lrw == 0) return

        ! Process
        this%m_rwork = 0.0d0
        this%m_iwork = 0
    end subroutine

! ------------------------------------------------------------------------------
    module subroutine oirk_init_workspace(this, liw, lrw, err)
        ! Arguments
        class(ode_irk), intent(inout) :: this
        integer(int32), intent(in) :: liw, lrw
        class(errors), intent(inout) :: err

        ! Local Variables
        integer(int32) :: flag

        ! Process
        flag = 0
        if (allocated(this%m_iwork)) then
            if (size(this%m_iwork) < liw) then
                deallocate(this%m_iwork)
                allocate(this%m_iwork(liw), stat = flag)
            end if
        else
            allocate(this%m_iwork(liw), stat = flag)
        end if

        if (flag == 0) then
            if (allocated(this%m_rwork)) then
                if (size(this%m_rwork) < lrw) then
                    deallocate(this%m_rwork)
                    allocate(this%m_rwork(lrw), stat = flag)
                end if
            else
                allocate(this%m_rwork(lrw), stat = flag)
            end if
        end if

        ! Error Checking
        if (flag /= 0) then
            call err%report_error("oirk_init_workspace", &
                "Insufficient memory available.", INT_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Initialization
        call this%reset()
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oirk_get_newton_iter_count(this) result(rst)
        class(ode_irk), intent(in) :: this
        integer(int32) :: rst
        rst = this%m_newtonIter
    end function

! --------------------
    module subroutine oirk_set_newton_iter_count(this, x)
        class(ode_irk), intent(inout) :: this
        integer(int32), intent(in) :: x
        this%m_newtonIter = x
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oirk_get_index1(this) result(rst)
        class(ode_irk), intent(in) :: this
        integer(int32) :: rst
        rst = this%m_index1
    end function

! --------------------
    module subroutine oirk_set_index1(this, x)
        class(ode_irk), intent(inout) :: this
        integer(int32), intent(in) :: x
        if (x < 0) then
            this%m_index1 = 0
        else
            this%m_index1 = x
        end if
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oirk_get_index2(this) result(rst)
        class(ode_irk), intent(in) :: this
        integer(int32) :: rst
        rst = this%m_index2
    end function

! --------------------
    module subroutine oirk_set_index2(this, x)
        class(ode_irk), intent(inout) :: this
        integer(int32), intent(in) :: x
        if (x < 0) then
            this%m_index2 = 0
        else
            this%m_index2 = x
        end if
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oirk_get_index3(this) result(rst)
        class(ode_irk), intent(in) :: this
        integer(int32) :: rst
        rst = this%m_index3
    end function

! --------------------
    module subroutine oirk_set_index3(this, x)
        class(ode_irk), intent(inout) :: this
        integer(int32), intent(in) :: x
        if (x < 0) then
            this%m_index3 = 0
        else
            this%m_index3 = x
        end if
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oirk_get_controller(this) result(rst)
        class(ode_irk), intent(in) :: this
        integer(int32) :: rst
        rst = this%m_controller
    end function

! --------------------
    module subroutine oirk_set_controller(this, x)
        class(ode_irk), intent(inout) :: this
        integer(int32), intent(in) :: x
        if (x /= INT_PREDICTIVE_STEP_SIZE_CONTROLLER .or. &
            x /= INT_CLASSICAL_STEP_SIZE_CONTROLLER) &
        then
            this%m_controller = INT_PREDICTIVE_STEP_SIZE_CONTROLLER
        else
            this%m_controller = x
        end if
    end subroutine

! ------------------------------------------------------------------------------
    pure module function oirk_get_newton_tol(this) result(rst)
        class(ode_irk), intent(in) :: this
        real(real64) :: rst
        rst = this%m_newtonTol
    end function

! --------------------
    module subroutine oirk_set_newton_tol(this, x)
        class(ode_irk), intent(inout) :: this
        real(real64), intent(in) :: x
        this%m_newtonTol = x
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

! ------------------------------------------------------------------------------
end submodule
