! integral_ode_rk45.f90

submodule (integral_core) integral_ode_rk45
contains
! ------------------------------------------------------------------------------
    module function ork45_step(this, fcn, x, y, xout, rtol, atol, err) result(brk)
        ! Arguments
        class(ode_rk45), intent(inout) :: this
        class(ode_helper), intent(inout) :: fcn
        real(real64), intent(inout) :: x
        real(real64), intent(inout), dimension(:) :: y
        real(real64), intent(in) :: xout
        real(real64), intent(in), dimension(:) :: rtol, atol
        class(errors), intent(inout), optional, target :: err
        logical :: brk

        ! Local Variables
        integer(int32) :: neq, lrw, liw, itol, iout, lrw, liw, nrdens, idid, &
            ipar(1)
        real(real64) :: rpar(1)
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
        nrdens = neq
        lrw = 8 * neq + 5 * nrdens + 21
        liw = nrdens + 21
        if (this%get_provide_all_output()) then
            iout = 2
        else
            iout = 0
        end if

        ! Set up the workspace arrays
        call this%init_workspace(liw, lrw, errmgr)
        if (errmgr%has_error_occurred()) return
        this%m_rwork(6) = this%get_max_step_size()
        this%m_iwork(1) = this%get_iteration_limit()
        this%m_iwork(5) = nrdens

        ! Compute the solution
        call DOPRI5(neq, doprifcn, x, y, xout, rtol, atol, itol, dopriout, &
            iout, this%m_rwork, lrw, this%m_iwork, liw, rpar, ipar, idid)
        select case (idid)
        case (-1)
            call errmgr%report_error("ork45_step", &
                "Invalid input encountered.", INT_INVALID_INPUT_ERROR)
            return
        case (-2)
            call errmgr%report_error("ork45_step", &
                "The step could not converge in the allowed number of " // &
                "iterations.  Try increasing the number of allowed " // &
                "iterations.", INT_CONVERGENCE_ERROR)
            return
        case (-3)
            call errmgr%report_error("ork45_step", &
                "The step size has become too small.", &
                INT_STEP_SIZE_TOO_SMALL_ERROR)
            return
        case (-4)
            call errmgr%report_error("ork45_step", &
                "The problem appears too stiff for this integrator.  " // &
                "Try an integrator better suited to stiff problems.", &
                INT_IMPROPER_INTEGRATOR_ERROR)
            return
        end select
    contains
        subroutine doprifcn(n_, x_, y_, f_, rpar_, ipar_)
            integer(int32), intent(in) :: n_
            real(real64), intent(in) :: x_, y_(n_)
            real(real64), intent(out) :: f_(n_)
            real(real64), intent(inout) :: rpar_(*)
            integer(int32), intent(inout) :: ipar_(*)
            call fcn%evaluate_ode(x_, y_, f_)
        end subroutine

        subroutine dopriout(nr_, xold_, x_, y_, n_, con_, icomp_, nd_, rpar_, &
                ipar_, irtrn_)
            integer(int32), intent(in) :: nr_, n_, nd_, icomp_(nd_)
            real(real64), intent(in) :: xold_, x_, y_(n_), con_(5*nd_)
            real(real64), intent(inout) :: rpar_(*)
            integer(int32), intent(inout) :: ipar_(*)
            integer(int32), intent(out) :: irtrn_ ! set to < 0 to force exit

            if (associated(this%m_collector)) call this%m_collector(x_, y_)
        end subroutine
    end function

! ------------------------------------------------------------------------------
    module subroutine ork45_reset_integrator(this)
        ! Arguments
        class(ode_rk45), intent(inout) :: this

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
    module subroutine ork45_init_workspace(this, liw, lrw, err)
        ! Arguments
        class(ode_rk45), intent(inout) :: this
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
            call err%report_error("ork45_init_workspace", &
                "Insufficient memory available.", INT_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Initialization
        call this%reset()
    end subroutine

! ------------------------------------------------------------------------------
end submodule
