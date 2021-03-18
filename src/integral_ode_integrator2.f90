! integral_ode_integrator2.f90

submodule (integral_core) integral_ode_integrator2
contains
! ------------------------------------------------------------------------------
    module function oi2_integrate(this, fcnobj, x, y, err) result(rst)
        ! Arguments
        class(ode_integrator2), intent(inout) :: this
        class(ode_helper), intent(inout) :: fcnobj
        real(real64), intent(in), dimension(:) :: x, y
        class(errors), intent(inout), optional, target :: err
        real(real64), allocatable, dimension(:,:) :: rst

        ! Local Variables
        logical :: brk
        integer(int32) :: i, n, neqn, ncols, nbuffer, flag
        real(real64) :: xi
        real(real64), allocatable, dimension(:) :: ytemp
        real(real64), allocatable, dimension(:,:) :: buffer
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Set up the error handling
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        
        ! Initialization
        n = size(x)
        neqn = fcnobj%get_equation_count()
        ncols = neqn + 1
        nbuffer = this%get_min_buffer_size()
        this%m_collector => store_data
        call this%initialize_integrator(fcnobj, x, y, errmgr)
        if (errmgr%has_error_occurred()) return

        ! Allocate buffer space
        allocate(ytemp(neqn), stat = flag)
        if (flag == 0 .and. this%get_provide_all_output()) &
            allocate(buffer(nbuffer, ncols), stat = flag)
        if (flag == 0 .and. .not.this%get_provide_all_output()) &
            allocate(rst(size(x), ncols), stat = flag)
        if (flag /= 0) then
            call errmgr%report_error("oi2_integrate", &
                "Insufficient memory available.", INT_OUT_OF_MEMORY_ERROR)
            return
        end if

        ! Integration Process:
        ! Only a single call is needed.  The subroutine store_data handles
        ! collecting the output from the integrator

        ! Store the initial conditions
        xi = x(1)
        ytemp = y
        
        ! Perform the integration
        if (this%get_provide_all_output()) then
            ! Store the initial conditions
            buffer(1,1) = xi
            buffer(1,2:ncols) = ytemp
            i = 1

            brk = this%step(fcnobj, xi, ytemp, x(n), this%m_rtol, this%m_atol, &
                errmgr)
            if (errmgr%has_error_occurred()) go to 100

            ! Trim the output buffer
        100 continue
            rst = buffer(1:i,:)
        else
            ! Only return output at the appropriate spots
            do i = 2, size(x)
                ! Store the initial conditions
                rst(1,1) = xi
                rst(1,2:ncols) = ytemp

                ! Take the integration step.  The solution is returned in 
                ! xi and ytemp - this way there's no need for an update of
                ! these terms
                brk = this%step(fcnobj, xi, ytemp, x(i), this%m_rtol, &
                    this%m_atol, errmgr)
                if (errmgr%has_error_occurred()) return
                rst(i,1) = xi
                rst(i,2:ncols) = ytemp
                if (brk) exit
            end do
        end if
        
    contains
        ! Routine for storing data when all output is requested.
        subroutine store_data(x_, y_)
            ! Arguments
            real(real64), intent(in) :: x_
            real(real64), intent(in), dimension(:) :: y_

            ! Local Variables
            integer(int32) :: nOld, nNew
            real(real64), allocatable, dimension(:,:) :: copy

            ! Store the results and index the buffer tracking variable
            i = i + 1
            if (size(buffer, 1) < i) then
                ! Update the size of the buffer
                if (allocated(copy)) deallocate(copy)
                copy = buffer
                nOld = size(buffer, 1)
                nNew = nOld + this%get_min_buffer_size()
                deallocate(buffer)
                allocate(buffer(nNew, ncols))
                buffer(1:nOld,:) = copy
            end if

            ! Store the results
            buffer(i,1) = x_
            buffer(i,2:ncols) = y_
        end subroutine
    end function

! ------------------------------------------------------------------------------
end submodule
