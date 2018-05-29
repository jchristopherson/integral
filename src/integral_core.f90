! integral_core.f90

module integral_core
    use, intrinsic :: iso_fortran_env, only : int32, real64
    implicit none

    interface
        !> @brief Defines a function of one variable to be integrated.
        !!
        !! @param[in] x The value of the independent variable at which to
        !!  evaluate the function.
        !! @return The value of the function at @p x.
        function integral_fcn(x) result(f)
            use, intrinsic :: iso_fortran_env, only : real64
            real(real64), intent(in) :: x
            real(real64) :: f
        end function
    end interface
end module
