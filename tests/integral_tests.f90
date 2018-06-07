! integral_tests.f90

program tests
    use finite_interval_tests
    use ode_tests

    ! Local Variables
    logical :: test_result, overall
    overall = .true.

    ! Tests
    test_result = integral_test_1()
    if (.not.test_result) overall = .false.

    test_result = integral_test_2()
    if (.not.test_result) overall = .false.

    test_result = ode_step_test_1()
    if (.not.test_result) overall = .false.

    test_result = ode_step_test_2()
    if (.not.test_result) overall = .false.

    test_result = ode_test_1()
    if (.not.test_result) overall = .false.

    test_result = ode_test_2()
    if (.not.test_result) overall = .false.

    ! Output
    if (overall) then
        print '(A)', "INTEGRAL TEST STATUS: PASS"
        call exit(0)
    else
        print '(A)', "FERROR TEST STATUS: FAIL"
        call exit(1)
    end if
end program
