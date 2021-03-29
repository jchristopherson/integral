# Define the source files
set(integral_test_sources
    integral_tests.f90
    finite_interval_tests.f90
    ode_tests.f90
    infinite_interval_tests.f90
)

# Build the executable
add_executable(integral_test ${integral_test_sources})
target_link_libraries(integral_test integral)

# Run the tests
add_custom_command(
    OUTPUT integral_api_test
    DEPENDS ${integral_test_sources}
    COMMAND integral_test
)
add_custom_target(run_integral_tests ALL DEPENDS integral_api_test)

# Be sure to copy FERROR on Windows systems
if (WIN32)
    get_target_property(ferror_LibLocation ferror LOCATION)
    add_custom_command(TARGET integral_test POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${ferror_LibLocation} $<TARGET_FILE_DIR:integral_test>
    )
endif()
