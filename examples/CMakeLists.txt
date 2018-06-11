# Link in dependencies necessary to generate plots
find_package(fplot)
find_package(ferror 1.3.0)
find_package(collections)

# Van Der Pol Example
add_executable(vanderpol vanderpol.f90)
target_link_libraries(vanderpol integral fplot)

# Lorenz Example
add_executable(lorenz lorenz.f90)
target_link_libraries(lorenz integral fplot)

# Bouncing Ball Example
add_executable(bouncing_ball bouncing_ball.f90)
target_link_libraries(bouncing_ball integral fplot)

# Integral Example 1
add_executable(integral_example_1 integral_example_1.f90)
target_link_libraries(integral_example_1 integral)

# Pendulum Example
add_executable(pendulum pendulum.f90)
target_link_libraries(pendulum integral fplot)

# ------------------------------------------------------------------------------
if (WIN32)
    # Copy shared libraries necessary to generate plots to the runtime directory
    get_target_property(ferror_LibLocation ferror LOCATION)
    add_custom_command(TARGET vanderpol POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${ferror_LibLocation} $<TARGET_FILE_DIR:vanderpol>
    )

    get_target_property(fplot_LibLocation fplot LOCATION)
    add_custom_command(TARGET vanderpol POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${fplot_LibLocation} $<TARGET_FILE_DIR:vanderpol>
    )

    get_target_property(collections_LibLocation collections LOCATION)
    add_custom_command(TARGET vanderpol POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${collections_LibLocation} $<TARGET_FILE_DIR:vanderpol>
    )
endif()
