# Link in dependencies necessary to generate plots
find_package(fplot)
find_package(ferror 1.3.0)
find_package(fcore)

# Van Der Pol Example
add_executable(vanderpol vanderpol.f90)
target_link_libraries(vanderpol integral fplot fcore)

# Lorenz Example
add_executable(lorenz lorenz.f90)
target_link_libraries(lorenz integral fplot fcore)

# Bouncing Ball Example
add_executable(bouncing_ball bouncing_ball.f90)
target_link_libraries(bouncing_ball integral fplot fcore)

# Integral Example 1
add_executable(integral_example_1 integral_example_1.f90)
target_link_libraries(integral_example_1 integral)

# Pendulum Example
add_executable(pendulum pendulum.f90)
target_link_libraries(pendulum integral fplot fcore)

# Infinite Limit Integration Example 1
add_executable(infinite_integral_example_1 infinite_integral_example_1.f90)
target_link_libraries(infinite_integral_example_1 integral)

# Infinite Limit Integration Example 2
add_executable(infinite_integral_example_2 infinite_integral_example_2.f90)
target_link_libraries(infinite_integral_example_2 integral)

# Runge-Kutta Example 1
add_executable(rk_example_1 rk_example_1.f90)
target_link_libraries(rk_example_1 integral fplot fcore)
