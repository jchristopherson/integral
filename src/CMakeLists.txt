# Define sources
set(integral_sources
    integral_core.f90
    integral_quadpack.f90
)

# Define the library
add_library(integral ${integral_sources})
target_link_libraries(integral ${FERROR_LIBRARIES})
