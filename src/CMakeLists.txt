# Define sources
set(integral_sources
)

# Define the library
add_library(integral ${integral_sources})
target_link_libraries(integral ${FERROR_LIBRARIES})
