# Define sources
set(integral_sources
    integral_core.f90
    integral_integrator_base.f90
    integral_adaptive_integrator.f90
    integral_nonadaptive_integrator.f90
    integral_ode_helper.f90
    integral_ode_integrator.f90
    integral_ode_auto.f90
)

# Define the QUADPACK source files
# REF:
# - http://www.netlib.org/quadpack/
# - https://people.sc.fsu.edu/~jburkardt/f_src/quadpack/quadpack.html
set(quadpack_sources
    QUADPACK/dqagie.f
    QUADPACK/dqags.f
    QUADPACK/dqagse.f
    QUADPACK/dqawc.f
    QUADPACK/dqawce.f
    QUADPACK/dqawf.f
    QUADPACK/dqawfe.f
    QUADPACK/dqawo.f
    QUADPACK/dqawoe.f
    QUADPACK/dqaws.f
    QUADPACK/dqawse.f
    QUADPACK/dqc25c.f
    QUADPACK/dqc25f.f
    QUADPACK/dqc25s.f
    QUADPACK/dqcheb.f
    QUADPACK/dqelg.f
    QUADPACK/dqk15.f
    QUADPACK/dqk15i.f
    QUADPACK/dqk15w.f
    QUADPACK/dqk21.f
    QUADPACK/dqk31.f
    QUADPACK/dqk41.f
    QUADPACK/dqk51.f
    QUADPACK/dqk61.f
    QUADPACK/dqmomo.f
    QUADPACK/dqng.f
    QUADPACK/dqpsrt.f
    QUADPACK/dqwgtc.f
    QUADPACK/dqwgtc.f
    QUADPACK/dqwgtf.f
    QUADPACK/dqwgts.f
    QUADPACK/dqag.f
    QUADPACK/dqage.f
    QUADPACK/dqagi.f
    QUADPACK/dqagp.f
    QUADPACK/dqagpe.f
    QUADPACK/xerror.f90
    QUADPACK/d1mach.f
    QUADPACK/dgtsl.f
)

# Include the QUADPACK sources
list(APPEND integral_sources ${quadpack_sources})

# Define the ODEPACK sources
set(odepack_sources
    ODEPACK/opkda1.f
    ODEPACK/opkda2.f
    ODEPACK/opkdmain.f
)

# Include the ODEPACK sources
list(APPEND integral_sources ${odepack_sources})

# Define the library
add_library(integral ${integral_sources})
target_link_libraries(integral ${FERROR_LIBRARIES})

# ------------------------------------------------------------------------------
# INSTALLATION INSTRUCTIONS
# ------------------------------------------------------------------------------
# Define target information
set_property(TARGET integral PROPERTY VERSION ${INTEGRAL_VERSION})
set_property(TARGET integral PROPERTY SOVERSION ${INTEGRAL_MAJOR_VERSION})
set_property(TARGET integral PROPERTY INTERFACE_integral_MAJOR_VERSION ${INTEGRAL_MAJOR_VERSION})
set_property(TARGET integral APPEND PROPERTY COMPATIBLE_INTERFACE_STRING integral_MAJOR_VERSION)

# Locate the "include" directories
set(integral_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include)

# Define the installation instructions
install(TARGETS integral EXPORT integralTargets
    RUNTIME DESTINATION integral/bin
    LIBRARY DESTINATION integral/lib
    ARCHIVE DESTINATION integral/lib
    CONFIGURATIONS Release
    INCLUDES DESTINATION integral/include
)
install(
    DIRECTORY ${PROJECT_SOURCE_DIR}/include
    DESTINATION ${CMAKE_INSTALL_PREFIX}/integral
    COMPONENT Devel
)

# Include the documentation
install(DIRECTORY ${PROJECT_SOURCE_DIR}/doc/html DESTINATION integral/doc)

# Define the version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/integralConfigVersion.cmake"
    VERSION ${INTEGRAL_VERSION}
    COMPATIBILITY AnyNewerVersion
)

export(EXPORT integralTargets
    FILE "${CMAKE_BINARY_DIR}/integralTargets.cmake"
)

# Define the configuration file
configure_file(
    "${PROJECT_SOURCE_DIR}/integralConfig.cmake.in"
    "${CMAKE_BINARY_DIR}/integralConfig.cmake"
    COPYONLY
)

set(ConfigPackageLocation integral/lib/cmake/integral)
install(
    EXPORT integralTargets
    FILE integralTargets.cmake
    DESTINATION ${ConfigPackageLocation}
)
install(
    FILES
        "${CMAKE_BINARY_DIR}/integralConfig.cmake"
        "${CMAKE_BINARY_DIR}/integralConfigVersion.cmake"
    DESTINATION
        ${ConfigPackageLocation}
    COMPONENT
        Devel
)
