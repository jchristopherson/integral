# Define sources
set(integral_sources
    integral_core.f90
    integral_integrator_base.f90
    integral_adaptive_integrator.f90
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
