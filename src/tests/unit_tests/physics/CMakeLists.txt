SET(TESTS
    utest_CH_slider_pend
    utest_CH_double_pend
    utest_CH_shafts
    utest_CH_compute_contact
    utest_CH_assembly
    utest_CH_composite_inertia
    utest_CH_prismatic_spring
    utest_CH_actuators
    utest_CH_actuator_models
)

if(CH_ENABLE_MODULE_PARDISO_MKL)
  set(TESTS ${TESTS}
    utest_CH_gears
  )
endif()

set(LIBS Chrono_core)
if(CH_ENABLE_MODULE_PARDISO_MKL)
    list(APPEND LIBS Chrono_pardisomkl)
endif()


MESSAGE(STATUS "Add unit test programs for PHYSICS module")
build_utests(YES "${TESTS}" "physics" "${LIBS}")

add_executable(utest_CH_actuator_models_isolated utest_CH_actuator_models_isolated.cpp)
set_target_properties(utest_CH_actuator_models_isolated PROPERTIES FOLDER tests)
set_property(TARGET utest_CH_actuator_models_isolated PROPERTY
             INCLUDE_DIRECTORIES "${CMAKE_SOURCE_DIR}/src/chrono/physics/actuator_models")
add_test(utest_CH_actuator_models_isolated ${PROJECT_BINARY_DIR}/bin/utest_CH_actuator_models_isolated)
