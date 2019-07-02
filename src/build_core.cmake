function(build_core CORE)
  set(EXE_NAME ${CORE}_model)
  set(NAMELIST_SUFFIX ${CORE})

  # Map the ESM component corresponding to each MPAS core
  if (CORE STREQUAL "ocean")
    set(COMPONENT "ocn")
  elseif(CORE STREQUAL "landice")
    set(COMPONENT "glc")
  elseif(CORE STREQUAL "seaice")
    set(COMPONENT "ice")
  endif()

  # build_options.mk stuff handled here
  if (CORE STREQUAL "ocean")
    list(APPEND CPPDEFS "-DCORE_OCEAN")
    list(APPEND INCLUDES "${CMAKE_BINARY_DIR}/core_ocean/BGC" "${CMAKE_BINARY_DIR}/core_ocean/shared" "${CMAKE_BINARY_DIR}/core_ocean/analysis_members" "${CMAKE_BINARY_DIR}/core_ocean/cvmix" "${CMAKE_BINARY_DIR}/core_ocean/mode_forward" "${CMAKE_BINARY_DIR}/core_ocean/mode_analysis" "${CMAKE_BINARY_DIR}/core_ocean/mode_init")

  elseif (CORE STREQUAL "seaice")
    list(APPEND CPPDEFS "-DCORE_SEAICE" "-Dcoupled" "-DCCSMCOUPLED")
    list(APPEND INCLUDES "${CMAKE_BINARY_DIR}/core_seaice/column" "${CMAKE_BINARY_DIR}/core_seaice/shared" "${CMAKE_BINARY_DIR}/core_seaice/analysis_members" "${CMAKE_BINARY_DIR}/core_seaice/model_forward")

  elseif (CORE STREQUAL "landice")
    list(APPEND CPPDEFS "-DCORE_LANDICE")
    list(APPEND INCLUDES "${CMAKE_BINARY_DIR}/core_landice/shared" "${CMAKE_BINARY_DIR}/core_landice/analysis_members" "${CMAKE_BINARY_DIR}/core_landice/mode_forward")

    #
    # Check if building with LifeV, Albany, and/or PHG external libraries
    #

    if (LIFEV)
      # LifeV can solve L1L2 or FO
      list(APPEND CPPDEFS "-DLIFEV" "-DUSE_EXTERNAL_L1L2" "-DUSE_EXTERNAL_FIRSTORDER" "-DMPAS_LI_BUILD_INTERFACE")
    endif()

    # Albany can only solve FO at present
    if (ALBANY)
      list(APPEND CPPDEFS "-DUSE_EXTERNAL_FIRSTORDER" "-DMPAS_LI_BUILD_INTERFACE")
    endif()

    if (LIFEV AND ALBANY)
      message(FATAL "Compiling with both LifeV and Albany is not allowed at this time.")
    endif()

    # PHG currently requires LifeV
    if (PHG AND NOT LIFEV)
      message(FATAL "Compiling with PHG requires LifeV at this time.")
    endif()

    # PHG can only Stokes at present
    if (PHG)
      list(APPEND CPPDEFS "-DUSE_EXTERNAL_STOKES" "-DMPAS_LI_BUILD_INTERFACE")
    endif()
  endif()

  add_library(${COMPONENT})
  target_compile_definitions(${COMPONENT} PRIVATE ${CPPDEFS})
  target_include_directories(${COMPONENT} PRIVATE ${INCLUDES})

  # Gather sources

  # externals
  set(RAW_SOURCES external/ezxml/ezxml.c)

  # framework
  list(APPEND RAW_SOURCES
    framework/mpas_kind_types.F
    framework/mpas_framework.F
    framework/mpas_timer.F
    framework/mpas_timekeeping.F
    framework/mpas_constants.F
    framework/mpas_attlist.F
    framework/mpas_hash.F
    framework/mpas_sort.F
    framework/mpas_block_decomp.F
    framework/mpas_block_creator.F
    framework/mpas_dmpar.F
    framework/mpas_abort.F
    framework/mpas_decomp.F
    framework/mpas_threading.F
    framework/mpas_io.F
    framework/mpas_io_streams.F
    framework/mpas_bootstrapping.F
    framework/mpas_io_units.F
    framework/mpas_stream_manager.F
    framework/mpas_stream_list.F
    framework/mpas_forcing.F
    framework/mpas_c_interfacing.F
    framework/random_id.c
    framework/pool_hash.c
    framework/mpas_derived_types.F
    framework/mpas_domain_routines.F
    framework/mpas_field_routines.F
    framework/mpas_pool_routines.F
    framework/xml_stream_parser.c
    framework/regex_matching.c
    framework/mpas_field_accessor.F
    framework/mpas_log.F
  )

  # operators
  list(APPEND RAW_SOURCES
    operators/mpas_vector_operations.F
    operators/mpas_matrix_operations.F
    operators/mpas_tensor_operations.F
    operators/mpas_rbf_interpolation.F
    operators/mpas_vector_reconstruction.F
    operators/mpas_spline_interpolation.F
    operators/mpas_tracer_advection_helpers.F
    operators/mpas_tracer_advection_mono.F
    operators/mpas_tracer_advection_std.F
    operators/mpas_geometry_utils.F
  )

  set(COMMON_RAW_SOURCES ${RAW_SOURCES})

  set(CORE_BLDDIR ${CMAKE_BINARY_DIR}/core_${CORE})
  if (NOT EXISTS ${CORE_BLDDIR})
    file(MAKE_DIRECTORY ${CORE_BLDDIR})
  endif()

  set(CORE_INPUT_DIR ${CORE_BLDDIR}/default_inputs)
  if (NOT EXISTS ${CORE_INPUT_DIR})
    file(MAKE_DIRECTORY ${CORE_INPUT_DIR})
  endif()

  # Make .inc files
  add_custom_command (
    OUTPUT ${CORE_BLDDIR}/Registry_processed.xml
    COMMAND cpp -P -traditional ${CPPDEFS} -Uvector
    ${CMAKE_CURRENT_SOURCE_DIR}/core_${CORE}/Registry.xml > Registry_processed.xml
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/core_${CORE}/Registry.xml
    WORKING_DIRECTORY ${CORE_BLDDIR}
    )

  set(INC_DIR ${CORE_BLDDIR}/inc)
  if (NOT EXISTS ${INC_DIR})
    file(MAKE_DIRECTORY ${INC_DIR})
  endif()

  add_custom_command(
    OUTPUT ${INC_DIR}/core_variables.inc
    COMMAND ${CMAKE_BINARY_DIR}/mpas-source/src/parse < ${CORE_BLDDIR}/Registry_processed.xml
    DEPENDS parse ${CORE_BLDDIR}/Registry_processed.xml
    WORKING_DIRECTORY ${INC_DIR}
  )

  include(${CMAKE_CURRENT_SOURCE_DIR}/core_${CORE}/${CORE}.cmake)

  # Disable qsmp for some files
  if (FFLAGS MATCHES ".*-qsmp.*")
    foreach(DISABLE_QSMP_FILE IN LISTS DISABLE_QSMP)
      get_filename_component(SOURCE_EXT ${DISABLE_QSMP_FILE} EXT)
      string(REPLACE "${SOURCE_EXT}" ".f90" SOURCE_F90 ${DISABLE_QSMP_FILE})
      set_property(SOURCE ${CMAKE_BINARY_DIR}/${SOURCE_F90} APPEND_STRING PROPERTY COMPILE_FLAGS " -nosmp")
    endforeach()
  endif()

  # Run all .F files through cpp to generate the f90 file
  foreach(ITEM IN LISTS INCLUDES)
    list(APPEND INCLUDES_I "-I${ITEM}")
  endforeach()

  list(GET CORES 0 FIRST_CORE)
  foreach(RAW_SOURCE_FILE IN LISTS RAW_SOURCES)
    get_filename_component(SOURCE_EXT ${RAW_SOURCE_FILE} EXT)
    if ( (SOURCE_EXT STREQUAL ".F" OR SOURCE_EXT STREQUAL ".F90") AND NOT RAW_SOURCE_FILE IN_LIST NO_PREPROCESS)
      string(REPLACE "${SOURCE_EXT}" ".f90" SOURCE_F90 ${RAW_SOURCE_FILE})
      get_filename_component(DIR_RELATIVE ${SOURCE_F90} DIRECTORY)
      set(DIR_ABSOLUTE ${CMAKE_BINARY_DIR}/${DIR_RELATIVE})
      if (NOT EXISTS ${DIR_ABSOLUTE})
        file(MAKE_DIRECTORY ${DIR_ABSOLUTE})
      endif()
      if (CORE STREQUAL ${FIRST_CORE} OR NOT RAW_SOURCE_FILE IN_LIST COMMON_RAW_SOURCES)
        add_custom_command (
          OUTPUT ${CMAKE_BINARY_DIR}/${SOURCE_F90}
          COMMAND cpp -P -traditional ${CPPDEFS} ${INCLUDES_I} -Uvector
          ${CMAKE_CURRENT_SOURCE_DIR}/${RAW_SOURCE_FILE} > ${CMAKE_BINARY_DIR}/${SOURCE_F90}
          DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${RAW_SOURCE_FILE} ${INC_DIR}/core_variables.inc)
      endif()
      list(APPEND SOURCES ${CMAKE_BINARY_DIR}/${SOURCE_F90})
    else()
      list(APPEND SOURCES ${RAW_SOURCE_FILE})
    endif()
  endforeach()

  target_sources(${COMPONENT} PRIVATE ${SOURCES})

endfunction(build_core)
