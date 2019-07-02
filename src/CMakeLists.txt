include(${CASEROOT}/Macros.cmake)

#
# General setup
#

if (USE_ESMF_LIB)
  set(ESMFDIR "esmf")
else()
  set(ESMFDIR "noesmf")
endif()

set(CMAKE_C_COMPILER ${MPICC})
set(CMAKE_CXX_COMPILER ${MPICXX})
set(CMAKE_Fortran_COMPILER ${MPIFC})
set(CMAKE_EXE_LINKER_FLAGS "${LDFLAGS}")
set(CMAKE_VERBOSE_MAKEFILE TRUE)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../../lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../../lib)

# Set up CPPDEFS
set(FILE_OFFSET "-DOFFSET64BIT")
set(CPPDEFS "-DMPAS_NO_LOG_REDIRECT" "-DMPAS_NO_ESMF_INIT" "-DMPAS_ESM_SHR_CONST" "-DMPAS_PERF_MOD_TIMERS" "${MODEL_FORMULATION}" "${FILE_OFFSET}" "${ZOLTAN_DEFINE}" "-D_MPI" "-DMPAS_NAMELIST_SUFFIX=${NAMELIST_SUFFIX}" "-DMPAS_EXE_NAME=${EXE_NAME}")
if (DEBUG)
  list(APPEND CPPDEFS "-DMPAS_DEBUG")
endif()
if (compile_threaded)
  list(APPEND CPPDEFS "-DMPAS_OPENMP")
endif()

set(INCLUDES "${INSTALL_SHAREDPATH}/include" "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/csm_share" "${NETCDF_PATH}/include" "${INSTALL_SHAREDPATH}/pio" "${PNETCDF_PATH}/include" "${CMAKE_CURRENT_SOURCE_DIR}/external/ezxml" "${CMAKE_BINARY_DIR}/framework" "${CMAKE_BINARY_DIR}/operators")

set(CMAKE_Fortran_FLAGS "${FFLAGS}")
set(CMAKE_C_FLAGS "${CFLAGS}")
set(CMAKE_CXX_FLAGS "${CXXFLAGS}")

# Make build tools
set(CMAKE_C_COMPILER ${SCC})

add_executable(streams_gen tools/input_gen/streams_gen.c tools/input_gen/test_functions.c external/ezxml/ezxml.c)
add_executable(namelist_gen tools/input_gen/namelist_gen.c tools/input_gen/test_functions.c external/ezxml/ezxml.c)
add_executable(parse tools/registry/parse.c tools/registry/dictionary.c tools/registry/gen_inc.c tools/registry/fortprintf.c tools/registry/utility.c external/ezxml/ezxml.c)

foreach(EXEITEM streams_gen namelist_gen parse)
  target_compile_definitions(${EXEITEM} PRIVATE ${CPPDEFS})
  target_compile_options(${EXEITEM} PRIVATE "-Uvector")
  target_include_directories(${EXEITEM} PRIVATE ${INCLUDES})
endforeach()

set(CMAKE_C_COMPILER ${MPICC})

# Macro for handling nl and st gen
macro(handle_nl_st_gen)
  foreach(NL_GEN_ARG IN LISTS NL_GEN_ARGS)
    separate_arguments(ITEMS UNIX_COMMAND "${NL_GEN_ARG}")
    list(GET ITEMS 0 ITEM)
    list(APPEND INPUTS ${ITEM})
    add_custom_command(
      OUTPUT ${CORE_INPUT_DIR}/${ITEM}
      COMMAND ${CMAKE_BINARY_DIR}/namelist_gen ${CORE_BLDDIR}/Registry_processed.xml ${NL_GEN_ARG}
      DEPENDS namelist_gen ${CORE_BLDDIR}/Registry_processed.xml
      WORKING_DIRECTORY ${CORE_INPUT_DIR}
    )
  endforeach()

  foreach(ST_GEN_ARG IN LISTS ST_GEN_ARGS)
    separate_arguments(ITEMS UNIX_COMMAND "${ST_GEN_ARG}")
    list(GET ITEMS 0 ITEM)
    list(APPEND INPUTS ${ITEM})
    add_custom_command(
      OUTPUT ${CORE_INPUT_DIR}/${ITEM}
      COMMAND ${CMAKE_BINARY_DIR}/streams_gen ${CORE_BLDDIR}/Registry_processed.xml ${ST_GEN_ARG}
      DEPENDS streams_gen ${CORE_BLDDIR}/Registry_processed.xml
      WORKING_DIRECTORY ${CORE_INPUT_DIR}
    )
  endforeach()

  foreach(INPUT IN LISTS INPUTS)
    add_custom_command(
      OUTPUT ${CORE_BLDDIR}/${INPUT}
      COMMAND ${CMAKE_COMMAND} -E copy ${CORE_INPUT_DIR}/${INPUT} ${CORE_BLDDIR}/${INPUT}
      DEPENDS ${CORE_INPUT_DIR}/${INPUT}
      WORKING_DIRECTORY ${CORE_BLDDIR}
    )
  endforeach()
endmacro(handle_nl_st_gen)

include(${CMAKE_CURRENT_SOURCE_DIR}/build_core.cmake)

foreach(CORE IN LISTS CORES)
  build_core(${CORE})
endforeach()
