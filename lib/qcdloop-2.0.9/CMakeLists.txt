cmake_minimum_required (VERSION 3.0.2)

# Disable the use of RPATHS - we probably are not
# that interested in relocatable binaries and it
# causes a lot of trouble on macOS
SET(CMAKE_SKIP_RPATH TRUE)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

project(qcdloop)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(VERSION "\"2.0.9\"")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -fext-numeric-literals")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopenmp -fext-numeric-literals")
else()
    message("Unsupported C++ compiler")
endif()
set(CMAKE_EXE_LINKER_FLAGS_DEBUG "-fsanitize=address" CACHE STRING "debug linker flags" FORCE)
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix "${prefix}")
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")

configure_file(
  "${PROJECT_SOURCE_DIR}/src/qcdloop/config.h.in"
  "${PROJECT_SOURCE_DIR}/src/qcdloop/config.h"
  )

configure_file(
  "${PROJECT_SOURCE_DIR}/src/qcdloop-config.in"
  "${PROJECT_BINARY_DIR}/src/qcdloop-config"
  )

configure_file(
  "${PROJECT_SOURCE_DIR}/src/qcdloop.pc.in"
  "${PROJECT_SOURCE_DIR}/src/qcdloop.pc"
  )

set(QUADMATH_NAMES ${QUADMATH_NAMES} libquadmath.so quadmath)
find_library(QUADMATH_LIBRARY
  NAMES ${QUADMATH_NAMES}
  PATHS /usr/lib64/atlas /usr/lib/atlas
        /usr/lib64 /usr/lib /usr/local/lib64
        /usr/local/lib /usr/x86_64-linux-gnu/*
        /usr/lib/gcc/x86_64-linux-gnu/*
        /usr/lib/gcc/x86_64-redhat-linux/*
)

if(QUADMATH_LIBRARY)
  set(QUADMATH_LIBRARIES ${QUADMATH_LIBRARY})
  set(QUADMATH_FOUND "YES")
else(QUADMATH_LIBRARY)
  set(QUADMATH_FOUND "NO")
endif(QUADMATH_LIBRARY)

if(QUADMATH_FOUND)
  if(NOT QUADMATH_FIND_QUIETLY)
    message(STATUS "Found QUADMATH: ${QUADMATH_LIBRARIES}")
  endif(NOT QUADMATH_FIND_QUIETLY)
else(QUADMATH_FOUND)
  if(QUADMATH_FIND_REQUIRED)
    message(FATAL_ERROR "Could not find QuadMath")
  endif(QUADMATH_FIND_REQUIRED)
endif(QUADMATH_FOUND)

# libqcdloop configuration
include_directories(src/qcdloop src)
FILE(GLOB_RECURSE Headers "src/qcdloop/*.h")
add_library(objlib OBJECT src/box.cc
                         src/bubble.cc
                         src/cache.cc
                         src/qcdloop.cc
                         src/tadpole.cc
                         src/tools.cc
                         src/topology.cc
                         src/triangle.cc
                         src/types.cc
                         src/wrapper.cc
                         ${Headers}
)
set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1)

add_library(qcdloop_shared SHARED $<TARGET_OBJECTS:objlib>)
add_library(qcdloop_static STATIC $<TARGET_OBJECTS:objlib>)

set_target_properties(qcdloop_shared PROPERTIES OUTPUT_NAME qcdloop)
set_target_properties(qcdloop_static PROPERTIES OUTPUT_NAME qcdloop)

target_link_libraries(qcdloop_shared)
target_link_libraries(qcdloop_static)
#install(FILES ${CMAKE_BINARY_DIR}/src/qcdloop-config
  #DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ
  #GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
#install(FILES ${PROJECT_SOURCE_DIR}/src/qcdloop.pc DESTINATION lib/pkgconfig)
#install(DIRECTORY src/qcdloop DESTINATION include)
install(TARGETS qcdloop_shared DESTINATION lib)
install(TARGETS qcdloop_static DESTINATION lib)

# enable disable fortran/wrapper (to avoid name conflicts with ql1.x)
option(ENABLE_FORTRAN_WRAPPER "Enable fortran wrapper" ON)
if(ENABLE_FORTRAN_WRAPPER)
  add_definitions(-DQL_NAMES)
endif()

# qcdloop examples
option(ENABLE_EXAMPLES "Compile test examples" OFF)
if(ENABLE_EXAMPLES)
  add_executable(cache_test examples/cache_test.cc)
  target_link_libraries(cache_test qcdloop ${QUADMATH_LIBRARIES})
  set_target_properties(cache_test PROPERTIES LINK_FLAGS "-Wl,-rpath,./")

  add_executable(cmass_test examples/cmass_test.cc)
  target_link_libraries(cmass_test qcdloop ${QUADMATH_LIBRARIES})
  set_target_properties(cmass_test PROPERTIES LINK_FLAGS "-Wl,-rpath,./")

  add_executable(trigger_test examples/trigger_test.cc)
  target_link_libraries(trigger_test qcdloop ${QUADMATH_LIBRARIES})
  set_target_properties(trigger_test PROPERTIES LINK_FLAGS "-Wl,-rpath,./")
endif()
