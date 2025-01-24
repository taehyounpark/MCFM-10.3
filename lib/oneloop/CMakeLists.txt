cmake_minimum_required(VERSION 3.13)
enable_language(Fortran)
project(oneloop)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -finit-local-zero -ffixed-line-length-none -ffree-line-length-512 -fopenmp")
else()
    message( FATAL_ERROR "Unsupported Fortran compiler ${CMAKE_Fortran_COMPILER_ID}")
endif()

add_library(oneloop STATIC)

target_sources(oneloop PRIVATE
    avh_olo.f90
)
