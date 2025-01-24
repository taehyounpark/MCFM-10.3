cmake_minimum_required(VERSION 3.13)
enable_language(Fortran)
project(vvamp)

if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -finit-local-zero -ffixed-line-length-none -ffree-line-length-512 -fopenmp")
else()
    message( FATAL_ERROR "Unsupported Fortran compiler ${CMAKE_Fortran_COMPILER_ID}")
endif()

add_library(vvamp STATIC)

target_sources(vvamp PRIVATE
    mod_vvamp.f
    mod_vvamp_eval1.f
    mod_vvamp_eval2.f
    mod_vvamp_eval3.f
    mod_vvamp_eval4.f
    mod_vvamp_eval5.f
    mod_vvamp_eval6.f
    mod_vvamp_eval7.f
    mod_vvamp_eval8.f
    mod_vvamp_eval9.f
    mod_vvamp_eval10.f
    mod_vvamp_eval11.f
    mod_vvamp_eval12.f
)

