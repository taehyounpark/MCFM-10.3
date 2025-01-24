# John Campbell, April 2022
# Simple configuration script to use correct compiler
# and flags (passed by CMake as argument)
#!/bin/sh
sed -e "s:FC = gfortran:FC = $2:g" -e "s:FFLAGS = -O:FFLAGS = ${*:3} -fPIC:g" $1/oneloop/Config.orig > $1/oneloop/Config
