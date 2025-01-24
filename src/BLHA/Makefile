# Makefile for MCFM C++ interface.
# Copyright (C) 2021 John M Campbell, Stefan Hoeche, Christian T Preuss

PREFIX = $(PWD)/../..
MCFMHOME = $(realpath $(PREFIX))
UNAME := $(shell uname -s)

CXXFLAGS=-g -O2 -fPIC -DPIC -fopenmp
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname,
LIB_SUFFIX=.so
OMP_FLAG=-fopenmp
ifeq ($(UNAME),Darwin)
  CXX_SHARED=-dynamiclib
  CXX_SONAME=-Wl,-dylib_install_name,
  LIB_SUFFIX=.dylib
endif

QLLIBS=$(PREFIX)/qcdloop-2.0.5/local/lib
QDLIBS=$(PREFIX)/qd-2.3.22/local/lib
LIBFLAGS=-lqcdloop -lqdmod -lqd -lquadmath

OBJECTS=Flavor_Map.o CXX_Interface.o \
	$(filter-out $(PREFIX)/obj/lhapdf_fortran.o, \
	  $(wildcard $(PREFIX)/obj/*.o)) \
	$(wildcard $(PREFIX)/obj/*.cpo) \
	$(filter-out $(PREFIX)/TensorReduction/pv/obj/pvfndd.o, \
	  $(wildcard $(PREFIX)/TensorReduction/*/obj/*.o)) \
	$(wildcard $(PREFIX)/TensorReduction/*/*/obj/*.o) 
HEADERS=CXX_Interface.h CXX_Wrapper.h Flavor_Map.h

# OpenLoops.
OLDIR=$(realpath $(PREFIX)/../OpenLoops)
OLLIBS=$(OLDIR)/lib
OLINCLUDE=$(OLDIR)/include

# Recola.
RCLDIR=$(realpath $(PREFIX)/../Recola2)
RCLLIBS=$(RCLDIR)/lib
RCLINCLUDE=$(RCLDIR)/include

# Recola.
RCLHEFTDIR=$(subst /install,_heft/install,$(RCLDIR))
RCLHEFTLIBS=$(RCLHEFTDIR)/lib
RCLHEFTINCLUDE=$(RCLHEFTDIR)/include

# MadLoop.
ML5DIR=$(realpath $(PREFIX)/../MadGraph)

all: install test

CXX_Interface.o: $(HEADERS) $(filter-out Flavor_Map.cxx,$(patsubst .cxx,.o,$(wildcard *.cxx)))
	$(CXX) $(CXXFLAGS) -I$(PREFIX)/include/ -c $(basename $@).cxx

%.o: $(HEADERS) %.cxx
	$(CXX) $(CXXFLAGS) -I$(PREFIX)/include/ -c $(basename $@).cxx

libMCFM$(LIB_SUFFIX): $(OBJECTS)
	@$(CXX) $(CXX_SHARED) $(OMP_FLAG) $(CXX_SONAME)$(MCFMHOME)/lib/$@ -o $@ $(OBJECTS) \
	-L$(QLLIBS) -L$(QDLIBS) $(LIBFLAGS) -lgfortran

install-headers: $(HEADERS) 
	test -d $(PREFIX)/include/MCFM || mkdir -p $(PREFIX)/include/MCFM/
	rsync -u $(HEADERS) $(PREFIX)/include/MCFM/

install: install-headers libMCFM$(LIB_SUFFIX)
	test -d $(PREFIX)/lib || mkdir -p $(PREFIX)/share/ $(PREFIX)/lib/
	rsync -u libMCFM$(LIB_SUFFIX) $(PREFIX)/lib/
	rsync -u $(PREFIX)/Bin/process.DAT $(PREFIX)/share/

clean:
	rm -f *.o test oltest libMCFM$(LIB_SUFFIX)
	rm -rf $(PREFIX)/share/ $(PREFIX)/lib/ $(PREFIX)/include/MCFM/ $(PWD)/lib

test.o: test.cxx params.cxx libMCFM$(LIB_SUFFIX)
	$(CXX) $(CXXFLAGS) -I$(PREFIX)/include -c $<

test: test.o
	$(CXX) $(CXXFLAGS) -o test $< -L$(PREFIX)/lib -lMCFM
	rsync -u $(PREFIX)/Bin/process.DAT .

oltest.o: oltest.cxx params.cxx install
	$(CXX) $(CXXFLAGS) -I$(PREFIX)/include -I$(OLINCLUDE) -c $<

oltest: oltest.o
	mkdir -p $(PWD)/lib
	rsync -u $(OLLIBS)/*$(LIB_SUFFIX) $(PWD)/lib/
	$(CXX) $(CXXFLAGS) -o $@ $< -L$(PREFIX)/lib -Wl,-rpath $(PREFIX)/lib -lMCFM \
	-L$(OLLIBS) -Wl,-rpath $(OLLIBS) -lopenloops
	rsync -u $(PREFIX)/Bin/process.DAT .

rcltest.o: rcltest.cxx params.cxx install
	$(CXX) $(CXXFLAGS) -I$(PREFIX)/include -I$(RCLINCLUDE) -c $<

rcltest: rcltest.o
	$(CXX) $(CXXFLAGS) -o $@ $< -L$(PREFIX)/lib -Wl,-rpath $(PREFIX)/lib -lMCFM \
	-L$(RCLLIBS) -Wl,-rpath $(RCLLIBS) -lrecola
	rsync -u $(PREFIX)/Bin/process.DAT .

rcltest_heft.o: rcltest.cxx params.cxx install
	$(CXX) $(CXXFLAGS) -I$(PREFIX)/include -I$(RCLHEFTINCLUDE) -o $@ -c $<

rcltest_heft: rcltest_heft.o
	$(CXX) $(CXXFLAGS) -o $@ $< -L$(PREFIX)/lib -Wl,-rpath $(PREFIX)/lib -lMCFM \
	-L$(RCLHEFTLIBS) -Wl,-rpath $(RCLHEFTLIBS) -lrecola
	rsync -u $(PREFIX)/Bin/process.DAT .

mltest.o: mltest.cxx params.cxx install
	$(CXX) $(CXXFLAGS) -DLIB_SUFFIX=\"$(LIB_SUFFIX)\" -I$(PREFIX)/include -c $<

mltest: mltest.o
	$(CXX) $(CXXFLAGS) -o $@ $< -L$(PREFIX)/lib -Wl,-rpath $(PREFIX)/lib -lMCFM -ldl
	sed -e's|ML5DIR|$(ML5DIR)|1' makeloops.in > makeloops && chmod ugo+x makeloops
	rsync -u $(PREFIX)/Bin/process.DAT .
