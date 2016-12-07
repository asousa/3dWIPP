IDIR =include
EIGEN_DIR=lib/eigen/
# MATLAB_DIR=/usr/local/MATLAB/R2012a/extern/include
CC=c++

CFLAGS=-I$(IDIR) -I$(EIGEN_DIR) -fopenmp
# MATLAB_FLAGS = -L /usr/local/MATLAB/R2012a/bin/glnxa64 -leng -lm -lmx -lmex -lmat -lut -Wl,-rpath=/usr/local/MATLAB/R2012a/bin/glnxa64

# compiled module directory
ODIR =build
# Libraries
LDIR =/shared/users/asousa/WIPP/3dWIPP/lib

LIBS=-lgfortran -lxformd -lgeopackd
	
	
# output binary directory
BDIR =bin
# source files here
SRC_DIR=src

# Dependencies (header files)
_DEPS = wipp.h consts.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))


# Objects to build
sources = \
	wipp_main.cpp \
	wipp_fileutils.cpp \
	math_utils.cpp \
	wipp_methods.cpp \
	wipp_legacy_methods.cpp \
	bmodel.cpp \
	coord_transforms.cpp \

_OBJ = ${sources:.cpp=.o}	
# _OBJ = wipp_main.o \
# 	   wipp_fileutils.o \
# 	   math_utils.o \
# 	   wipp_methods.o \
# 	   wipp_legacy_methods.o \
# 	   bmodel.o \
# 	   coord_transforms.o \

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_FLUX_OBJ = flux_main.o \
			math_utils.o \
			wipp_fileutils.o \
			coord_transforms.o
FLUX_OBJ = $(patsubst %,$(ODIR)/%,$(_FLUX_OBJ))

XFORM = $(LDIR)/xform_double
IRBEM = $(LDIR)/irbem-code
# GEOPACK = lib/geopack
# Rules for making individual objects (from .cpp files in src/)
$(ODIR)/%.o: $(SRC_DIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -L$(LDIR)

# Rule to link everything together + generate executable
wipp: $(OBJ) $(LDIR)/libxformd.a $(LDIR)/libgeopackd.a $(LDIR)/libfoust.a
	$(CC) $(CFLAGS) $(OBJ) -L $(LDIR) $(LIBS) -o $(BDIR)/$@

flux: $(FLUX_OBJ)
	$(CC) $(CFLAGS) $(FLUX_OBJ) -L $(LDIR) $(LIBS) -o $(BDIR)/$@

# $(ODIR)/sm2geo.o: $(SRC_DIR)/sm2geo.c

# 	gcc -c -o $@ $< $(CFLAGS)


$(LDIR)/libfoust.a: $(LDIR)/foust/*.f95 
# Fortran code from Forrest Foust
	gfortran -O3 -o $(LDIR)/foust/types.o -c $(LDIR)/foust/types.f95 
	gfortran -O3 -o $(LDIR)/foust/constants.o -c $(LDIR)/foust/constants.f95 
	gfortran -O3 -o $(LDIR)/foust/util.o -c $(LDIR)/foust/util.f95 
	gfortran -O3 -o $(LDIR)/foust/bmodel_dipole.o -c $(LDIR)/foust/bmodel_dipole.f95
	gfortran -O3 -o $(LDIR)/foust/bmodel.o -c $(LDIR)/foust/bmodel.f95 

	ar rc $(LDIR)/libfoust.a $(LDIR)/foust/types.o $(LDIR)/foust/constants.o $(LDIR)/foust/util.o $(LDIR)/foust/bmodel_dipole.o $(LDIR)/foust/bmodel.o

# Tyganenko's geopack library (IGRF and coordinate transforms)
$(LDIR)/libgeopackd.a:
	gfortran -O3 -o lib/geopack/geopack-2008_dp.o -c $(LDIR)/geopack/Geopack-2008_dp.for
	gfortran -O3 -o lib/geopack/TS05_aka_TS04.o -c $(LDIR)/geopack/TS05_aka_TS04.for
	ar rc $(LDIR)/libgeopackd.a $(LDIR)/geopack/Geopack-2008_dp.o  $(LDIR)/geopack/TS05_aka_TS04.o

# Legacy coordinate transforms, used in raytracer
$(LDIR)/libxformd.a:
	$(MAKE) -C $(XFORM)


# IRBEM utility library:
$(LDIR)/liboneradesp.a:
	$(MAKE) -C $(IRBEM) OS=linux64 ENV=gnu64 all
	mv $(IRBEM)/source/liboneradesp_linux_x86_64.a $(LDIR)/liboneradesp.a

# Safety rule for any file named "clean"
.PHONY: clean

# Purge the build and bin directories
clean:
	rm -f $(ODIR)/*
	rm -f $(BDIR)/*
	rm -f $(LDIR)/libxformd.a
	rm -f $(LDIR)/liboneradesp.a
	rm -f $(LDIR)/libgeopackd.a
	rm -f $(LDIR)/libfoust.a

	$(MAKE) -C $(XFORM) clean
	$(MAKE) -C $(IRBEM) clean