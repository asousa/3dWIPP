IDIR =include
EIGEN_DIR=lib/eigen/
# MATLAB_DIR=/usr/local/MATLAB/R2012a/extern/include
CC=c++

CFLAGS=-I$(IDIR) -I$(EIGEN_DIR)
# MATLAB_FLAGS = -L /usr/local/MATLAB/R2012a/bin/glnxa64 -leng -lm -lmx -lmex -lmat -lut -Wl,-rpath=/usr/local/MATLAB/R2012a/bin/glnxa64

# compiled module directory
ODIR =build
# Libraries
LDIR =/shared/users/asousa/WIPP/3dWIPP/lib

	
	
# output binary directory
BDIR =bin
# source files here
SRC_DIR=src

# Dependencies (header files)
_DEPS = wipp.h consts.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))


# Objects to build
_OBJ = wipp_main.o wipp_fileutils.o math_utils.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


# XFORM = lib/xform_double
IRBEM = lib/irbem-code
# GEOPACK = lib/geopack
# Rules for making individual objects
$(ODIR)/%.o: $(SRC_DIR)/%.cpp $(DEPS)

	$(CC) -c -o $@ $< $(CFLAGS) -L$(LDIR) 

# Rule to link everything together + generate executable
wipp: $(OBJ) $(LDIR)/liboneradesp.a
# wipp: $(OBJ) $(LDIR)/libxformd.a $(LDIR)/liboneradesp.a
	# $(MAKE) -C $(XFORM)
	$(CC) $(CFLAGS) $(OBJ) -L $(LDIR) -loneradesp -lgfortran  -o $(BDIR)/$@

# $(ODIR)/sm2geo.o: $(SRC_DIR)/sm2geo.c

# 	gcc -c -o $@ $< $(CFLAGS)

$(LDIR)/libxformd.a:
	$(MAKE) -C $(XFORM)

$(LDIR)/liboneradesp.a:
	$(MAKE) -C $(IRBEM) OS=linux64 ENV=gfortran64 all
	mv $(IRBEM)/source/liboneradesp_linux_x86_64.a $(LDIR)/liboneradesp.a

# Safety rule for any file named "clean"
.PHONY: clean

# Purge the build and bin directories
clean:
	rm -f $(ODIR)/*
	rm -f $(BDIR)/*
	rm -f $(LDIR)/libxformd.a
	rm -f $(LDIR)/libirbem.so
	$(MAKE) -C $(XFORM) clean
	$(MAKE) -C $(IRBEM) clean