IDIR =include
EIGEN_DIR=lib/eigen/
MATLAB_DIR=/usr/local/MATLAB/R2012a/extern/include
CC=c++

CFLAGS=-I$(IDIR) -I$(EIGEN_DIR) -I$(MATLAB_DIR)
MATLAB_FLAGS = -L /usr/local/MATLAB/R2012a/bin/glnxa64 -leng -lm -lmx -lmex -lmat -lut -Wl,-rpath=/usr/local/MATLAB/R2012a/bin/glnxa64

# compiled module directory
ODIR =build
# Libraries
LDIR =lib

	
	
# output binary directory
BDIR =bin
# source files here
SRC_DIR=src

# Dependencies (header files)
_DEPS = wipp.h consts.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))


# Objects to build
_OBJ = wipp_main.o wipp_fileutils.o damping_ngo.o damping_foust.o math_utils.o \
	   kp_to_pp.o polyfit.o psd_model.o integrand.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


# XFORM = lib/xform_double
# Rules for making individual objects
$(ODIR)/%.o: $(SRC_DIR)/%.cpp $(DEPS)

	$(CC) -c -o $@ $< $(CFLAGS) -I$(EIGEN_DIR) -L$(LDIR) 

# Rule to link everything together + generate executable
wipp: $(OBJ) libxformd.a $(ODIR)/gauss_legendre.o
	# $(MAKE) -C $(XFORM)
	$(CC) $(CFLAGS) $(OBJ) $(ODIR)/gauss_legendre.o -L $(LDIR) -lxformd -lgfortran $(MATLAB_FLAGS) -o $(BDIR)/$@

$(ODIR)/gauss_legendre.o: $(SRC_DIR)/gauss_legendre.c $(IDIR)/gauss_legendre.h

	gcc -c -o $@ $< $(CFLAGS)

libxformd.a:
	$(MAKE) -C lib/xform_double


# Safety rule for any file named "clean"
.PHONY: clean

# Purge the build and bin directories
clean:
	rm -f $(ODIR)/*
	rm -f $(BDIR)/*
	$(MAKE) -C $(XFORM) clean