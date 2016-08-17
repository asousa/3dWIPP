IDIR =include
EIGEN_DIR=lib/eigen/
CC=c++

CFLAGS=-I$(IDIR) -I$(EIGEN_DIR)
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
_OBJ = wipp_main.o wipp_fileutils.o damping_ngo.o math_utils.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))


# XFORM = lib/xform_double
# Rules for making individual objects
$(ODIR)/%.o: $(SRC_DIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -I$(EIGEN_DIR) -L$(LDIR) 

# Rule to link everything together + generate executable
wipp: $(OBJ) libxformd.a
	# $(MAKE) -C $(XFORM)
	$(CC) $(CFLAGS) $(OBJ) -L $(LDIR) -lxformd -lgfortran -o $(BDIR)/$@

libxformd.a:
	$(MAKE) -C lib/xform_double


# Safety rule for any file named "clean"
.PHONY: clean

# Purge the build and bin directories
clean:
	rm -f $(ODIR)/*
	rm -f $(BDIR)/*
	$(MAKE) -C $(XFORM) clean