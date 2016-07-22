IDIR =include
CC=g++
CFLAGS=-I$(IDIR)

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
_OBJ = wipp_main.o wipp_fileutils.o 
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

# Rules for making individual objects
$(ODIR)/%.o: $(SRC_DIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

# Rule to link everything together + generate executable
wipp: $(OBJ)
	$(CC) -o $(BDIR)/$@ $^ $(CFLAGS)

# Safety rule for any file named "clean"
.PHONY: clean

# Purge the build and bin directories
clean:
	rm -f $(ODIR)/*
	rm -f $(BDIR)/*