# Compiler and flags
CC = gcc
FC = gfortran

# Detect system and set paths dynamically
UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    GSL_PREFIX ?= /usr/local
    MINPACK_DIR ?= $(shell find /usr/local -name "*minpack*" 2>/dev/null | head -1)
    ifeq ($(MINPACK_DIR),)
        MINPACK_DIR = $(shell pwd)/../lib
    endif
    # Look for Fortran module files
    MINPACK_MOD_DIR ?= /Users/School/Documents/Fortran/MINPACK
    ifeq ($(wildcard $(MINPACK_MOD_DIR)/*.mod),)
        MINPACK_MOD_DIR = $(MINPACK_DIR)
    endif
endif

CFLAGS = -Wall -fPIC -I$(GSL_PREFIX)/include
FFLAGS = -fPIC -I$(MINPACK_MOD_DIR)
LDFLAGS = -L$(GSL_PREFIX)/lib -L$(MINPACK_DIR) -Wl,-rpath,@loader_path/../lib

LDLIBS = -lgsl -lgslcblas -lgfortran -lminpack

# Source files
C_SOURCES = ../src/native/RCM_solv.c
FORTRAN_SOURCES = ../src/native/nifco.f90

# Object files
C_OBJS = $(C_SOURCES:.c=.o)
FORTRAN_OBJS = $(FORTRAN_SOURCES:.f90=.o)

# Shared library name
TARGET = ../lib/RCM_solver.so

# Default target
all: $(TARGET)

# Build shared library
$(TARGET): $(C_OBJS) $(FORTRAN_OBJS)
	$(CC) -shared -o $@ $^ $(LDFLAGS) $(LDLIBS)

# Compile C source files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Compile Fortran source files
%.o: %.f90
	$(FC) -fPIC -c $< -o $@ $(FFLAGS) 

# Clean up build artifacts
clean:
	rm -f $(C_OBJS) $(FORTRAN_OBJS) $(TARGET)

# Phony targets
.PHONY: all clean
