
# Fortran 90 complier to use:
COMPILER = gfortran

# Target Architecture (Otherwise use -xHost)
# 	Skylake:	-xCORE-AVX512
# 	KNL: 			-xMIC-AVX512
TARGET_ARCH = #-xMIC-AVX512

# Libraries to Link
LDFLAGS = -L/share/software/user/open/fftw/3.3.10/lib

# Include Directories
INCLUDES = -I/share/software/user/open/fftw/3.3.10/include

# Compiler Options
FCFLAGS = -O3 -DSHARED_MEMORY
#FCFLAGS = -O0 -g -check all -traceback -debug extended # Debug


# --- # --- # --- # --- # --- # --- # --- # --- #

LDFLAGS += -lfftw3 -lgfortran

# FCFLAGS += $(TARGET_ARCH) -DHDF5
FCFLAGS += $(INCLUDES)
ifeq ($(COMPILER),ifort)
  FCFLAGS += -fpp -ipo -xHost
else ifeq ($(COMPILER),gfortran)
  FCFLAGS += -cpp -march=native -flto
endif

COMPILER = h5pfc

SOURCES = domain.o fft.o flow.o parameters.o tools.o

diablo: diablo.f90 $(SOURCES)
	$(COMPILER) $(FCFLAGS) $^ -o $@ $(LDFLAGS)

%.o: %.f90
	$(COMPILER) $(FCFLAGS) -c $<



domain.o: domain.f90 grid_def grid_mpi

fft.o: fft.f90 domain.o

flow.o: flow.f90 fft.o domain.o parameters.o tools.o \
				diagnostics.f90 les.f90 phdf5.f90 solver.f90

parameters.o: parameters.f90 fft.o domain.o

tools.o: tools.f90 domain.o parameters.o


.PHONY: clean
clean:
	-rm -f *.o *.mod diablo
