# Disable all of make's built-in rules (similar to Fortran's implicit none)
MAKEFLAGS += --no-builtin-rules --no-builtin-variables
# configuration
FC := gfortran
LD := $(FC)
RM := rm -f
prefix := ../bin/

# list of all source files
SRCS := cmplx_roots_sg.f90 constants.f90 pade_coeffs.f90
RAM_TYPED := ram_typed.f
PROG_TYPED := ram_typed

RAM := ../original_ram/ram15.f
PROG_RAM := ram

RAM_V1 := ram_v1.f
PROG_RAM_V1 := ram_v1

OBJS := $(addsuffix .o, $(SRCS))

.PHONY: all clean

all: $(PROG_RAM) $(PROG_TYPED) $(PROG_RAM_V1)

$(PROG_RAM): $(RAM)
	$(LD) -o $(prefix)$@ $^

$(PROG_RAM_V1): $(RAM_V1) $(OBJS)
	$(LD) -o $(prefix)$@ $^

$(PROG_TYPED): $(RAM_TYPED) $(OBJS)
	$(LD) -o $(prefix)$@ $^

$(OBJS): %.o: %
	$(FC) -c -o $@ $<

ram_typed.f.o: pade_coeffs.f90.o
pade_coeffs.f90.o: cmplx_roots_sg.f90.o constants.f90.o

$(OBJS): $(MAKEFILE_LIST)

clean:
	$(RM) $(filter %.o, $(OBJS)) $(wildcard *.mod) $(prefix)$(PROG_TYPED) $(prefix)$(PROG_RAM)
