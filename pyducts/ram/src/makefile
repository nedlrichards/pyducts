# Disable all of make's built-in rules (similar to Fortran's implicit none)
MAKEFLAGS += --no-builtin-rules --no-builtin-variables
# configuration
FC := gfortran
LD := $(FC)
RM := rm -f
prefix := ../bin/

# list of all source files
#SRCS := cmplx_roots_sg.f90 constants.f90 pade_coeffs.f90 ram_v3.f90 ram_v4.f90 ram_io.f90
SRCS := cmplx_roots_sg.f90 constants.f90 pade_coeffs.f90 ram_v3.f90 ram_io.f90

RAM := ../original_ram/ram15.f
PROG_RAM := ram

RAM_V1 := ram_v1.f
PROG_RAM_V1 := ram_v1

RAM_V2 := ram_v2.f90
PROG_RAM_V2 := ram_v2

RAM_V3 := run_ram_v1.f90
PROG_RAM_V3 := ram_v3

RAM_V4 := run_ram_v2.f90
PROG_RAM_V4 := ram_v4

OBJS := $(addsuffix .o, $(SRCS))

.PHONY: all clean

all: $(PROG_RAM) $(PROG_RAM_V1) $(PROG_RAM_V2) $(PROG_RAM_V4)
#all: $(PROG_RAM) $(PROG_RAM_V1) $(PROG_RAM_V2) $(PROG_RAM_V3) $(PROG_RAM_V4)

$(PROG_RAM): $(RAM)
	$(LD) -o $(prefix)$@ $^

$(PROG_RAM_V1): $(RAM_V1) $(OBJS)
	$(LD) -o $(prefix)$@ $^

$(PROG_RAM_V2): $(RAM_V2) $(OBJS)
	$(LD) -o $(prefix)$@ $^

#$(PROG_RAM_V3): $(RAM_V3) $(OBJS)
	#$(LD) -o $(prefix)$@ $^

$(PROG_RAM_V4): $(RAM_V4) $(OBJS)
	$(LD) -o $(prefix)$@ $^

$(OBJS): %.o: %
	$(FC) -c -o $@ $<

run_ram_v2.f90.o: ram_v3.f90.o ram_io.f90 pade_coeffs.f90.o constants.f90.o
ram_io.f90.o: constants.f90.o
run_ram_v1.f90.o: ram_v3.f90.o pade_coeffs.f90.o constants.f90.o
ram_v4.f90.o: ram_io.f90.o pade_coeffs.f90.o constants.f90.o
ram_v3.f90.o: pade_coeffs.f90.o constants.f90.o
ram_v2.f90.o: pade_coeffs.f90.o constants.f90.o
ram_v1.f.o: pade_coeffs.f90.o
pade_coeffs.f90.o: cmplx_roots_sg.f90.o constants.f90.o

$(OBJS): $(MAKEFILE_LIST)

clean:
	$(RM) $(filter %.o, $(OBJS)) $(wildcard *.mod) $(prefix)$(PROG_RAM) $(prefix)$(PROG_RAM_V1) $(prefix)$(PROG_RAM_V2) $(prefix)$(PROG_RAM_V3)
