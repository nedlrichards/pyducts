# Disable all of make's built-in rules (similar to Fortran's implicit none)
MAKEFLAGS += --no-builtin-rules --no-builtin-variables
# configuration
FC := gfortran
LD := $(FC)
RM := rm -f
# list of all source files

PADE := cmplx_roots_sg.o constants.o pade_coeffs.o
PADETEST := pade_test.o
PROG := ../bin/pade_test

RAMOBJ := ../ram_original/ram15.o
RAMPROG := ../bin/ram

RAMTOBJ := ram_typed.o
RAMTPROG := ../bin/ram_typed

all: $(PROG) $(RAMPROG) $(RAMTPROG)

$(PROG): $(PADE) $(PADETEST)
	gfortran -o $@ $^
$(RAMPROG): $(RAMOBJ)
	gfortran -o $@ $^
$(RAMTPROG): $(PADE) $(RAMTOBJ)
	gfortran -o $@ $^

$(PADE): %.o: %.f90
	gfortran -c -o $@ $<
$(PADETEST): %.o: %.f90
	gfortran -c -o $@ $<
$(RAMOBJ): %.o: %.f
	gfortran -c -o $@ $<
$(RAMTOBJ): %.o: %.f
	gfortran -c -o $@ $<
