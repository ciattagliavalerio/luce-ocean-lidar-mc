# ==========================================================================
#  Makefile — Monte Carlo code
#
#  make          builds the executable ./mc_water (optimised)
#  make debug    rebuilds with runtime checks (out-of-range indices, NaN, ...)
#  make clean    removes the executable and intermediate files (build/)
# ==========================================================================

FC     = gfortran
FFLAGS = -O2
DEBUG  = -O0 -g -fcheck=all -fbacktrace -ffpe-trap=invalid,zero,overflow

# Code version (git describe) written in the header of the .DAT files;
# MC.f90 goes through the preprocessor (-cpp) to receive the GIT_VERSION macro.
GITVER = $(shell git describe --tags --dirty --always 2>/dev/null || echo unknown)
CPPDEF = -cpp -DGIT_VERSION="'$(GITVER)'"

SRC    = src
BUILD  = build
EXE    = mc_water

all: $(EXE)

# The module is compiled first: MC.f90 depends on mc_module.mod
$(BUILD)/mc_module.o: $(SRC)/mc_module.f90 | $(BUILD)
	$(FC) $(FFLAGS) -J$(BUILD) -c $< -o $@

# MC.o is also rebuilt when the git state changes (version in the output)
$(BUILD)/MC.o: $(SRC)/MC.f90 $(BUILD)/mc_module.o $(wildcard .git/index)
	$(FC) $(FFLAGS) $(CPPDEF) -J$(BUILD) -c $< -o $@

$(EXE): $(BUILD)/mc_module.o $(BUILD)/MC.o
	$(FC) $(FFLAGS) $^ -o $@

$(BUILD):
	mkdir -p $@

debug: FFLAGS = $(DEBUG)
debug: clean all

clean:
	rm -rf $(BUILD) $(EXE)

.PHONY: all debug clean
