AUTOMAKE_OPTIONS = foreign

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
bin_PROGRAMS = PLUMEMoM
bindir=@bindir@

# list of sources for the binary
PLUMEMoM_SOURCES = \
		moments.f90 \
		plume.f90 \
                variables.f90 \
		meteo.f90 \
		particles.f90 \
                mixture.f90 \
                solver_rise.f90 \
                inpout.f90 \
		rise.f90 \
                plume_model.f90

FCFLAGS         = -g -O3 -ffast-math
FCDEBUGFLAGS    = -Wall -g

LIBS            = -lz -lm -llapack

FC              = gfortran
MPFC            = gfortran
LINKER          = gfortran

