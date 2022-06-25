
AUTOMAKE_OPTIONS = foreign

###############################################################################
# THE PROGRAMS TO BUILD
###############################################################################

# the program to build (the names of the final binaries)
bin_PROGRAMS = MAMMA
bindir=@bindir@

# list of sources for the binary
MAMMA_SOURCES = \
                geometry.f90 \
                parameters.f90 \
                melts_fit.f90 \
		complexify.f90 \
		constitutive.f90 \
		equations.f90 \
                init.f90 \
                inpout.f90 \
                steady_solver.f90 \
                c_main.f90 

FCFLAGS         = -Wunused -O3 -funroll-all-loops
FCDEBUGFLAGS    = -Wall -g

LIBS            = -lm -llapack -lblas

FC              = gfortran
MPFC            = gfortran
LINKER          = gfortran
