# include Machine_PGF90
# include Machine_intel
include Machine_GFORTRAN
# include Machine_colima

PROG = MAMMA

SRCS = complexify constitutive geometry init inpout parameters steady_solver

OBJS = complexify.o constitutive.o geometry.o init.o inpout.o parameters.o steady_solver.o

all: MAMMA.x

MAMMA.x: $(OBJS) c_main.o 
	 $(LINKER) -o MAMMA.x $(MPFFLAGS) $(LINKFLAGS) $(OBJS) c_main.o $(LIBS)

clean:
	rm -f *.o *.mod *.a core core.* *.stb *~

include dependencies

