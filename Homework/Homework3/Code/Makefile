FC = gfortran
LD = gfortran
LDFLAGS = 
FFLAGS = -O3
F90FLAGS = -O3

OBJECTS = lglnodes.o gauss_approx.o gauss_approx_p2.o trapappx.o trapappx_p2.o

.PHONY: clean

all: gquad.x gquad2.x trappx.x trappx2.x

gquad.x: lglnodes.o gauss_approx.o
	$(LD) lglnodes.o gauss_approx.o -o gquad.x

gquad2.x: lglnodes.o gauss_approx_p2.o
	$(LD) lglnodes.o gauss_approx_p2.o -o gquad2.x

trappx.x: trapappx.o
	$(LD) trapappx.o -o trappx.x

trappx2.x: trapappx_p2.o
	$(LD) trapappx_p2.o -o trappx2.x

%.o : %.f90
	$(FC) $(F90FLAGS) -c $<

%.o : %.f
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) *.x *.eps *.txt *~
