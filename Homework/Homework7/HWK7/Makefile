# Makefile
# 
# 
# 
# 
# 

FC = mpif90
LD = mpif90
LDFLAGS =
FFLAGS = 
F90FLAGS = -Wall -fbounds-check -g -fbacktrace
EX = ./homework7.x

OBJECTS = xycoord.o homework7.o

.PHONY: clean

$(EX): $(OBJECTS) 
	$(LD) $(OBJECTS) -o $(EX) $(LDFLAGS) 

%.o : %.f90
	$(FC) $(F90FLAGS) -c $<

%.o : %.f
	$(FC) $(FFLAGS) -c $<

%.mod : %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f *.o *.txt *.x
