# Makefile for Homework 6
# 
# 
# 
# 
# 

FC = gfortran
LD = gfortran
LDFLAGS = -llapack -lblas
FFLAGS = 
F90FLAGS = -O3 -g -fbounds-check -fbacktrace 
EX = ./homework6.x


OBJECTS = homework6.o rhside.o printdble.o 
OBJECTS2 = getmatrix.o rhside.o printdble.o 
MODULES = 

.PHONY: clean

# Compile, run, process and open.
$(EX): $(MODULES) $(OBJECTS) 
	$(LD) $(LDFLAGS) $(OBJECTS) -o $(EX)

getmat.x: $(MODULES) $(OBJECTS2) 
	$(LD) $(LDFLAGS) $(OBJECTS2) -o getmat.x

%.o : %.f90
	$(FC) $(F90FLAGS) -c $<

%.o : %.f
	$(FC) $(FFLAGS) -c $<

%.mod : %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) $(OBJECTS2) $(MODULES) $(EX) *.txt
