# Run test (test) or cleanup (clean). Paolo 
# Test : dvr rot dip spt 
F90      =ifort
DEBUG    = #-openmp -check all -traceback
F90_OPTS = -O0  -traceback -check all
LIBS= -qmkl=sequential  # -llapack 
REL_PATH =source/
OBJ =obj/

null:
	@echo Please specify linelist or clean

linelist:  
	${F90} ${F90_OPTS} input.f90 linelist.f90  -o linelist.x 


clean:
	@-rm  *.x *mod 



