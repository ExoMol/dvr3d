F90      =ifort
DEBUG    = #-openmp -check all -traceback
F90_OPTS = -O0  -traceback -check all -qopenmp
LIBS= -qmkl=parallel  # -llapack 
REL_PATH =source/
OBJ =source

test:
	make dvr.out rot.out dip.out
	
null:
	@echo Please specify test or clean


dvr.out:  ../source/dvr3drjz.f90
	${F90} ${F90_OPTS} ../${OBJ}/dvr3drjz.f90 source/pot_h3p.f90  -o dvr.out ${LIBS} ; 

dip.out:  ../source/dipolez.f90
	${F90} ${F90_OPTS} ../${OBJ}/dipolez.f90 source/dip.f90  -o dip.out ${LIBS} ; 

rot.out: ../source/rotlev3z.f90  ../source/f02fjf.f 
	${F90} ${F90_OPTS}  ../${OBJ}/rotlev3z.f90 ../${OBJ}/f02fjf.f  -o rot.out ${LIBS} ; \


clean:
	@-rm fort.*  result.H3p* result96.H3p.dip   *out *mod 


