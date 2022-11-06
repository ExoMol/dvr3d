F90      =ifort
DEBUG    = #-openmp -check all -traceback
F90_OPTS = -O0  -traceback -check all -qopenmp
LIBS= -qmkl=parallel  # -llapack 
REL_PATH =source/
OBJ =source

null:
	@echo Please specify test or clean

test :
	make dvr.out rot.out dip.out spc.out \


dvr.out:  ../source/dvr3drjz.f90
	${F90} ${F90_OPTS} ../${OBJ}/dvr3drjz.f90 source/pot_h3p.f90  -o dvr.out ${LIBS} ; 

rot.out: ../source/rotlev3.f90  ../source/f02fjf.f 
	${F90} ${F90_OPTS}  ../${OBJ}/rotlev3.f90 ../${OBJ}/f02fjf.f  -o rot.out ${LIBS} ; \

spc.out : ../source/spectra.f90  ../source/dpsort.f
	${F90}  ${F90_OPTS} ../${OBJ}/spectra.f90 ../${OBJ}/dpsort.f  -o spc.out ${LIBS} ;  \

dip.out:  ../source/dipole3.f90
	${F90}  ${F90_OPTS} ../${OBJ}/dipole3.f90  source/dip_h3p.f  -o dip.out ${LIBS} ; \


clean:
	@-rm fort*  result.H3p* *out *mod 


