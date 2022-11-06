# Run test (test) or cleanup (clean). Paolo 
# Test : dvr rot dip spt 
F90      =ifort
DEBUG    = #-openmp -check all -traceback
F90_OPTS = -O0  -traceback -check all 
LIBS= -qmkl=sequential  # -llapack 
REL_PATH =source/
OBJ =obj/

null:
	@echo Please specify test or clean

test :
	make dvr.out rot.out dip.out spc.out \

dvr.out:  dvr3drjz.o
	${F90} ${F90_OPTS} ../${OBJ}/dvr3drjz.o source/pot_HCN.f  -o dvr.out ${LIBS} ; 

rot.out: rotlev3.o  f02fjf.o 
	${F90} ${F90_OPTS}  ../${OBJ}/rotlev3.o ../${OBJ}f02fjf.o  -o rot.out ${LIBS} ; \

spc.out : spectra.o  dpsort.o 
	${F90}  ${F90_OPTS} ../${OBJ}/spectra.o ../${OBJ}/dpsort.o  -o spc.out ${LIBS} ;  \

dip.out:  dipole3.o
	${F90}  ${F90_OPTS} ../${OBJ}/dipole3.o  source/dip_HCN.f  -o dip.out ${LIBS} ; \

%.o: 
	cd ../ ;  \
	make $@ ; 

clean:
	@-rm fort.* *_file result.HCN* *out *mod 


