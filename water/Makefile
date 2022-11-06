# Run test (test) or cleanup (clean). Paolo 
# Test : J0 J2 
#  J0 :   dvr  xpt 
#  J2 :   dvr  rot  dip  xpt
include ../var.dvr 

null:
	@echo Please specify test or clean

test :
	make dvr.out rot.out xpect.out dip.out  \


dvr.out:  dvr3drjz.o
	${F90} ${F90_OPTS} ../${OBJ}/dvr3drjz.o source/pot_water.f  -o dvr.out ${LIBS} ; 

rot.out: rotlev3b.o  f02fjf.o 
	${F90} ${F90_OPTS}  ../${OBJ}/rotlev3b.o ../${OBJ}f02fjf.o  -o rot.out ${LIBS} ; \

xpect.out : xpect3.o
	${F90}  ${F90_OPTS} ../source/xpect3.f90 source/prop_water.f90  -o xpect.out ${LIBS} ;  \

dip.out:  dipole3.o
	${F90}  ${F90_OPTS} ../${OBJ}/dipole3.o  source/dip_water.f  -o dip.out ${LIBS} ; \

%.o: 
	cd ../ ;  \
	make $@ ; 

clean:
	@-rm fort* ${scratch_files}  result* *out *mod 


