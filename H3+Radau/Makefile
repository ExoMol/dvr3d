# Run test (test) or cleanup (clean). Paolo 
# Test : J0 J1 dvr 
include ../var.dvr 
depo=scratch
scratch_files= dip_bra_J0 26_file 27_file dip_ket_J1

null:
	@echo Please specify test or clean

test :
	make J0 J1 dip J0big ; 

J0 : dvr.out 
	./dvr.out < jobs/dvr_J0.job > result.H3p_J0.dvr ; \
	cp -f fort.69 dip_bra_J0 ; 

J0big : dvr.out 
	./dvr.out < jobs/dvr_J0_big.job > result.H3p_J0_big.dvr ; \

J1 : 
	make J1_dvr J1_rot ; 

J1_dvr :  dvr.out 
	./dvr.out < jobs/dvr_J1_0.job > result.H3p_J1_0.dvr ; \
	cp fort.26 26_file  ;  \
	./dvr.out < jobs/dvr_J1_1.job > result.H3p_J1_1.dvr ; \
	cp fort.26 27_file  ;  

J1_rot : rot.out 26_file  27_file 
	@cp -f 26_file  fort.26 ; \
	cp -f 27_file  fort.27 ; \
	./rot.out < jobs/rot.job > result.H3p_J1.rot ; \
	cp -f fort.8 dip_ket_J1 ; 

dip :   dip_ket_J1  dip_bra_J0 dip.out 
	@cp dip_ket_J1 fort.11 ; \
	cp  dip_bra_J0 fort.12 ; \
	./dip.out < jobs/dip.job > result.H3p.dip ; \
	mv fort.96 result96.H3p.dip ; 

${scratch_files}  : 
	cp ${depo}/$@   ./  ;  

dvr.out:  dvr3drjz.o
	${F90} ${F90_OPTS} ../${OBJ}/dvr3drjz.o source/pot_h3p.f90  -o dvr.out ${LIBS} ; 

dip.out:  dipolez.o
	${F90} ${F90_OPTS} ../${OBJ}/dipolez.o source/dip.f90  -o dip.out ${LIBS} ; 

rot.out: rotlev3z.o  f02fjf.o 
	${F90} ${F90_OPTS}  ../${OBJ}/rotlev3z.o ../${OBJ}f02fjf.o  -o rot.out ${LIBS} ; \

%.o: 
	cd ../ ;  \
	make $@ ; 

clean:
	@-rm fort.* ${scratch_files} result.H3p* result96.H3p.dip   *out *mod 


