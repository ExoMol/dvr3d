# Run test (test) or cleanup (clean). Paolo 
# Test : J0 J2 
#  J0 :   dvr  xpt 
#  J2 :   dvr  rot  dip  xpt
include ../var.dvr 
depo=scratch/
scratch_files= 26_file 26_even 26_odd  fort-e.8 fort-o.8 fort-e.9  fort-o.9

null:
	@echo Please specify test or clean

test :
	make J0 J2 \

J0 :   
	make J0_dvr ; \
	make J0_xpt ; 

J2 :  
	make J2_dvr ; \
	make J2_rot ; \
	make J2_dip ; \
	make J2_xpt ; 

J0_dvr :  dvr.out 
	./dvr.out < jobs/dvr0.job > result-0.water.dvr ; \
	cp fort.26 26_file ; 

J0_xpt :  xpect.out 26_file 
	cp 26_file fort.11 ; \
	./xpect.out < jobs/xpect0.job > result-0.water.xpect ; \

J2_dvr :  dvr.out 
	./dvr.out < jobs/dvr-even.job > result-even.water.dvr ; 
	cp fort.26 26_even ; \
	./dvr.out < jobs/dvr-odd.job > result-odd.water.dvr ; \
	cp fort.26 26_odd ; 

J2_rot : rot.out 26_even 26_odd 
	@cp -f 26_even fort.26 ; \
	cp -f 26_even fort.4 ; \
	./rot.out < jobs/rot.job > result-even.water.rot ; 
	cp fort.8 fort-e.8 ; \
	cp fort.9 fort-e.9 ; 
	cp -f 26_odd fort.26 ; \
	cp -f 26_odd fort.4 ; \
	./rot.out < jobs/rot.job > result-odd.water.rot ; \
	cp fort.8 fort-o.8 ; \
	cp fort.9 fort-o.9 ; 

J2_xpt : xpect.out 
	./xpect.out < jobs/xpect2.job > result-2.water.xpect ; \

J2_dip :  fort-e.8 fort-o.9 dip.out 
	cp fort-e.8 fort.11 ; \
	cp fort-o.9 fort.12 ; \
	./dip.out < jobs/dip.job > result.water.dip ; \

${scratch_files}  :
	cp ${depo}/$@   ./  ;

dvr.out:  dvr3drjz.o
	${F90} ${F90_OPTS} ../${OBJ}/dvr3drjz.o source/pot_water.f  -o dvr.out ${LIBS} ; 

rot.out: rotlev3b.o  f02fjf.o 
	${F90} ${F90_OPTS}  ../${OBJ}/rotlev3b.o ../${OBJ}f02fjf.o  -o rot.out ${LIBS} ; \

xpect.out : xpect3.o
	${F90}  ${F90_OPTS} ../${OBJ}/xpect3.o source/prop_water.f90  -o xpect.out ${LIBS} ;  \

dip.out:  dipole.o
	${F90}  ${F90_OPTS} ../${OBJ}/dipole.o  source/dip_water.f  -o dip.out ${LIBS} ; \

%.o: 
	cd ../ ;  \
	make $@ ; 

clean:
	@-rm fort* ${scratch_files}  result* *out *mod 


