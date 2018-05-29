# Run test (test) or cleanup (clean). Paolo 
# Test : dvr rot dip spt 
include ../var.dvr 
depo=scratch
scratch_files= fort13_h3p_j0j1  fort26_h3p_j0_odd fort26_h3p_j1_even  fort8_h3p_j1_even fort9_h3p_j1_even  


null:
	@echo Please specify test or clean

test :
	make J0 J1 dip spt  \

J0 : dvr.out 
	./dvr.out < jobs/dvrJ0.job > result.H3p.J0.dvr  ; \
	cp fort.26 fort26_h3p_j0_odd ; 

J1 : 
	make J1_dvr J1_rot 

J1_dvr :  dvr.out 
	./dvr.out < jobs/dvrJ1.job > result.H3p.J1even.dvr  ; \
	cp fort.26 fort26_h3p_j1_even 

J1_rot : rot.out fort26_h3p_j1_even
	@cp -f fort26_h3p_j1_even fort.26 ; \
	./rot.out < jobs/rot.job > result.H3p.J1even.rot ; \
	cp fort.8 fort8_h3p_j1_even  ; \
	cp fort.9 fort9_h3p_j1_even  ; 

dip :  fort26_h3p_j0_odd fort8_h3p_j1_even  dip.out 
	@cp fort26_h3p_j0_odd fort.11  ;  \
	cp  fort8_h3p_j1_even fort.12  ;  \
	./dip.out < jobs/dip.job > result.H3p.dip ;  \

spt : spc.out fort13_h3p_j0j1
	@cp fort13_h3p_j0j1 fort.13 ; \
	./spc.out < jobs/spec.job > result.H3p.spec ; 

%_file : 
	@cp ${depo}/$@  ./ 

dvr.out:  dvr3drjz.o
	${F90} ${F90_OPTS} ../${OBJ}/dvr3drjz.o source/pot_h3p.f90  -o dvr.out ${LIBS} ; 

rot.out: rotlev3.o  f02fjf.o 
	${F90} ${F90_OPTS}  ../${OBJ}/rotlev3.o ../${OBJ}f02fjf.o  -o rot.out ${LIBS} ; \

spc.out : spectra.o  dpsort.o
	${F90}  ${F90_OPTS} ../${OBJ}/spectra.o ../${OBJ}/dpsort.o  -o spc.out ${LIBS} ;  \

dip.out:  dipole.o
	${F90}  ${F90_OPTS} ../${OBJ}/dipole.o  source/dip_h3p.f  -o dip.out ${LIBS} ; \

%.o: 
	cd ../ ;  \
	make $@ ; 

${scratch_files}  : 
	cp ${depo}/$@   ./  ;  

clean:
	@-rm fort* *_file result.H3p* *out *mod 


