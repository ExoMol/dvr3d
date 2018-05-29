# Run test (test) or cleanup (clean). Paolo 
# Test : dvr rot dip spt 
include ../var.dvr 
depo=scratch

null:
	@echo Please specify test or clean

test :
	make dvr rot dip spt \

dvr :  dvr.out 
	./dvr.out < jobs/dvr.job > result.HCN.dvr ; \

rot : rot.out 26_file 
	@cp -f 26_file fort.26 ; \
	./rot.out < jobs/rot.job > result.HCN.rot ; \
	cp fort.8 8_file ; \
	cp fort.9 9_file ; 

dip :  8_file  9_file  dip.out 
	@cp 8_file fort.11 ; \
	cp 9_file fort.12 ; \
	./dip.out < jobs/dip.job > result.HCN.dip ; \

spt : spc.out 
	./spc.out < jobs/spec.job > result.HCN.spec ; \

%_file : 
	@cp ${depo}/$@  ./ 

dvr.out:  dvr3drjz.o
	${F90} ${F90_OPTS} ../${OBJ}/dvr3drjz.o source/pot_HCN.f  -o dvr.out ${LIBS} ; 

rot.out: rotlev3.o  f02fjf.o 
	${F90} ${F90_OPTS}  ../${OBJ}/rotlev3.o ../${OBJ}f02fjf.o  -o rot.out ${LIBS} ; \

spc.out : spectra.o  dpsort.o 
	${F90}  ${F90_OPTS} ../${OBJ}/spectra.o ../${OBJ}/dpsort.o  -o spc.out ${LIBS} ;  \

dip.out:  dipole.o
	${F90}  ${F90_OPTS} ../${OBJ}/dipole.o  source/dip_HCN.f  -o dip.out ${LIBS} ; \

%.o: 
	cd ../ ;  \
	make $@ ; 

clean:
	@-rm fort.* *_file result.HCN* *out *mod 


