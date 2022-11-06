include var.dvr 
EXEC = dvr3d rot3b rot3 rot3z xpct dipx dipz spec
OBJS = dipj0dvr.o dipole.o dipolez.o rotlev3b.o spectra.o xpect3.o dvr3drjz.o rotlev3.o rotlev3b.o dipd_JT.o f02fjf.o 
TEST_DIR = water  HCN  H3+Radau   H3+Jacobi  

###$(EXEC): \
###$(OBJ)class_wavefile_to_kblock.o                \
###$(OBJ)class_kblock_to_elements.o                \
###$(OBJ)class_input.o                             \
###$(OBJ)class_complex_h.o                        \
###$(OBJ)class_curvs.o                             \
###$(OBJ)class_sort.o                              \
###$(OBJ)class_stash.o                             \
###$(OBJ)class_res_consistency.o                   \
###$(OBJ)class_res_search.o                        \
###$(OBJ)main_resonance.o                         
###	$(F90) $(DEBUG) $(F90_OPTS) -o $@ $? $(LIBS)

test : ${TEST_DIR} 

water : 
	cd water; \
	make rot.out ; 

HCN : 
	cd HCN; \
	make test ; 

H3+Radau : 
	cd H3+Radau ; \
	make test ; 

H3+Jacobi : 
	cd H3+Jacobi ; \
	make test ; 

.PHONY: ${TEST_DIR}

objects : $(OBJS)

%.o : ${REL_PATH}%.f90 
	cd $(OBJ) ; \
	$(F90) $(DEBUG) $(F90_OPTS) -c ../$< -o $@

%.o : ${REL_PATH}%.f 
	cd $(OBJ) ; \
	$(F90) $(DEBUG) $(F90_OPTS) -c ../$< -o $@

dvr3d: ${OBJ}dvr3drjz.o pot.f90 
	$(F90) $(DEBUG) $(F90_OPTS) $(LIBS) $(OBJ)dvr3drjz.o pot.f90 -o dvr3d 	     

rot3: rotlev3.o f02fjf.o
	$(F90) $(DEBUG) $(F90_OPTS) $(LIBS) $(OBJ)rotlev3.o $(OBJ)f02fjf.o -o rot3 	     

rot3b: rotlev3b.o f02fjf.o
	$(F90) $(DEBUG) $(F90_OPTS) $(LIBS) $(OBJ)rotlev3b.o $(OBJ)f02fjf.o  -o rot3b  	     

rot3z: rotlev3z.o f02fjf.o
	$(F90) $(DEBUG) $(F90_OPTS) $(LIBS) $(OBJ)rotlev3z.o $(OBJ)f02fjf.o  -o rot3z  	     

dipx: dipole.o dipd.f90 
	$(F90) $(DEBUG) $(F90_OPTS) $(LIBS) $(OBJ)dipole.o  -o dipx  	     

dipz: dipolez.o dipd.f90 
	$(F90) $(DEBUG) $(F90_OPTS) $(LIBS) $(OBJ)dipole.o  -o dipz  	     

spec: spectra.o dpsort.o
	$(F90) $(DEBUG) $(F90_OPTS) $(LIBS) $(OBJ)spectra.o $(OBJ)dpsort.o -o spec  	     

xpct: xpect3.o prop.f90 
	$(F90) $(DEBUG) $(F90_OPTS) $(LIBS) $(OBJ)xpect3.o prop.f90 -o xpct  	     

dipd.f90:
	@echo "The user must supply a dipd.f90 routine in the working dir" ;\
	-exit 1 ; 

prop.f90:
	@echo "The user must supply a prop.f90 routine in the working dir"

pot.f90:
	@echo "The user must supply a pot.f90 routine in the working dir"

.PHONY: pot.f90 prop.f90 dipd.f dipd.f90 

clean : 
	@for DIR in $(TEST_DIR); do  \
	cd $$DIR; make clean; cd ..; done 
	-rm $(OBJ)*.o $(OBJ)*.mod $(EXEC) ; \

purge :
	rm *~ */*~ */*/*~ */*/*/*~ 



