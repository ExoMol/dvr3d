F90      =ifort
DEBUG    = #-openmp -check all -traceback
F90_OPTS = -O0 -xP -traceback -check all -openmp #-xP
LIBS= -llapack #-L/homes/bruno/intel/mkl/9.0/lib/em64t -lmkl -lmkl_lapack64 #-lpthread 
REL_PATH =source/
OBJ =obj/
EXEC = res3d #_new_cap_build

$(EXEC): \
$(OBJ)class_wavefile_to_kblock.o                \
$(OBJ)class_kblock_to_elements.o                \
$(OBJ)class_input.o                             \
$(OBJ)class_complex_h.o                        \
$(OBJ)class_curvs.o                             \
$(OBJ)class_sort.o                              \
$(OBJ)class_stash.o                             \
$(OBJ)class_res_consistency.o                   \
$(OBJ)class_res_search.o                        \
$(OBJ)main_resonance.o                         
	$(F90) $(DEBUG) $(F90_OPTS) -o $@ $? $(LIBS)

$(OBJ)class_input.o                             : $(REL_PATH)class_input.f90
	cd $(OBJ) ;$(F90) $(DEBUG) $(F90_OPTS) -c ../$? 	     

$(OBJ)class_wavefile_to_kblock.o : $(REL_PATH)cap/xpect/class_wavefile_to_kblock.f90 
	cd $(OBJ) ;$(F90) $(DEBUG) $(F90_OPTS) -c ../$? 	     

$(OBJ)class_kblock_to_elements.o : $(REL_PATH)cap/xpect/class_kblock_to_elements.f90 
	cd $(OBJ) ;$(F90) $(DEBUG) $(F90_OPTS) -c ../$? 	     

$(OBJ)class_complex_h.o                        : $(REL_PATH)cap/class_complex_h.f90 
	cd $(OBJ) ;$(F90) $(DEBUG) $(F90_OPTS) -c ../$? 	     

$(OBJ)class_curvs.o                             : $(REL_PATH)curvs/class_curvs.f90 
	cd $(OBJ) ;$(F90) $(DEBUG) $(F90_OPTS) -c ../$? 	     

$(OBJ)class_sort.o                              : $(REL_PATH)utils/class_sort.f90 
	cd $(OBJ) ;$(F90) $(DEBUG) $(F90_OPTS) -c ../$? 	     

$(OBJ)class_stash.o                             : $(REL_PATH)utils/class_stash.f90
	cd $(OBJ) ;$(F90) $(DEBUG) $(F90_OPTS) -c ../$? 	     

$(OBJ)class_res_consistency.o                   : $(REL_PATH)search/class_res_consistency.f90 
	cd $(OBJ) ;$(F90) $(DEBUG) $(F90_OPTS) -c ../$? 	     

$(OBJ)class_res_search.o                        : $(REL_PATH)search/class_res_search.f90
	cd $(OBJ) ;$(F90) $(DEBUG) $(F90_OPTS) -c ../$? 	     

$(OBJ)main_resonance.o                          : $(REL_PATH)main_resonance.f90 
	cd $(OBJ) ;$(F90) $(DEBUG) $(F90_OPTS) -c ../$? 

clean : 
	rm $(OBJ)*.o $(OBJ)*.mod $(EXEC)
purge :
	rm *~ */*~ */*/*~ */*/*/*~ 
