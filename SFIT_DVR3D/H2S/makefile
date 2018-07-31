goal:   dvr3d.x dvr3d_sfit.x rotlev.x xpect3.x

tarball:
	tar cf dvr_rot.tar makefile dvr3_rot.f90 rotlev3v.f90  f02fjf.f  h2o_ai.f dvr_pol.f90 wfread.f90



PLAT = _o3_0706
F90    = ifort 

FFLAGS = -O3  -ip  -traceback -static -openmp
LIB    =  -mkl=parallel

#LDIR = /opt/intel/composerxe-2011.4.191/mkl/lib/intel64/
#LAPACK =  -Wl,--start-group $(LDIR)libmkl_intel_lp64.a   \
# $(LDIR)libmkl_intel_thread.a $(LDIR)libmkl_core.a   -Wl,--end-group -lguide -lpthread

###############################################################################


OBJ    =   f02fjf.o  prop.o

xpect3.x:	xpect3.o pot_h2s.o prop.o pes_h2s_pcl.o
	$(F90) -o xpect3.x xpect3.o pot_h2s.o prop.o pes_h2s_pcl.o  $(FFLAGS) $(LIB)

dvr3d_sfit.x:	dvr3d_sfit.o pot_h2s.o pes_h2s_pcl.o
	$(F90) -o j-dvr3d_sfit.x dvr3d_sfit.o pot_h2s.o pes_h2s_pcl.o  $(FFLAGS) $(LIB)

rotlev.x:	f02fjf.o rotlev3b.o pes_h2s_pcl.o
	$(F90) -o rotlev.x rotlev3b.o f02fjf.o pes_h2s_pcl.o $(FFLAGS) $(LIB)

dvr3d.x:	 dvr3d.o  pot_h2s.o pes_h2s_pcl.o 
	$(F90) -o dvr3d.x dvr3d.o  pot_h2s.o pes_h2s_pcl.o $(FFLAGS) $(LIB)

pot_h2s.o:	pot_h2s.f
	$(F90) -c pot_h2s.f  $(FFLAGS)

pes_h2s_pcl.o:	pes_h2s_pcl.f90  pot_h2s.f
	$(F90) -c pes_h2s_pcl.f90  $(FFLAGS)

dvr3d_sfit.o:	dvr3d_sfit.f90
	$(F90) -c dvr3d_sfit.f90 $(FFLAGS)

prop.o:	prop.f90
	$(F90) -c prop.f90 $(FFLAGS)

xpect3.o:	xpect3.f90
	$(F90) -c xpect3.f90 $(FFLAGS)

f02fjf.o:	f02fjf.f
	$(F90) -c f02fjf.f $(FFLAGS)

rotlev3b.o:	rotlev3b.f90
	$(F90) -c rotlev3b.f90 $(FFLAGS)

dvr3d.o:	$(OBJ) dvr3drjz.f90
	$(F90) -c dvr3drjz.f90 -o dvr3d.o $(FFLAGS)

clean:
	rm *.o *.mod

