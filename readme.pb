DVR3D: a program suite for the calculation of rotation-vibration
spectra of triatomic molecules

The files are stuctured into directory as follows: 

source : original fortran files
obj    : scratch directory containing the object files from compilation
water, 
HCN, 
H3+Jacobi, 
H3+Radau : test directories (details below).


All operations can be controlloed from the Makefile.
In particular : 

make x  (x is one of the test directories) will run the associated test.
make test will run all four tests serially

make code will compile the code, where code is one of : 

dvr3d     ...   dvr3drjz.f90    pot.f90
rot3      ...   rotlev3.f90 
rot3b     ...   rotlev3b.f90 
rot3z     ...   rotlev3z.f90 
dipx      ...   dipole.f90      dipd.f90
dipz      ...   dipolez.f90     dipd.f90
spec      ...   spectra.f90
xpct      ...   xpect3.f90      prop.f90

Please note that the user need supplying the files in the third column
in the compiling directories.


Details of test runs : 

All test runs can be controlled by the Makefile provided in each 
directory. The Makefile will also compile the relevant codes.
All test directories are organised into four subdir : 

jobs :     contains inputs for the jobs
results :  contains the test outcome files for comparaison
source  :  contains the specific property routines 
scratch : contains scratch files 

In all dir make test will run the test, and make clean will clean the dir.
More specifically, 

water : 

J0 : 
     J0_dvr   : vibrational level calculation
     J0_xpt   : expectation (more)

J2 : 
     J2_dvr   : first step in rotational calculation (J=2)
     J2_rot   : second step
     J2_dip   : J=2->J=2 transition
     J2_xpt     expectation (more)


HCN : 

dvr : First step of rot calc (J=2)
rot : second step
dip : J=2->J=2 transition
spt : spectrum at 2000 K


H3+Jacobi : 

J0  : vibrational level calculation

J1  : 
      J1_dvr : first step in rot calc J=1
      J1_rot : second step

dip : 
      J0->J1 dipole intensities calculation

spt : spectrum at 1400 K 	


H3+Radau : 

J0  : vibrational level calculation

J1  : 
      J1_dvr : first step in rot calc J=1
      J1_rot : second step

dip : 
      J0->J1 dipole intensities calculation
