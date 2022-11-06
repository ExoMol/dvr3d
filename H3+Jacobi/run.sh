#!/bin/bash
make clean

make dvr.out 
./dvr.out < jobs/dvrJ0.job > result.H3p.J0.dvr  
cp fort.26 fort26_h3p_j0_odd  
./dvr.out < jobs/dvrJ1.job > result.H3p.J1even.dvr  
cp fort.26 fort26_h3p_j1_even 

make rot.out 
cp -f fort26_h3p_j1_even fort.26 
./rot.out < jobs/rot.job > result.H3p.J1even.rot 

cp fort.8 fort8_h3p_j1_even  
cp fort.9 fort9_h3p_j1_even  

cp fort26_h3p_j0_odd fort.11  
cp  fort8_h3p_j1_even fort.12  

make dip.out 
./dip.out < jobs/dip.job > result.H3p.dip 
make spc.out 
./spc.out < jobs/spec.job > result.H3p.spec 
