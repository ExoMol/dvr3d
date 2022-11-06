#!/bin/bash 

make clean 

make dvr.out 
./dvr.out < jobs/dvr_J0.job > result.H3p_J0.dvr 
cp -f fort.69 dip_bra_J0 
./dvr.out < jobs/dvr_J0_big.job > result.H3p_J0_big.dvr
./dvr.out < jobs/dvr_J1_0.job > result.H3p_J1_0.dvr
cp fort.26 fort.26.tmp
./dvr.out < jobs/dvr_J1_1.job > result.H3p_J1_1.dvr
cp fort.26 fort.27  
mv fort.26.tmp fort.26

make rot.out 
./rot.out < jobs/rot.job > result.H3p_J1.rot
cp -f fort.8 dip_ket_J1


cp dip_ket_J1 fort.11
cp  dip_bra_J0 fort.12

make dip.out 
./dip.out < jobs/dip.job > result.H3p.dip
#mv fort.96 result96.H3p.dip

