#!/bin/bash 

make clean 

make dvr.out
./dvr.out < jobs/dvr0.job > result-0.water.dvr 
cp fort.26 fort.11

make xpect.out 
./xpect.out < jobs/xpect0.job > result-0.water.xpect 

./dvr.out < jobs/dvr-even.job > result-even.water.dvr 
cp fort.26 26_even 
./dvr.out < jobs/dvr-odd.job > result-odd.water.dvr
cp fort.26 26_odd  

make rot.out 
cp -f 26_even fort.26
cp -f 26_even fort.4
./rot.out < jobs/rot.job > result-even.water.rot  
cp fort.8 fort-e.8 
cp fort.9 fort-e.9  
cp -f 26_odd fort.26
cp -f 26_odd fort.4 
./rot.out < jobs/rot.job > result-odd.water.rot 
cp fort.8 fort-o.8 
cp fort.9 fort-o.9  
./xpect.out < jobs/xpect2.job > result-2.water.xpect 

make dip.out 
cp fort-e.8 fort.11
cp fort-o.9 fort.12
./dip.out < jobs/dip.job > result.water.dip 
