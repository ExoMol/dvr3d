make clean 


make dvr.out 
./dvr.out < jobs/dvr.job > result.HCN.dvr
make rot.out 
./rot.out < jobs/rot.job > result.HCN.rot 
cp fort.8 fort.11 
cp fort.9 fort.12 
make dip.out 
./dip.out < jobs/dip.job > result.HCN.dip 
#make spc.out 
#./spc.out < jobs/spec.job > result.HCN.spec

