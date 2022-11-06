#!/bin/bash 

for i in {1..100}
do

freq=$(sed "${i}q;d"  tmp2| awk '{print $4}' ) 
abscoef=$(sed "${i}q;d"  tmp2| awk '{print $3}' )
init=$(sed "${i}q;d"  tmp2| awk '{print $1}' )
final=$(sed "${i}q;d"  tmp2| awk '{print $2}' )
e_init=$(sed "${init}q;d"  SiO2_J6_all.states | awk '{print $2}' )
e_final=$(sed "${final}q;d"  SiO2_J6_all.states| awk '{print $2}' )

#line=$(grep "$freq" dipole.dat.com)

ofreq=$(grep "$freq" dipole.data.com | awk '{print $5}')
oabscoef=$(grep "$freq" dipole.data.com | awk '{print $10}')
oe_init=$(grep "$freq" dipole.data.com | awk '{print $3}')
oe_final=$(grep "$freq" dipole.data.com | awk '{print $4}')

echo "From DVR3D:    from Eini= $oe_init cm-1 to Efin= $oe_final cm-1. Transition freq: $ofreq cm-1 with abs coefficient of $oabscoef"
echo "From linelist: from Eini= $e_init cm-1 to Efin= $e_final cm-1. Transition freq: $freq cm-1 with abs coefficient of $abscoef"   
echo " "
#echo $abscoef $oabscoef
done 
