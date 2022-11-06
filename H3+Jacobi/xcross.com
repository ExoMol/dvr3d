Temperature  100 
Range 0.0  10000.0

Npoints 2001

absorption
(gaussian)
(HWHM 10. (cm-1))


output tmp 

States SiO2_J6_all.states 
Transitions  file_sorted.trans 

