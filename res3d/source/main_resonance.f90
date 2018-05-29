program main
use class_input
use class_complex_h
use class_curvs
use class_res_search
use class_res_consistency
implicit none
type(input) :: in
write(6,*)
write(6,*)
write(6,"(10X,'Program Resonance, version December, 2008')")
write(6,"(10X,'Written by B. C. Silva, P. Barletta and J. J. Munro')")
write(6,*)
write(6,*)
call read_input(in)

call cap_matrix_build(in)
!
call curvs_build(in)

call res_search(in)
write(6,*)
write(6,*)
write(6,"(10X,'Done!!')")
end program main


    