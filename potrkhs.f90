module RKHSPES_MOD
contains
SUBROUTINE RKHSPES(X,ENER)
use RKHS      
implicit none

!SOME VARIABLES ARE DEFINED HERE
REAL(KIND(0D0))    :: X(:)        ! COORDINATE ARRAY (EVEN WHEN D = 1, THIS NEEDS TO BE AN ARRAY)
REAL(KIND(0D0))    :: ENER        ! WILL STORE THE VALUE OF THE MODEL FUNCTION 
TYPE(KERNEL),SAVE  :: K           ! THE KERNEL TYPE IS NEEDED TO SET UP AND EVALUATE A RKHS MODEL
LOGICAL, SAVE :: IS_INITIALIZED = .FALSE.

IF (.NOT.IS_INITIALIZED) THEN
!CALL K%READ_GRID("ABINITIO")
!======================================================================================================================
!CALL K%K1D(1)%INIT(TAYLOR_SPLINE_N2_KERNEL)   ! CHOOSE ONE-DIMENSIONAL KERNEL FOR DIMENSION 1
!call K%k1d(2)%init(RECIPROCAL_POWER_N2_M4_KERNEL)       ! choose one-dimensional kernel for dimension 2
!======================================================================================================================
!CALL K%CALCULATE_COEFFICIENTS_FAST() ! ONCE ALL KERNELS ARE INITIALIZED, THE COEFFICIENTS CAN BE CALCULATED (USES TENSOR PRODUCT)
!CALL K%CALCULATE_SUMS() ! THIS CALL WILL CALCULATE THE COMPLETE LOOKUP TABLE
!call K%save_to_file("POT.KERNEL")
CALL K%LOAD_FROM_FILE("POT.KERNEL")
IS_INITIALIZED = .TRUE.
END IF
CALL K%EVALUATE_FAST(X,ENER)


END SUBROUTINE
end module RKHSPES_MOD

SUBROUTINE POTV( V, SMLR, CAPR, CT )
USE RKHSPES_MOD
IMPLICIT NONE
REAL*8, DIMENSION ( : ) :: R ( 2 )
REAL*8 :: V
REAL*8 :: CAPR, SMLR, CT
REAL*8, PARAMETER ::  PI=ACOS(-1.0D0)
!SOME VARIABLES ARE DEFINED HERE
!R(1)=(1.0D0-COS(PI-ACOS(CT)))/2.0D0
R(1)=(1.0D0-CT)/2.0D0
R(2)=CAPR
CALL RKHSPES(R,V)
IF (V>0.0D0)V=0.0D0
END SUBROUTINE


