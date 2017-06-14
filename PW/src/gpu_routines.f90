
! Copyright (C) 2002-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------=!
#ifdef USE_CUDA
! gpu_routines module contains subroutines which wrap problematic CUF kernels in cdiaghg.F90.
! This is a temporary workaround.
module gpu_routines
!=----------------------------------------------------------------------=!
contains

  SUBROUTINE cufMemcpy2D( dst, dpitch, src, spitch, n, m )
    USE kinds
    USE cudafor
    IMPLICIT NONE
    COMPLEX(DP), DEVICE  :: dst(1:dpitch, *)
    COMPLEX(DP), DEVICE, INTENT(IN)  :: src(1:spitch, *)
    INTEGER, INTENT(IN) :: dpitch, spitch, n, m
    INTEGER :: i, j

    !$cuf kernel do(2) <<<*,*>>>
    DO j = 1, m
      DO i = 1, n
       dst(i,j) = src(i,j)
      END DO
    END DO
  END SUBROUTINE cufMemcpy2D

  SUBROUTINE copy_diag( dst, src, spitch, n)
    USE kinds
    USE cudafor
    IMPLICIT NONE
    REAL(DP), DEVICE  :: dst(*) 
    COMPLEX(DP), DEVICE, INTENT(IN)  :: src(1:spitch, *)
    INTEGER, INTENT(IN) :: spitch, n
    INTEGER :: i

    !$cuf kernel do(1) <<<*,*>>>
    DO i = 1, n
      dst(i) = DBLE( src(i,i) )
    END DO
  END SUBROUTINE copy_diag

  SUBROUTINE restore_upper_tri( a, lda, diag, n)
    USE kinds
    USE cudafor
    IMPLICIT NONE
    REAL(DP), DEVICE, INTENT(IN)  :: diag(*) 
    COMPLEX(DP), DEVICE  :: a(1:lda, *)
    INTEGER, INTENT(IN) :: lda, n
    INTEGER :: i, j

    !$cuf kernel do(1) <<<*,*>>>
    DO i = 1, n
      a(i,i) = DCMPLX( diag(i), 0.0_DP)
      DO j = i + 1, n
        a(i,j) = DCONJG( a(j,i) )
      END DO
      ! JR Zeroing out upper tri isn't necessary
      !DO j = n + 1, ldh
      !  s(j,i) = ( 0.0_DP, 0.0_DP )
      !END DO
    END DO

  END SUBROUTINE restore_upper_tri




end module gpu_routines
#endif
