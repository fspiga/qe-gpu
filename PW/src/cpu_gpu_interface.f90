
! Copyright (C) 2002-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#ifdef USE_CUDA
!=----------------------------------------------------------------------=!
MODULE cpu_gpu_interface
!=----------------------------------------------------------------------=!
  USE cudafor
  IMPLICIT NONE

  INTERFACE add_vuspsi
     SUBROUTINE add_vuspsi_cpu( lda, n, m, hpsi )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m
        COMPLEX(DP) :: hpsi(:,:)
     END SUBROUTINE add_vuspsi_cpu

     SUBROUTINE add_vuspsi_gpu( lda, n, m, hpsi )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m
        COMPLEX(DP) :: hpsi(:,:)
        ATTRIBUTES( DEVICE ) :: hpsi
     END SUBROUTINE add_vuspsi_gpu
  END INTERFACE

  INTERFACE vloc_psi_k
     SUBROUTINE vloc_psi_k_cpu( lda, n, m, psi, v, hpsi  )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m
        COMPLEX(DP) :: psi(:,:), hpsi(:,:)
        REAL(DP) :: v(:)
     END SUBROUTINE vloc_psi_k_cpu

     SUBROUTINE vloc_psi_k_gpu( lda, n, m, psi, v, hpsi  )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m
        COMPLEX(DP) :: psi(:,:), hpsi(:,:)
        REAL(DP) :: v(:)
        ATTRIBUTES( DEVICE ) :: psi, v, hpsi
     END SUBROUTINE vloc_psi_k_gpu
  END INTERFACE

  INTERFACE h_psi
     SUBROUTINE h_psi_cpu( lda, n, m, psi, hpsi  )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m
!pgi$ ignore_tkr(r) psi,hpsi
        COMPLEX(DP) :: psi(*), hpsi(*)
     END SUBROUTINE h_psi_cpu

     SUBROUTINE h_psi_gpu( lda, n, m, psi, hpsi  )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m
!pgi$ ignore_tkr(r) psi,hpsi
        COMPLEX(DP) :: psi(*), hpsi(*)
        ATTRIBUTES( DEVICE ) :: psi, hpsi
     END SUBROUTINE h_psi_gpu
   END INTERFACE

  INTERFACE s_psi
     SUBROUTINE s_psi_cpu( lda, n, m, psi, spsi  )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m
!pgi$ ignore_tkr(r) psi,spsi
        COMPLEX(DP) :: psi(*), spsi(*)
     END SUBROUTINE s_psi_cpu

     SUBROUTINE s_psi_gpu( lda, n, m, psi, spsi  )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m
!pgi$ ignore_tkr(r) psi,spsi
        COMPLEX(DP) :: psi(*), spsi(*)
        ATTRIBUTES( DEVICE ) :: psi, spsi
     END SUBROUTINE s_psi_gpu
   END INTERFACE


!=----------------------------------------------------------------------=!
END MODULE cpu_gpu_interface
!=----------------------------------------------------------------------=!
#endif
