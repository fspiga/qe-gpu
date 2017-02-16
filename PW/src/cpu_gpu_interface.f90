
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

   INTERFACE cdiaghg
      SUBROUTINE cdiaghg_cpu( n, m, h, s, ldh, e, v )
         USE kinds,            ONLY : DP
         INTEGER, INTENT(IN) :: n, m, ldh
         COMPLEX(DP), INTENT(INOUT) :: h(:,:), s(:,:)
         REAL(DP), INTENT(OUT) :: e(:)
         COMPLEX(DP), INTENT(OUT) :: v(:,:)
      END SUBROUTINE cdiaghg_cpu

      SUBROUTINE cdiaghg_gpu( n, m, h, s, ldh, e, v )
         USE kinds,            ONLY : DP
         INTEGER, INTENT(IN) :: n, m, ldh
         COMPLEX(DP), INTENT(INOUT) :: h(:,:), s(:,:)
         REAL(DP), INTENT(OUT) :: e(:)
         COMPLEX(DP), INTENT(OUT) :: v(:,:)
         ATTRIBUTES( DEVICE ) :: h, s, e, v
      END SUBROUTINE cdiaghg_gpu
   END INTERFACE

   INTERFACE rotate_wfc_k
      SUBROUTINE rotate_wfc_k_cpu( npwx, npw, nstart, nbnd, npol, psi, overlap, evc, e )
         USE kinds,            ONLY : DP
         INTEGER :: npw, npwx, nstart, nbnd, npol
         LOGICAL :: overlap
         COMPLEX(DP) :: psi(:,:), evc(:,:)
         REAL(DP) :: e(nbnd)
      END SUBROUTINE rotate_wfc_k_cpu

      SUBROUTINE rotate_wfc_k_gpu( npwx, npw, nstart, nbnd, npol, psi, overlap, evc, e )
         USE kinds,            ONLY : DP
         INTEGER :: npw, npwx, nstart, nbnd, npol
         LOGICAL :: overlap
         COMPLEX(DP) :: psi(:,:), evc(:,:)
         REAL(DP) :: e(nbnd)
         ATTRIBUTES( DEVICE ) :: psi, evc, e
      END SUBROUTINE rotate_wfc_k_gpu
   END INTERFACE
!=----------------------------------------------------------------------=!
END MODULE cpu_gpu_interface
!=----------------------------------------------------------------------=!
#endif
