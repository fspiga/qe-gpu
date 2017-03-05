
! Copyright (C) 2002-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!=----------------------------------------------------------------------=!
MODULE cpu_gpu_interface
!=----------------------------------------------------------------------=!
#ifdef USE_CUDA
  USE cudafor
#endif
  IMPLICIT NONE

  INTERFACE add_vuspsi
     SUBROUTINE add_vuspsi_cpu( lda, n, m, hpsi )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m
        COMPLEX(DP) :: hpsi(:,:)
     END SUBROUTINE add_vuspsi_cpu

#ifdef USE_CUDA
     SUBROUTINE add_vuspsi_gpu( lda, n, m, hpsi )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m
        COMPLEX(DP) :: hpsi(:,:)
        ATTRIBUTES( DEVICE ) :: hpsi
     END SUBROUTINE add_vuspsi_gpu
#endif
  END INTERFACE

  INTERFACE vloc_psi_k
     SUBROUTINE vloc_psi_k_cpu( lda, n, m, psi, v, hpsi  )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m
        COMPLEX(DP) :: psi(:,:), hpsi(:,:)
        REAL(DP) :: v(:)
     END SUBROUTINE vloc_psi_k_cpu

#ifdef USE_CUDA
     SUBROUTINE vloc_psi_k_gpu( lda, n, m, psi, v, hpsi  )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m
        COMPLEX(DP) :: psi(:,:), hpsi(:,:)
        REAL(DP) :: v(:)
        ATTRIBUTES( DEVICE ) :: psi, v, hpsi
     END SUBROUTINE vloc_psi_k_gpu
#endif
  END INTERFACE

  INTERFACE h_psi
     SUBROUTINE h_psi_cpu( lda, n, m, psi, hpsi  )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m
!pgi$ ignore_tkr(r) psi,hpsi
        COMPLEX(DP) :: psi(*), hpsi(*)
     END SUBROUTINE h_psi_cpu

#ifdef USE_CUDA
     SUBROUTINE h_psi_gpu( lda, n, m, psi, hpsi  )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m
!pgi$ ignore_tkr(r) psi,hpsi
        COMPLEX(DP) :: psi(*), hpsi(*)
        ATTRIBUTES( DEVICE ) :: psi, hpsi
     END SUBROUTINE h_psi_gpu
#endif
   END INTERFACE

  INTERFACE s_psi
     SUBROUTINE s_psi_cpu( lda, n, m, psi, spsi  )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m
!pgi$ ignore_tkr(r) psi,spsi
        COMPLEX(DP) :: psi(*), spsi(*)
     END SUBROUTINE s_psi_cpu

#ifdef USE_CUDA
     SUBROUTINE s_psi_gpu( lda, n, m, psi, spsi  )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m
!pgi$ ignore_tkr(r) psi,spsi
        COMPLEX(DP) :: psi(*), spsi(*)
        ATTRIBUTES( DEVICE ) :: psi, spsi
     END SUBROUTINE s_psi_gpu
#endif
   END INTERFACE

  INTERFACE g_psi
     SUBROUTINE g_psi_cpu( lda, n, m, npol, psi, e  )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m, npol
!pgi$ ignore_tkr(r) psi,e
        COMPLEX(DP) :: psi(*)
        REAL(DP) :: e(*)
     END SUBROUTINE g_psi_cpu

#ifdef USE_CUDA
     SUBROUTINE g_psi_gpu( lda, n, m, npol, psi, e  )
        USE kinds, ONLY: DP
        INTEGER  :: lda, n, m, npol
!pgi$ ignore_tkr(r) psi,e
        COMPLEX(DP) :: psi(*)
        REAL(DP) :: e(*)
        ATTRIBUTES( DEVICE ) :: psi, e
     END SUBROUTINE g_psi_gpu
#endif
   END INTERFACE


   INTERFACE cdiaghg
      SUBROUTINE cdiaghg_cpu( n, m, h, s, ldh, e, v )
         USE kinds,            ONLY : DP
         INTEGER, INTENT(IN) :: n, m, ldh
         COMPLEX(DP), INTENT(INOUT) :: h(:,:), s(:,:)
         REAL(DP), INTENT(OUT) :: e(:)
         COMPLEX(DP), INTENT(OUT) :: v(:,:)
      END SUBROUTINE cdiaghg_cpu

#ifdef USE_CUDA
      SUBROUTINE cdiaghg_gpu( n, m, h, s, ldh, e, v )
         USE kinds,            ONLY : DP
         INTEGER, INTENT(IN) :: n, m, ldh
         COMPLEX(DP), INTENT(INOUT) :: h(:,:), s(:,:)
         REAL(DP), INTENT(OUT) :: e(:)
         COMPLEX(DP), INTENT(OUT) :: v(:,:)
         ATTRIBUTES( DEVICE ) :: h, s, e, v
      END SUBROUTINE cdiaghg_gpu
#endif
   END INTERFACE

   INTERFACE rotate_wfc_k
      SUBROUTINE rotate_wfc_k_cpu( npwx, npw, nstart, nbnd, npol, psi, overlap, evc, e )
         USE kinds,            ONLY : DP
         INTEGER :: npw, npwx, nstart, nbnd, npol
         LOGICAL :: overlap
         COMPLEX(DP) :: psi(:,:), evc(:,:)
         REAL(DP) :: e(nbnd)
      END SUBROUTINE rotate_wfc_k_cpu

#ifdef USE_CUDA
      SUBROUTINE rotate_wfc_k_gpu( npwx, npw, nstart, nbnd, npol, psi, overlap, evc, e )
         USE kinds,            ONLY : DP
         INTEGER :: npw, npwx, nstart, nbnd, npol
         LOGICAL :: overlap
         COMPLEX(DP) :: psi(:,:), evc(:,:)
         REAL(DP) :: e(nbnd)
         ATTRIBUTES( DEVICE ) :: psi, evc, e
      END SUBROUTINE rotate_wfc_k_gpu
#endif
   END INTERFACE

   INTERFACE rotate_wfc
      SUBROUTINE rotate_wfc_cpu( npwx, npw, nstart, gstart, nbnd, psi, npol, overlap, evc, e )
         USE kinds,            ONLY : DP
         INTEGER :: npw, npwx, nstart, nbnd, gstart, npol
         LOGICAL :: overlap
!pgi$ ignore_tkr(r) psi,evc,e
         COMPLEX(DP) :: psi(*), evc(*)
         REAL(DP) :: e(*)
      END SUBROUTINE rotate_wfc_cpu

#ifdef USE_CUDA
      SUBROUTINE rotate_wfc_gpu( npwx, npw, nstart, gstart, nbnd, psi, npol, overlap, evc, e )
         USE kinds,            ONLY : DP
         INTEGER :: npw, npwx, nstart, nbnd, gstart, npol
         LOGICAL :: overlap
!pgi$ ignore_tkr(r) psi,evc,e
         COMPLEX(DP) :: psi(*), evc(*)
         REAL(DP) :: e(*)
         ATTRIBUTES( DEVICE ) :: psi, evc, e
      END SUBROUTINE rotate_wfc_gpu
#endif
   END INTERFACE

!=----------------------------------------------------------------------=!
END MODULE cpu_gpu_interface
!=----------------------------------------------------------------------=!
