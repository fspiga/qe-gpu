!
! Copyright (C) 2009-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!---------------------------------------------------------------------------
SUBROUTINE compute_deff(deff, et)
!
!  This routine computes the effective value of the D-eS coefficients
!  which appear often in many expressions in the US or PAW case. 
!  This routine is for the collinear case.
!
USE kinds, ONLY : DP
USE ions_base, ONLY : nsp, nat, ityp
USE uspp, ONLY : deeq, qq, okvan
USE uspp_param, ONLY : nhm
USE lsda_mod, ONLY : current_spin
IMPLICIT NONE

INTEGER :: nt, na, is
REAL(DP), INTENT(OUT) :: deff(nhm, nhm, nat) 
REAL(DP), INTENT(IN) :: et

deff(:,:,:) = deeq(:,:,:,current_spin) 
IF (okvan) THEN
   DO nt = 1, nsp
      DO na = 1, nat
         IF ( ityp(na) == nt ) THEN
            deff(:,:,na) = deff(:,:,na) - et*qq(:,:,nt)
         END IF
      END DO
   END DO
ENDIF
RETURN
END SUBROUTINE compute_deff
!
#ifdef USE_CUDA
SUBROUTINE compute_deff_gpu(deff, et)
!
!  This routine computes the effective value of the D-eS coefficients
!  which appear often in many expressions in the US or PAW case. 
!  This routine is for the collinear case.
!
USE kinds, ONLY : DP
USE ions_base, ONLY : nsp, nat, ityp, ityp_d
USE uspp, ONLY : deeq_d, qq_d, okvan
USE uspp_param, ONLY : nhm
USE lsda_mod, ONLY : current_spin
IMPLICIT NONE

INTEGER :: nt, na, is
REAL(DP), DEVICE, INTENT(OUT) :: deff(nhm, nhm, nat)
REAL(DP), INTENT(IN) :: et
INTEGER :: i,j,k

!$cuf kernel do(3) <<<*,*>>>
DO k=1,nat
  DO j=1,nhm
    DO i=1,nhm
      deff(i,j,k) = deeq_d(i,j,k,current_spin)
    END DO
  END DO
END DO

IF (okvan) THEN
  !$cuf kernel do(3) <<<*,*>>>
  DO na = 1, nat
    DO j=1,nhm
      DO i=1,nhm
         deff(i,j,na) = deff(i,j,na) - et*qq_d(i,j,ityp_d(na))
      END DO
    END DO
  END DO
ENDIF
RETURN
END SUBROUTINE compute_deff_gpu
#endif
!
SUBROUTINE compute_deff_nc(deff, et)
!
!  This routine computes the effective value of the D-eS coefficients
!  which appears often in many expressions. This routine is for the
!  noncollinear case.
!
USE kinds, ONLY : DP
USE ions_base, ONLY : nsp, nat, ityp
USE spin_orb,  ONLY : lspinorb
USE noncollin_module, ONLY : noncolin, npol
USE uspp, ONLY : deeq_nc, qq, qq_so, okvan
USE uspp_param, ONLY : nhm
USE lsda_mod, ONLY : nspin
IMPLICIT NONE

INTEGER :: nt, na, is, js, ijs
COMPLEX(DP), INTENT(OUT) :: deff(nhm, nhm, nat, nspin) 
REAL(DP), INTENT(IN) :: et

deff=deeq_nc
IF (okvan) THEN
   DO nt = 1, nsp
      DO na = 1, nat
         IF ( ityp(na) == nt ) THEN
            IF (lspinorb) THEN
               deff(:,:,na,:) = deff(:,:,na,:) - et * qq_so(:,:,:,nt)
            ELSE
               ijs=0
               DO is=1,npol
                  DO js=1,npol
                     ijs=ijs+1
                     IF (is==js) deff(:,:,na,ijs)=deff(:,:,na,ijs)-et*qq(:,:,nt)
                  END DO
               END DO
            END IF
         END IF
      END DO
   END DO
ENDIF

RETURN
END SUBROUTINE compute_deff_nc

