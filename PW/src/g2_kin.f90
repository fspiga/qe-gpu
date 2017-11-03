!
! Copyright (C) 2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE g2_kin ( ik )
  !----------------------------------------------------------------------------
  !
  ! ... Calculation of kinetic energy - includes the case of the modified
  ! ... kinetic energy functional for variable-cell calculations
  !
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : tpiba2 
  USE gvecw,                ONLY : ecfixed, qcutz, q2sigma
#ifdef USE_CUDA
  USE klist,                ONLY : xk, ngk, igk_k=>igk_k_d
  USE wvfct,                ONLY : g2kin=>g2kin_d
  USE gvect,                ONLY : g=>g_d
#else
  USE klist,                ONLY : xk, ngk, igk_k
  USE wvfct,                ONLY : g2kin
  USE gvect,                ONLY : g
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT (IN) :: ik
  !
  ! ... local variables
  !
  INTEGER :: ig, npw,i
  REAL(DP):: xk1,xk2,xk3
#ifndef USE_CUDA
  REAL(DP), EXTERNAL :: qe_erf  !we use the erf intrinsic for the GPU path and the QE erf on the cpu
#endif
  !
  !
  npw = ngk(ik)

  xk1 = xk(1,ik)
  xk2 = xk(2,ik)
  xk3 = xk(3,ik)

#ifdef USE_CUDA
!$cuf kernel do(1) <<<*,*>>>
#else
!$omp parallel do
#endif
  do i=1,npw
  g2kin(i) = ( ( xk1 + g(1,igk_k(i,ik)) )*( xk1 + g(1,igk_k(i,ik)) ) + &
               ( xk2 + g(2,igk_k(i,ik)) )*( xk2 + g(2,igk_k(i,ik)) ) + &
               ( xk3 + g(3,igk_k(i,ik)) )*( xk3 + g(3,igk_k(i,ik)) ) ) * tpiba2
  !
  end do
#ifndef USE_CUDA
!$omp end parallel do
#endif

  IF ( qcutz > 0.D0 ) THEN
     !
#ifdef USE_CUDA
!$cuf kernel do(1) <<<*,*>>>
     DO ig = 1, npw
        !
        g2kin(ig) = g2kin(ig) + qcutz * &
             ( 1.D0 + erf( ( g2kin(ig) - ecfixed ) / q2sigma ) )
        !
     END DO
#else
!$omp parallel do
     DO ig = 1, npw
        !
        g2kin(ig) = g2kin(ig) + qcutz * &
             ( 1.D0 + qe_erf( ( g2kin(ig) - ecfixed ) / q2sigma ) )
        !
     END DO
!$omp end parallel do
#endif
     !
  END IF

  RETURN
  !
END SUBROUTINE g2_kin
