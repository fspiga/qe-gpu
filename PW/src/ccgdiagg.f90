  !
  ! Copyright (C) 2001-2007 Quantum ESPRESSO group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  !
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
  ! define __VERBOSE to print a message after each eigenvalue is computed
  !
#ifdef USE_CUDA
  SUBROUTINE cgDdot( n, A, B, res )
    !
    USE kinds, ONLY : DP
    USE cudafor
    USE cublas
    !
    IMPLICIT NONE
    !
    INTEGER, INTENT(IN) :: n
    !
!!!!!pgi$ ignore_tkr A
    REAL(DP), DEVICE, INTENT(IN) :: A(n), B(n)
    !
    REAL(DP), INTENT(OUT) :: res
    !
    res = cublasDdot( n, A, 1, B, 1 )
    !
    RETURN
    !
  END SUBROUTINE cgDdot
#endif

  !----------------------------------------------------------------------------
#ifdef USE_CUDA
  SUBROUTINE ccgdiagg( npwx, npw, nbnd, npol, psi, psi_d, e, e_d, btype, &
       precondition, ethr, maxter, reorder, notconv, avg_iter )
#else
    SUBROUTINE ccgdiagg( npwx, npw, nbnd, npol, psi, e, btype, precondition, &
         ethr, maxter, reorder, notconv, avg_iter )
#endif
      !----------------------------------------------------------------------------
      !
      ! ... "poor man" iterative diagonalization of a complex hermitian matrix
      ! ... through preconditioned conjugate gradient algorithm
      ! ... Band-by-band algorithm with minimal use of memory
      ! ... Calls h_1psi and s_1psi to calculate H|psi> and S|psi>
      ! ... Works for generalized eigenvalue problem (US pseudopotentials) as well
      !
      USE constants,        ONLY : pi
      USE kinds,            ONLY : DP
      USE mp_bands,         ONLY : intra_bgrp_comm
      USE mp,               ONLY : mp_sum
      USE cpu_gpu_interface, ONLY : s_psi, h_psi
#if defined(__VERBOSE)
      USE io_global, ONLY : stdout
#endif
      !
#ifdef USE_CUDA
      USE cudafor
      USE gpu_routines
#endif
      !
      IMPLICIT NONE
      !
      ! ... I/O variables
      !
      INTEGER,     INTENT(IN)    :: npwx, npw, nbnd, npol, maxter
      INTEGER,     INTENT(IN)    :: btype(nbnd)
      REAL(DP),    INTENT(IN)    :: precondition(npwx*npol), ethr
      COMPLEX(DP), INTENT(INOUT) :: psi(npwx*npol,nbnd)
      REAL(DP),    INTENT(INOUT) :: e(nbnd)
      INTEGER,     INTENT(OUT)   :: notconv
      REAL(DP),    INTENT(OUT)   :: avg_iter
#ifdef USE_CUDA
      COMPLEX(DP), DEVICE, INTENT(INOUT) ::psi_d(npwx*npol,nbnd)
      REAL(DP), DEVICE, INTENT(INOUT) :: e_d(nbnd)
#endif
      !
      ! ... local variables
      !
      INTEGER                  :: i, j, m, iter, moved
      COMPLEX(DP), ALLOCATABLE :: hpsi(:), spsi(:), lagrange(:), &
           g(:), cg(:), scg(:), ppsi(:), g0(:)
      REAL(DP)                 :: psi_norm, a0, b0, gg0, gamma, gg, gg1, &
           cg0, e0, es(2)
      REAL(DP)                 :: theta, cost, sint, cos2t, sin2t
      LOGICAL                  :: reorder
      INTEGER                  :: kdim, kdmx, kdim2, istat, ierr
      REAL(DP)                 :: empty_ethr, ethr_m, ddot_temp
      !
#ifdef USE_CUDA
      COMPLEX(DP),DEVICE, ALLOCATABLE :: hpsi_d(:), spsi_d(:), lagrange_d(:), &
           g_d(:), cg_d(:), scg_d(:), ppsi_d(:), g0_d(:)
      REAL(DP), DEVICE, ALLOCATABLE  :: precondition_d(:)
      REAL(DP), DEVICE  :: psi_norm_d, es_d, gamma_d, cost_d, sint_d, cg0_d
#endif
      !
      ! ... external functions
      !
#ifndef  USE_CUDA
      REAL (DP), EXTERNAL :: ddot
#endif
      !
      !
      CALL start_clock( 'ccgdiagg' )
      !
      empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
      !
      IF ( npol == 1 ) THEN
         !
         kdim = npw
         kdmx = npwx
         !
      ELSE
         !
         kdim = npwx * npol
         kdmx = npwx * npol
         !
      END IF
      !
      kdim2 = 2 * kdim
      !
      ALLOCATE( spsi( kdmx ) )
      ALLOCATE( scg(  kdmx ) )
      ALLOCATE( hpsi( kdmx ) )
      ALLOCATE( g(    kdmx ) )
      ALLOCATE( cg(   kdmx ) )
      ALLOCATE( g0(   kdmx ) )
      ALLOCATE( ppsi( kdmx ) )
      !
      ALLOCATE( lagrange( nbnd ) )
      !
#ifdef USE_CUDA
      ALLOCATE(  hpsi_d(kdmx), STAT=ierr )
      IF( ierr /= 0 ) &
      CALL errore( ' ccgdiagg ',' cannot allocate hpsi_d ', ABS(ierr) )
      ALLOCATE(  spsi_d(kdmx), STAT=ierr )
      IF( ierr /= 0 ) &
      CALL errore( ' ccgdiagg ',' cannot allocate spsi_d ', ABS(ierr) )
      ALLOCATE(  g_d(kdmx), STAT=ierr )
      IF( ierr /= 0 ) &
      CALL errore( ' ccgdiagg ',' cannot allocate g_d ', ABS(ierr) )
      ALLOCATE(  cg_d(kdmx), STAT=ierr )
      IF( ierr /= 0 ) &
      CALL errore( ' ccgdiagg ',' cannot allocate cg_d ', ABS(ierr) )
      ALLOCATE(  scg_d(kdmx), STAT=ierr )
      IF( ierr /= 0 ) &
      CALL errore( ' ccgdiagg ',' cannot allocate scg_d ', ABS(ierr) )
      ALLOCATE(  ppsi_d(kdmx), STAT=ierr )
      IF( ierr /= 0 ) &
      CALL errore( ' ccgdiagg ',' cannot allocate ppsi_d ', ABS(ierr) )
      ALLOCATE(  g0_d(kdmx), STAT=ierr )
      IF( ierr /= 0 ) &
      CALL errore( ' ccgdiagg ',' cannot allocate g0_d ', ABS(ierr) )
      ALLOCATE(  precondition_d(kdmx), STAT=ierr )
      IF( ierr /= 0 ) &
      CALL errore( ' ccgdiagg ',' cannot allocate precondition_d ', ABS(ierr) )
      !
      ALLOCATE(  lagrange_d(nbnd), STAT=ierr )
      IF( ierr /= 0 ) &
      CALL errore( ' ccgdiagg ',' cannot allocate lagrange_d ', ABS(ierr) )
#endif
      avg_iter = 0.D0
      notconv  = 0
      moved    = 0
#ifdef USE_CUDA
      precondition_d = precondition
#endif
      !
      ! ... every eigenfunction is calculated separately
      !
      cgiter: DO m = 1, nbnd
         !
         !
         IF ( btype(m) == 1 ) THEN
            !
            ethr_m = ethr
            !
         ELSE
            !
            ethr_m = empty_ethr
            !
         END IF
         !
         spsi     = ZERO
         scg      = ZERO
         hpsi     = ZERO
         g        = ZERO
         cg       = ZERO
         g0       = ZERO
         ppsi     = ZERO
         lagrange = ZERO
#ifdef USE_CUDA
         spsi_d     = ZERO
         scg_d      = ZERO
         hpsi_d     = ZERO
         g_d        = ZERO
         cg_d       = ZERO
         g0_d       = ZERO
         ppsi_d     = ZERO
         lagrange_d = ZERO
#endif
         !
         ! ... calculate S|psi>
         !
#ifdef USE_CUDA
        !  CALL s_1psi( npwx, npw, psi(1,m), spsi ) !Original
         CALL s_psi(npwx, npw, 1, psi_d(1,m), spsi_d)
         !
         ! ... orthogonalize starting eigenfunction to those already calculated
         !
         istat = cublasZgemv(cublasH, 2, kdim, m, ONE, &
              psi_d, kdmx, spsi_d, 1, ZERO, lagrange_d, 1 )
         !
         CALL mp_sum( lagrange_d( 1:m ), intra_bgrp_comm )
         lagrange = lagrange_d
         !
         psi_norm_d = DBLE( lagrange(m) )
         psi_norm = psi_norm_d
         !
         DO j = 1, m - 1
            !
!$cuf kernel do(1) <<<*,*>>>
            DO i = 1, kdmx
               psi_d(i,m)  = psi_d(i,m) - lagrange_d(j) * psi_d(i,j)
            END DO
            !
            psi_norm = psi_norm - &
                 ( DBLE( lagrange(j) )**2 + AIMAG( lagrange(j) )**2 )
            !
         END DO
         !
         psi_norm = SQRT( psi_norm )
         psi_norm_d = psi_norm
         !
!$cuf kernel do(1) <<<*,*>>>
         DO i = 1, kdmx
           psi_d(i,m) = psi_d(i,m) / psi_norm_d
         END DO
         !
         ! ... calculate starting gradient (|hpsi> = H|psi>) ...
         !
         !FIX============================================
         !  CALL h_1psi( npwx, npw, psi(1,m), hpsi, spsi )
         CALL h_psi(npwx, npw, 1, psi_d(1,m), hpsi_d)
         CALL s_psi(npwx, npw, 1, psi_d(1,m), spsi_d)
         !===============================================
         !
         ! ... and starting eigenvalue (e = <y|PHP|y> = <psi|H|psi>)
         !
         ! ... NB:  ddot(2*npw,a,1,b,1) = REAL( zdotc(npw,a,1,b,1) )
         !
         CALL cgddot( kdim2, psi_d(1,m), hpsi_d, ddot_temp )
         e(m) = ddot_temp
         !
         CALL mp_sum( e(m), intra_bgrp_comm )
         !
         ! ... start iteration for this band
         !
         iterate: DO iter = 1, maxter
            !
            ! ... calculate  P (PHP)|y>
            ! ... ( P = preconditioning matrix, assumed diagonal )
            !
!$cuf kernel do(1) <<<*,*>>>
            DO i = 1, kdmx
               g_d(i)    = hpsi_d(i) / precondition_d(i)
               ppsi_d(i) = spsi_d(i) / precondition_d(i)
            END DO
            !
            ! ... ppsi is now S P(P^2)|y> = S P^2|psi>)
            !
            CALL cgddot( kdim2, spsi_d(1), g_d(1), ddot_temp)
            es(1) = ddot_temp
            CALL cgddot( kdim2, spsi_d(1), ppsi_d(1), ddot_temp)
            es(2) = ddot_temp
            !
            CALL mp_sum( es , intra_bgrp_comm )
            !
            es(1) = es(1) / es(2)
            es_d = es(1)
            !
!$cuf kernel do(1) <<<*,*>>>
            DO i = 1, kdmx
               g_d(i) = g_d(i) - es_d * ppsi_d(i)
            END DO
            !
            ! ... e1 = <y| S P^2 PHP|y> / <y| S S P^2|y> ensures that
            ! ... <g| S P^2|y> = 0
            ! ... orthogonalize to lowest eigenfunctions (already calculated)
            !
            ! ... scg is used as workspace
            !
            ! CALL s_1psi( npwx, npw, g(1), scg(1) ) !original
            CALL s_psi(npwx, npw, 1, g_d(1), scg_d(1))
            scg = scg_d
            !
            istat = cublasZgemv(cublasH, 2, kdim, (m -1) , ONE, &
                      psi_d, kdmx, scg_d, 1, ZERO, lagrange_d, 1 )
            !
            CALL mp_sum( lagrange_d( 1:m-1 ), intra_bgrp_comm )
            !
            DO j = 1, ( m - 1 )
               !
!$cuf kernel do(1) <<<*,*>>>
               DO i = 1, kdmx
                  g_d(i)   = g_d(i)   - lagrange_d(j) * psi_d(i,j)
                  scg_d(i) = scg_d(i) - lagrange_d(j) * psi_d(i,j)
               END DO
               !
            END DO
            !
            IF ( iter /= 1 ) THEN
               !
               ! ... gg1 is <g(n+1)|S|g(n)> (used in Polak-Ribiere formula)
               !
               g0_d = g0
               CALL cgddot( kdim2, g_d(1), g0_d(1), ddot_temp )
               gg1 = ddot_temp
               !
               CALL mp_sum( gg1, intra_bgrp_comm )
               !
            END IF
            !
            ! ... gg is <g(n+1)|S|g(n+1)>
            !
!$cuf kernel do(1) <<<*,*>>>
            DO i = 1, kdmx
               g0_d(i) = scg_d(i)
            END DO
            !
!$cuf kernel do(1) <<<*,*>>>
            DO i = 1, kdmx
               g0_d(i) = g0_d(i) * precondition_d(i)
            END DO
            !
            CALL cgddot( kdim2, g_d(1), g0_d(1), ddot_temp )
            gg = ddot_temp
            !
            CALL mp_sum( gg, intra_bgrp_comm )
            !
            IF ( iter == 1 ) THEN
               !
               ! ... starting iteration, the conjugate gradient |cg> = |g>
               !
               gg0 = gg
               !
!$cuf kernel do(1) <<<*,*>>>
               DO i = 1, kdmx
                  cg_d(i) = g_d(i)
               END DO
               !
            ELSE
               !
               ! ... |cg(n+1)> = |g(n+1)> + gamma(n) * |cg(n)>
               !
               ! ... Polak-Ribiere formula :
               !
               gamma = ( gg - gg1 ) / gg0
               gg0   = gg
               gamma_d = gamma
               !
!$cuf kernel do(1) <<<*,*>>>
               DO i = 1, kdmx
                  cg_d(i) = cg_d(i) * gamma_d
               END DO
!$cuf kernel do(1) <<<*,*>>>
               DO i = 1, kdmx
                  cg_d(i) = g_d(i) + cg_d(i)
               END DO
               !
               ! ... The following is needed because <y(n+1)| S P^2 |cg(n+1)>
               ! ... is not 0. In fact :
               ! ... <y(n+1)| S P^2 |cg(n)> = sin(theta)*<cg(n)|S|cg(n)>
               !
               psi_norm = gamma * cg0 * sint
               psi_norm_d = psi_norm
               !
!$cuf kernel do(1) <<<*,*>>>
               DO i = 1, kdmx
                  cg_d(i) = cg_d(i) - psi_norm_d * psi_d(i,m)
               END DO
               !
            END IF
            !
            ! ... |cg> contains now the conjugate gradient
            !
            ! ... |scg> is S|cg>
            !
            !FIX==============================================
            ! CALL h_1psi( npwx, npw, cg(1), ppsi(1), scg(1) )
            CALL h_psi(npwx, npw, 1, cg_d(1), ppsi_d(1))
            CALL s_psi(npwx, npw, 1, cg_d(1), scg_d(1))
            !================================================
            !
            CALL cgddot( kdim2, cg_d(1), scg_d(1), ddot_temp )
            cg0 = ddot_temp
            !
            CALL mp_sum(  cg0 , intra_bgrp_comm )
            !
            cg0 = SQRT( cg0 )
            !
            ! ... |ppsi> contains now HP|cg>
            ! ... minimize <y(t)|PHP|y(t)> , where :
            ! ...                         |y(t)> = cos(t)|y> + sin(t)/cg0 |cg>
            ! ... Note that  <y|P^2S|y> = 1, <y|P^2S|cg> = 0 ,
            ! ...           <cg|P^2S|cg> = cg0^2
            ! ... so that the result is correctly normalized :
            ! ...                           <y(t)|P^2S|y(t)> = 1
            !
            CALL cgddot( kdim2, psi_d(1,m), ppsi_d(1), ddot_temp)
            a0 = 2.D0 * ddot_temp / cg0
            !
            CALL mp_sum(  a0 , intra_bgrp_comm )
            !
            CALL cgddot( kdim2, cg_d(1), ppsi_d(1), ddot_temp )
            b0 = ddot_temp / cg0**2
            !
            CALL mp_sum(  b0 , intra_bgrp_comm )
            !
            e0 = e(m)
            !
            theta = 0.5D0 * ATAN( a0 / ( e0 - b0 ) )
            !
            cost = COS( theta )
            sint = SIN( theta )
            !
            cos2t = cost*cost - sint*sint
            sin2t = 2.D0*cost*sint
            !
            es(1) = 0.5D0 * (   ( e0 - b0 ) * cos2t + a0 * sin2t + e0 + b0 )
            es(2) = 0.5D0 * ( - ( e0 - b0 ) * cos2t - a0 * sin2t + e0 + b0 )
            !
            ! ... there are two possible solutions, choose the minimum
            !
            IF ( es(2) < es(1) ) THEN
               !
               theta = theta + 0.5D0 * pi
               !
               cost = COS( theta )
               sint = SIN( theta )
               !
            END IF
            !
            ! ... new estimate of the eigenvalue
            !
            e(m) = MIN( es(1), es(2) )
            !
            ! ... upgrade |psi>
            !
            cost_d = cost
            sint_d = sint
            cg0_d = cg0
!$cuf kernel do(1) <<<*,*>>>
            DO i = 1, kdmx
               psi_d(i,m) = cost_d * psi_d(i,m) + sint_d / cg0_d * cg_d(i)
            END DO
            !
            ! ... here one could test convergence on the energy
            !
            IF ( ABS( e(m) - e0 ) < ethr_m ) EXIT iterate
            !
            ! ... upgrade H|psi> and S|psi>
            !
!$cuf kernel do(1) <<<*,*>>>
            DO i = 1, kdmx
               spsi_d(i) = cost_d * spsi_d(i) + sint_d / cg0_d * scg_d(i)
            END DO
            !
!$cuf kernel do(1) <<<*,*>>>
            DO i = 1, kdmx
               hpsi_d(i) = cost_d * hpsi_d(i) + sint_d / cg0_d * ppsi_d(i)
            END DO
            !
         END DO iterate
         !
#if defined(__VERBOSE)
         IF ( iter >= maxter ) THEN
            WRITE(stdout,'("e(",i4,") = ",f12.6," eV  (not converged after ",i3,&
                 & " iterations)")') m, e(m)*13.6058, iter
         ELSE
            WRITE(stdout,'("e(",i4,") = ",f12.6," eV  (",i3," iterations)")') &
                 m, e(m)*13.6058, iter
         END IF
         FLUSH (stdout)
#endif
         IF ( iter >= maxter ) notconv = notconv + 1
         !
         avg_iter = avg_iter + iter + 1
         !
         ! ... reorder eigenvalues if they are not in the right order
         ! ... ( this CAN and WILL happen in not-so-special cases )
         !
         IF ( m > 1 .AND. reorder ) THEN
            !
            IF ( e(m) - e(m-1) < - 2.D0 * ethr_m ) THEN
               !
               ! ... if the last calculated eigenvalue is not the largest...
               !
               DO i = m - 2, 1, - 1
                  !
                  IF ( e(m) - e(i) > 2.D0 * ethr_m ) EXIT
                  !
               END DO
               !
               i = i + 1
               !
               moved = moved + 1
               !
               ! ... last calculated eigenvalue should be in the
               ! ... i-th position: reorder
               !
               e0 = e(m)
               !
!$cuf kernel do(1) <<<*,*>>>
               DO i = 1, kdmx
                  ppsi_d(i) = psi_d(i,m)
               END DO
               !
               DO j = m, i + 1, - 1
                  !
                  e(j) = e(j-1)
                  !
!$cuf kernel do(1) <<<*,*>>>
                  DO i = 1, kdmx
                     psi_d(i,j) = psi_d(i,j-1)
                  END DO
                  !
               END DO
               !
               e(i) = e0
               !
!$cuf kernel do(1) <<<*,*>>>
               DO i = 1, kdmx
                  psi_d(i,i) = ppsi_d(i)
               END DO
               !
               ! ... this procedure should be good if only a few inversions occur,
               ! ... extremely inefficient if eigenvectors are often in bad order
               ! ... ( but this should not happen )
               !
            END IF
            !
         END IF
         !
         ! FIX : NOT Sure about this COPY's
         psi = psi_d
         e_d = e
#else
         CALL s_1psi( npwx, npw, psi(1,m), spsi )
         !
         ! ... orthogonalize starting eigenfunction to those already calculated
         !
         CALL ZGEMV( 'C', kdim, m, ONE, psi, kdmx, spsi, 1, ZERO, lagrange, 1 )
         !
         CALL mp_sum( lagrange( 1:m ), intra_bgrp_comm )
         !
         psi_norm = DBLE( lagrange(m) )
         !
         DO j = 1, m - 1
            !
            psi(:,m)  = psi(:,m) - lagrange(j) * psi(:,j)
            !
            psi_norm = psi_norm - &
                 ( DBLE( lagrange(j) )**2 + AIMAG( lagrange(j) )**2 )
            !
         END DO
         !
         psi_norm = SQRT( psi_norm )
         !
         psi(:,m) = psi(:,m) / psi_norm
         !
         ! ... calculate starting gradient (|hpsi> = H|psi>) ...
         !
         CALL h_1psi( npwx, npw, psi(1,m), hpsi, spsi )
         !
         ! ... and starting eigenvalue (e = <y|PHP|y> = <psi|H|psi>)
         !
         ! ... NB:  ddot(2*npw,a,1,b,1) = REAL( zdotc(npw,a,1,b,1) )
         !
         e(m) = ddot( kdim2, psi(1,m), 1, hpsi, 1 )
         !
         CALL mp_sum( e(m), intra_bgrp_comm )
         !
         ! ... start iteration for this band
         !
         iterate: DO iter = 1, maxter
            !
            ! ... calculate  P (PHP)|y>
            ! ... ( P = preconditioning matrix, assumed diagonal )
            !
            g(:)    = hpsi(:) / precondition(:)
            ppsi(:) = spsi(:) / precondition(:)
            !
            ! ... ppsi is now S P(P^2)|y> = S P^2|psi>)
            !
            es(1) = ddot( kdim2, spsi(1), 1, g(1), 1 )
            es(2) = ddot( kdim2, spsi(1), 1, ppsi(1), 1 )
            !
            CALL mp_sum( es , intra_bgrp_comm )
            !
            es(1) = es(1) / es(2)
            !
            g(:) = g(:) - es(1) * ppsi(:)
            !
            ! ... e1 = <y| S P^2 PHP|y> / <y| S S P^2|y> ensures that
            ! ... <g| S P^2|y> = 0
            ! ... orthogonalize to lowest eigenfunctions (already calculated)
            !
            ! ... scg is used as workspace
            !
            CALL s_1psi( npwx, npw, g(1), scg(1) )
            !
            CALL ZGEMV( 'C', kdim, ( m - 1 ), ONE, psi, &
                 kdmx, scg, 1, ZERO, lagrange, 1  )
            !
            CALL mp_sum( lagrange( 1:m-1 ), intra_bgrp_comm )
            !
            DO j = 1, ( m - 1 )
               !
               g(:)   = g(:)   - lagrange(j) * psi(:,j)
               scg(:) = scg(:) - lagrange(j) * psi(:,j)
               !
            END DO
            !
            IF ( iter /= 1 ) THEN
               !
               ! ... gg1 is <g(n+1)|S|g(n)> (used in Polak-Ribiere formula)
               !
               gg1 = ddot( kdim2, g(1), 1, g0(1), 1 )
               !
               CALL mp_sum( gg1, intra_bgrp_comm )
               !
            END IF
            !
            ! ... gg is <g(n+1)|S|g(n+1)>
            !
            g0(:) = scg(:)
            !
            g0(:) = g0(:) * precondition(:)
            !
            gg = ddot( kdim2, g(1), 1, g0(1), 1 )
            !
            CALL mp_sum( gg, intra_bgrp_comm )
            !
            IF ( iter == 1 ) THEN
               !
               ! ... starting iteration, the conjugate gradient |cg> = |g>
               !
               gg0 = gg
               !
               cg(:) = g(:)
               !
            ELSE
               !
               ! ... |cg(n+1)> = |g(n+1)> + gamma(n) * |cg(n)>
               !
               ! ... Polak-Ribiere formula :
               !
               gamma = ( gg - gg1 ) / gg0
               gg0   = gg
               !
               cg(:) = cg(:) * gamma
               cg(:) = g + cg(:)
               !
               ! ... The following is needed because <y(n+1)| S P^2 |cg(n+1)>
               ! ... is not 0. In fact :
               ! ... <y(n+1)| S P^2 |cg(n)> = sin(theta)*<cg(n)|S|cg(n)>
               !
               psi_norm = gamma * cg0 * sint
               !
               cg(:) = cg(:) - psi_norm * psi(:,m)
               !
            END IF
            !
            ! ... |cg> contains now the conjugate gradient
            !
            ! ... |scg> is S|cg>
            !
            CALL h_1psi( npwx, npw, cg(1), ppsi(1), scg(1) )
            !
            cg0 = ddot( kdim2, cg(1), 1, scg(1), 1 )
            !
            CALL mp_sum(  cg0 , intra_bgrp_comm )
            !
            cg0 = SQRT( cg0 )
            !
            ! ... |ppsi> contains now HP|cg>
            ! ... minimize <y(t)|PHP|y(t)> , where :
            ! ...                         |y(t)> = cos(t)|y> + sin(t)/cg0 |cg>
            ! ... Note that  <y|P^2S|y> = 1, <y|P^2S|cg> = 0 ,
            ! ...           <cg|P^2S|cg> = cg0^2
            ! ... so that the result is correctly normalized :
            ! ...                           <y(t)|P^2S|y(t)> = 1
            !
            a0 = 2.D0 * ddot( kdim2, psi(1,m), 1, ppsi(1), 1 ) / cg0
            !
            CALL mp_sum(  a0 , intra_bgrp_comm )
            !
            b0 = ddot( kdim2, cg(1), 1, ppsi(1), 1 ) / cg0**2
            !
            CALL mp_sum(  b0 , intra_bgrp_comm )
            !
            e0 = e(m)
            !
            theta = 0.5D0 * ATAN( a0 / ( e0 - b0 ) )
            !
            cost = COS( theta )
            sint = SIN( theta )
            !
            cos2t = cost*cost - sint*sint
            sin2t = 2.D0*cost*sint
            !
            es(1) = 0.5D0 * (   ( e0 - b0 ) * cos2t + a0 * sin2t + e0 + b0 )
            es(2) = 0.5D0 * ( - ( e0 - b0 ) * cos2t - a0 * sin2t + e0 + b0 )
            !
            ! ... there are two possible solutions, choose the minimum
            !
            IF ( es(2) < es(1) ) THEN
               !
               theta = theta + 0.5D0 * pi
               !
               cost = COS( theta )
               sint = SIN( theta )
               !
            END IF
            !
            ! ... new estimate of the eigenvalue
            !
            e(m) = MIN( es(1), es(2) )
            !
            ! ... upgrade |psi>
            !
            psi(:,m) = cost * psi(:,m) + sint / cg0 * cg(:)
            !
            ! ... here one could test convergence on the energy
            !
            IF ( ABS( e(m) - e0 ) < ethr_m ) EXIT iterate
            !
            ! ... upgrade H|psi> and S|psi>
            !
            spsi(:) = cost * spsi(:) + sint / cg0 * scg(:)
            !
            hpsi(:) = cost * hpsi(:) + sint / cg0 * ppsi(:)
            !
         END DO iterate
         !
#if defined(__VERBOSE)
         IF ( iter >= maxter ) THEN
            WRITE(stdout,'("e(",i4,") = ",f12.6," eV  (not converged after ",i3,&
                 & " iterations)")') m, e(m)*13.6058, iter
         ELSE
            WRITE(stdout,'("e(",i4,") = ",f12.6," eV  (",i3," iterations)")') &
                 m, e(m)*13.6058, iter
         END IF
         FLUSH (stdout)
#endif
         IF ( iter >= maxter ) notconv = notconv + 1
         !
         avg_iter = avg_iter + iter + 1
         !
         ! ... reorder eigenvalues if they are not in the right order
         ! ... ( this CAN and WILL happen in not-so-special cases )
         !
         IF ( m > 1 .AND. reorder ) THEN
            !
            IF ( e(m) - e(m-1) < - 2.D0 * ethr_m ) THEN
               !
               ! ... if the last calculated eigenvalue is not the largest...
               !
               DO i = m - 2, 1, - 1
                  !
                  IF ( e(m) - e(i) > 2.D0 * ethr_m ) EXIT
                  !
               END DO
               !
               i = i + 1
               !
               moved = moved + 1
               !
               ! ... last calculated eigenvalue should be in the
               ! ... i-th position: reorder
               !
               e0 = e(m)
               !
               ppsi(:) = psi(:,m)
               !
               DO j = m, i + 1, - 1
                  !
                  e(j) = e(j-1)
                  !
                  psi(:,j) = psi(:,j-1)
                  !
               END DO
               !
               e(i) = e0
               !
               psi(:,i) = ppsi(:)
               !
               ! ... this procedure should be good if only a few inversions occur,
               ! ... extremely inefficient if eigenvectors are often in bad order
               ! ... ( but this should not happen )
               !
            END IF
            !
         END IF
         !
#endif
      END DO cgiter
      !
      avg_iter = avg_iter / DBLE( nbnd )
      !
      DEALLOCATE( lagrange )
      DEALLOCATE( ppsi )
      DEALLOCATE( g0 )
      DEALLOCATE( cg )
      DEALLOCATE( g )
      DEALLOCATE( hpsi )
      DEALLOCATE( scg )
      DEALLOCATE( spsi )
      !
#ifdef USE_CUDA
      DEALLOCATE( lagrange_d )
      DEALLOCATE( ppsi_d )
      DEALLOCATE( g0_d )
      DEALLOCATE( cg_d )
      DEALLOCATE( g_d )
      DEALLOCATE( hpsi_d )
      DEALLOCATE( scg_d )
      DEALLOCATE( spsi_d )
#endif
      CALL stop_clock( 'ccgdiagg' )
      !
      RETURN
      !
    END SUBROUTINE ccgdiagg
