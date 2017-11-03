!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE stres_us( ik, gk, sigmanlc )
  !----------------------------------------------------------------------------
  !
  ! nonlocal (separable pseudopotential) contribution to the stress
  ! NOTICE: sum of partial results over procs is performed in calling routine
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE constants,            ONLY : eps8
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE wvfct,                ONLY : npwx, nbnd, wg, et
  USE control_flags,        ONLY : gamma_only
  USE uspp_param,           ONLY : upf, lmaxkb, nh, newpseudo, nhm
  USE uspp,                 ONLY : nkb, vkb, qq, deeq, deeq_nc, qq_so
#ifdef USE_CUDA
  USE uspp,                 ONLY : vkb_d
#endif
  USE wavefunctions_module, ONLY : evc
  USE spin_orb,             ONLY : lspinorb
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : noncolin, npol
  USE mp_pools,             ONLY : me_pool, root_pool
  USE mp_bands,             ONLY : intra_bgrp_comm, me_bgrp, root_bgrp
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, &
                                   bec_type, becp, calbec
  USE mp,                   ONLY : mp_sum, mp_get_comm_null, mp_circular_shift_left 
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)    :: ik
  REAL(DP), INTENT(IN)   :: gk(3,npwx)
  REAL(DP), INTENT(INOUT):: sigmanlc(3,3)
  !
  REAL(DP), ALLOCATABLE  :: qm1(:)
  REAL(DP)               :: q
  INTEGER                :: npw, i
  !
  !
  IF ( nkb == 0 ) RETURN
  !
  IF ( lsda ) current_spin = isk(ik)
  npw = ngk(ik)

#ifdef USE_CUDA
  vkb = vkb_d
#endif

  IF ( nks > 1 ) CALL init_us_2( npw, igk_k(1,ik), xk(1,ik), vkb )
  !
  CALL allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm ) 
  CALL calbec( npw, vkb, evc, becp )
  !
  ALLOCATE( qm1( npwx ) )
  DO i = 1, npw
     q = SQRT( gk(1,i)**2 + gk(2,i)**2 + gk(3,i)**2 )
     IF ( q > eps8 ) THEN
        qm1(i) = 1.D0 / q
     ELSE
        qm1(i) = 0.D0
     END IF
  END DO
  !
  IF ( gamma_only ) THEN
     !
     CALL stres_us_gamma()
     !
  ELSE
     !
     CALL stres_us_k()
     !
  END IF
  !
  DEALLOCATE( qm1 )
  CALL deallocate_bec_type ( becp ) 
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE stres_us_gamma()
       !-----------------------------------------------------------------------
       ! 
       ! ... gamma version
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER                       :: na, np, ibnd, ipol, jpol, l, i, &
                                        ikb, jkb, ih, jh, ijkb0, ibnd_loc, &
                                        nproc, nbnd_loc, nbnd_begin, icyc
       INTEGER, EXTERNAL :: ldim_block, lind_block, gind_block
       REAL(DP)                 :: fac, xyz(3,3), evps, ddot
       REAL(DP), ALLOCATABLE    :: deff(:,:,:)
       COMPLEX(DP), ALLOCATABLE :: work1(:), work2(:), dvkb(:,:)
       ! dvkb contains the derivatives of the kb potential
       COMPLEX(DP)              :: ps
       ! xyz are the three unit vectors in the x,y,z directions
       DATA xyz / 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0 /
       !
       !
       IF( becp%comm /= mp_get_comm_null() ) THEN
          nproc   = becp%nproc
          nbnd_loc   = becp%nbnd_loc
          nbnd_begin = becp%ibnd_begin
          IF( ( nbnd_begin + nbnd_loc - 1 ) > nbnd ) nbnd_loc = nbnd - nbnd_begin + 1
       ELSE
          nproc = 1
          nbnd_loc = nbnd
          nbnd_begin = 1
       END IF

       ALLOCATE( work1( npwx ), work2( npwx ) ) 
       ALLOCATE( deff(nhm,nhm,nat) )
       !
       ! ... diagonal contribution - if the result from "calbec" are not 
       ! ... distributed, must be calculated on a single processor
       !
       evps = 0.D0
       IF ( nproc == 1 .AND. me_pool /= root_pool ) GO TO 100
       !
       DO ibnd_loc = 1, nbnd_loc
          ibnd = ibnd_loc + becp%ibnd_begin - 1 
          CALL compute_deff ( deff, et(ibnd,ik) )
          fac = wg(ibnd,ik)
          ijkb0 = 0
          DO np = 1, ntyp
             DO na = 1, nat
                IF ( ityp(na) == np ) THEN
                   DO ih = 1, nh(np)
                      ikb = ijkb0 + ih
                      evps = evps + fac * deff(ih,ih,na) * &
                                    ABS( becp%r(ikb,ibnd_loc) )**2
                      !
                      IF ( upf(np)%tvanp .OR. newpseudo(np) ) THEN
                         !
                         ! ... only in the US case there is a contribution 
                         ! ... for jh<>ih
                         ! ... we use here the symmetry in the interchange of 
                         ! ... ih and jh
                         !
                         DO jh = ( ih + 1 ), nh(np)
                            jkb = ijkb0 + jh
                            evps = evps + deff(ih,jh,na) * fac * 2.D0 * &
                                 becp%r(ikb,ibnd_loc) * becp%r(jkb,ibnd_loc)
                         END DO
                      END IF
                   END DO
                   ijkb0 = ijkb0 + nh(np)
                END IF
             END DO
          END DO
       END DO
       !
100    CONTINUE
       !
       ! ... non diagonal contribution - derivative of the bessel function
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ALLOCATE( dvkb( npwx, nkb ) )
       !
       CALL gen_us_dj( ik, dvkb )
       !
       DO icyc = 0, nproc -1
          !
          DO ibnd_loc = 1, nbnd_loc
             !  
             ibnd = ibnd_loc + becp%ibnd_begin - 1 
             CALL compute_deff ( deff, et(ibnd,ik) )
             work2(:) = (0.D0,0.D0)
             ijkb0 = 0
             DO np = 1, ntyp
                DO na = 1, nat
                   IF ( ityp(na) == np ) THEN
                      DO ih = 1, nh(np)
                         ikb = ijkb0 + ih
                         IF ( .NOT. ( upf(np)%tvanp .OR. newpseudo(np) ) ) THEN
                            ps = becp%r(ikb,ibnd_loc) * deff(ih,ih,na)
                         ELSE
                            !
                            ! ... in the US case there is a contribution 
                            ! ... also for jh<>ih
                            !
                            ps = (0.D0,0.D0)
                            DO jh = 1, nh(np)
                               jkb = ijkb0 + jh
                               ps = ps + becp%r(jkb,ibnd_loc) * deff(ih,jh,na)
                            END DO
                         END IF
                         CALL zaxpy( npw, ps, dvkb(1,ikb), 1, work2, 1 )
                      END DO
                      ijkb0 = ijkb0 + nh(np)
                   END IF
                END DO
             END DO
             !
             ! ... a factor 2 accounts for the other half of the G-vector sphere
             !
             DO ipol = 1, 3
                DO jpol = 1, ipol
                   DO i = 1, npw
                      work1(i) = evc(i,ibnd) * gk(ipol,i) * gk(jpol,i) * qm1(i)
                   END DO
                   sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - &
                        4.D0 * wg(ibnd,ik) * &
                        ddot( 2 * npw, work1, 1, work2, 1 )
                END DO
             END DO
          END DO
          IF ( nproc > 1 ) THEN
             CALL mp_circular_shift_left(becp%r, icyc, becp%comm)
             CALL mp_circular_shift_left(becp%ibnd_begin, icyc, becp%comm)
             CALL mp_circular_shift_left(nbnd_loc, icyc, becp%comm)
          END IF
       END DO
       !
       ! ... non diagonal contribution - derivative of the spherical harmonics
       ! ... (no contribution from l=0)
       !
       IF ( lmaxkb == 0 ) GO TO 10
       !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO ipol = 1, 3
          CALL gen_us_dy( ik, xyz(1,ipol), dvkb )
                
          DO icyc = 0, nproc -1
                
             DO ibnd_loc = 1, nbnd_loc
                ibnd = ibnd_loc + becp%ibnd_begin - 1 
                CALL compute_deff ( deff, et(ibnd,ik) )
                work2(:) = (0.D0,0.D0)
                ijkb0 = 0
                DO np = 1, ntyp
                   DO na = 1, nat
                      IF ( ityp(na) == np ) THEN
                         DO ih = 1, nh(np)
                            ikb = ijkb0 + ih
                            IF ( .NOT. ( upf(np)%tvanp .OR. newpseudo(np) ) ) THEN
                               ps = becp%r(ikb,ibnd_loc) * deff(ih,ih,na)
                            ELSE 
                               !
                               ! ... in the US case there is a contribution 
                               ! ... also for jh<>ih
                               !
                               ps = (0.D0,0.D0)
                               DO jh = 1, nh(np)
                                  jkb = ijkb0 + jh
                                  ps = ps + becp%r(jkb,ibnd_loc)*deff(ih,jh,na)
                               END DO
                            END IF
                            CALL zaxpy( npw, ps, dvkb(1,ikb), 1, work2, 1 )
                         END DO
                         ijkb0 = ijkb0 + nh(np)
                      END IF
                   END DO
                END DO
                !
                ! ... a factor 2 accounts for the other half of the G-vector sphere
                !
                DO jpol = 1, ipol
                   DO i = 1, npw
                      work1(i) = evc(i,ibnd) * gk(jpol,i)
                   END DO
                   sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - &
                        4.D0 * wg(ibnd,ik) * &
                        ddot( 2 * npw, work1, 1, work2, 1 )
                END DO
             END DO

             IF ( nproc > 1 ) THEN
                CALL mp_circular_shift_left(becp%r, icyc, becp%comm)
                CALL mp_circular_shift_left(becp%ibnd_begin, icyc, becp%comm)
                CALL mp_circular_shift_left(nbnd_loc, icyc, becp%comm)
             END IF

          ENDDO
       END DO

10     CONTINUE
       !
       DO l = 1, 3
          sigmanlc(l,l) = sigmanlc(l,l) - evps
       END DO
       !
       DEALLOCATE( dvkb )
       DEALLOCATE( deff, work2, work1 )
       !
       RETURN
       !
     END SUBROUTINE stres_us_gamma     
     !
     !
     !----------------------------------------------------------------------
     SUBROUTINE stres_us_k()
       !----------------------------------------------------------------------  
       !
       ! ... k-points version
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER                       :: na, np, ibnd, ipol, jpol, l, i, &
                                        ikb, jkb, ih, jh, ijkb0, is, js, ijs
       REAL(DP)                 :: fac, xyz (3, 3), evps, ddot
       COMPLEX(DP), ALLOCATABLE :: work1(:), work2(:), dvkb(:,:)
       COMPLEX(DP), ALLOCATABLE :: work2_nc(:,:)
       COMPLEX(DP), ALLOCATABLE :: deff_nc(:,:,:,:)
       REAL(DP), ALLOCATABLE :: deff(:,:,:)
       ! dvkb contains the derivatives of the kb potential
       COMPLEX(DP)              :: ps, ps_nc(2)
       ! xyz are the three unit vectors in the x,y,z directions
       DATA xyz / 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0 /
       !
       !
       if (noncolin) then
          ALLOCATE( work2_nc(npwx,npol) )
          ALLOCATE( deff_nc(nhm,nhm,nat,nspin) )
       else
          ALLOCATE( deff(nhm,nhm,nat) )
       endif
       !
       ALLOCATE( work1(npwx), work2(npwx) )
       !
       evps = 0.D0
       ! ... diagonal contribution
       !
       IF ( me_bgrp /= root_bgrp ) GO TO 100
       !
       ! ... the contribution is calculated only on one processor because
       ! ... partial results are later summed over all processors
       !
       DO ibnd = 1, nbnd
          fac = wg(ibnd,ik)
          IF (ABS(fac) < 1.d-9) CYCLE
          IF (noncolin) THEN
             CALL compute_deff_nc(deff_nc,et(ibnd,ik))
          ELSE
             CALL compute_deff(deff,et(ibnd,ik))
          ENDIF
          ijkb0 = 0
          DO np = 1, ntyp
             DO na = 1, nat
                IF ( ityp(na) == np ) THEN
                   DO ih = 1, nh(np)
                      ikb = ijkb0 + ih
                      IF (noncolin) THEN
                         ijs=0
                         DO is=1,npol
                            DO js=1,npol
                               ijs=ijs+1
                               evps=evps+fac*deff_nc(ih,ih,na,ijs)*   &
                                         CONJG(becp%nc(ikb,is,ibnd))* &
                                               becp%nc(ikb,js,ibnd)
                            END DO
                         END DO
                      ELSE
                         evps = evps+fac*deff(ih,ih,na)*ABS(becp%k(ikb,ibnd) )**2
                      END IF
                      !
                      IF ( upf(np)%tvanp .OR. newpseudo(np) ) THEN
                         !
                         ! ... only in the US case there is a contribution 
                         ! ... for jh<>ih
                         ! ... we use here the symmetry in the interchange of 
                         ! ... ih and jh
                         !
                         DO jh = ( ih + 1 ), nh(np)
                            jkb = ijkb0 + jh
                            IF (noncolin) THEN
                               ijs=0
                               DO is=1,npol
                                  DO js=1,npol
                                     ijs=ijs+1
                                     evps = evps+2.d0*fac&
                                            *DBLE(deff_nc(ih,jh,na,ijs)*      &
                                            (CONJG( becp%nc(ikb,is,ibnd) ) * &
                                                    becp%nc(jkb,js,ibnd))  )
                                  END DO
                               END DO
                            ELSE
                               evps = evps + deff(ih,jh,na) * fac * 2.D0 * &
                                     DBLE( CONJG( becp%k(ikb,ibnd) ) * &
                                                  becp%k(jkb,ibnd) )
                            END IF
                         END DO
                      END IF
                   END DO
                   ijkb0 = ijkb0 + nh(np)
                END IF
             END DO
          END DO
       END DO
       DO l = 1, 3
          sigmanlc(l,l) = sigmanlc(l,l) - evps
       END DO
       !
100    CONTINUE
       !
       ! ... non diagonal contribution - derivative of the bessel function
       !
       ALLOCATE( dvkb( npwx, nkb ) )
       !
       CALL gen_us_dj( ik, dvkb )
       !
       DO ibnd = 1, nbnd
          IF (noncolin) THEN
             work2_nc = (0.D0,0.D0)
             CALL compute_deff_nc(deff_nc,et(ibnd,ik))
          ELSE
             work2 = (0.D0,0.D0)
             CALL compute_deff(deff,et(ibnd,ik))
          ENDIF
          ijkb0 = 0
          DO np = 1, ntyp
             DO na = 1, nat
                IF ( ityp(na) == np ) THEN
                   DO ih = 1, nh(np)
                      ikb = ijkb0 + ih
                      IF ( .NOT. ( upf(np)%tvanp .OR. newpseudo(np) ) ) THEN
                         IF (noncolin) THEN
                            if (lspinorb) call errore('stres_us','wrong case',1)
                            ijs=0
                            ps_nc=(0.D0, 0.D0)
                            DO is=1,npol
                               DO js=1,npol
                                  ijs=ijs+1
                                  ps_nc(is)=ps_nc(is)+becp%nc(ikb,js,ibnd)* &
                                         deff_nc(ih,ih,na,ijs)
                               END DO
                            END DO
                         ELSE
                            ps = becp%k(ikb, ibnd) * deeq(ih,ih,na,current_spin)
                         ENDIF
                      ELSE
                         !
                         ! ... in the US case there is a contribution 
                         ! ... also for jh<>ih
                         !
                         ps = (0.D0,0.D0)
                         ps_nc = (0.D0,0.D0)
                         DO jh = 1, nh(np)
                            jkb = ijkb0 + jh
                            IF (noncolin) THEN
                               ijs=0
                               DO is=1,npol
                                  DO js=1,npol
                                     ijs=ijs+1
                                     ps_nc(is)=ps_nc(is)+becp%nc(jkb,js,ibnd)* &
                                           deff_nc(ih,jh,na,ijs)
                                  END DO
                               END DO
                            ELSE
                               ps = ps + becp%k(jkb,ibnd) * deff(ih,jh,na)
                            END IF
                         END DO
                      END IF
                      IF (noncolin) THEN
                         DO is=1,npol
                            CALL zaxpy(npw,ps_nc(is),dvkb(1,ikb),1,&
                                                      work2_nc(1,is),1)
                         END DO
                      ELSE
                         CALL zaxpy( npw, ps, dvkb(1,ikb), 1, work2, 1 )
                      END IF
                   END DO
                   ijkb0 = ijkb0 + nh(np)
                END IF
             END DO
          END DO
          DO ipol = 1, 3
             DO jpol = 1, ipol
                IF (noncolin) THEN
                   DO i = 1, npw
                      work1(i) = evc(i     ,ibnd)*gk(ipol,i)* &
                                                  gk(jpol,i)*qm1(i)
                      work2(i) = evc(i+npwx,ibnd)*gk(ipol,i)* &
                                                  gk(jpol,i)*qm1(i)
                   END DO
                   sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - &
                                   2.D0 * wg(ibnd,ik) * &
                                 ( ddot(2*npw,work1,1,work2_nc(1,1), 1) + &
                                   ddot(2*npw,work2,1,work2_nc(1,2), 1) )
                ELSE
                   DO i = 1, npw
                      work1(i) = evc(i,ibnd)*gk(ipol,i)*gk(jpol,i)*qm1(i)
                   END DO
                   sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - &
                                      2.D0 * wg(ibnd,ik) * &
                                      ddot( 2 * npw, work1, 1, work2, 1 )
                END IF
             END DO
          END DO
       END DO
       !
       ! ... non diagonal contribution - derivative of the spherical harmonics
       ! ... (no contribution from l=0)
       !
       IF ( lmaxkb == 0 ) GO TO 10
       !
       DO ipol = 1, 3
          CALL gen_us_dy( ik, xyz(1,ipol), dvkb )
          DO ibnd = 1, nbnd
             IF (noncolin) THEN
                work2_nc = (0.D0,0.D0)
                CALL compute_deff_nc(deff_nc,et(ibnd,ik))
             ELSE
                work2 = (0.D0,0.D0)
                CALL compute_deff(deff,et(ibnd,ik))
             ENDIF

             ijkb0 = 0
             DO np = 1, ntyp
                DO na = 1, nat
                   IF ( ityp(na) == np ) THEN
                      DO ih = 1, nh(np)
                         ikb = ijkb0 + ih
                         IF ( .NOT. ( upf(np)%tvanp .OR. newpseudo(np) ) ) THEN
                            IF (noncolin) THEN
                               ijs=0
                               ps_nc = (0.D0,0.D0)
                               DO is=1,npol
                                  DO js=1,npol
                                     ijs=ijs+1
                                     ps_nc(is)=ps_nc(is)+becp%nc(ikb,js,ibnd)* &
                                         deff_nc(ih,ih,na,ijs)
                                  END DO
                               END DO
                            ELSE
                               ps = becp%k(ikb,ibnd) * deeq(ih,ih,na,current_spin)
                            END IF
                         ELSE
                            !
                            ! ... in the US case there is a contribution 
                            ! ... also for jh<>ih
                            !
                            ps = (0.D0,0.D0)
                            ps_nc = (0.D0,0.D0)
                            DO jh = 1, nh(np)
                               jkb = ijkb0 + jh
                               IF (noncolin) THEN
                                  ijs=0
                                  DO is=1,npol
                                     DO js=1,npol
                                        ijs=ijs+1
                                        ps_nc(is)=ps_nc(is)+ &
                                               becp%nc(jkb,js,ibnd)* &
                                               deff_nc(ih,jh,na,ijs)
                                     END DO
                                  END DO
                               ELSE
                                  ps = ps + becp%k(jkb,ibnd) * deff(ih,jh,na)
                               END IF
                            END DO
                         END IF
                         IF (noncolin) THEN
                            DO is=1,npol
                               CALL zaxpy(npw,ps_nc(is),dvkb(1,ikb),1, &
                                          work2_nc(1,is),1)
                            END DO
                         ELSE
                            CALL zaxpy( npw, ps, dvkb(1,ikb), 1, work2, 1 )
                         END IF
                      END DO
                      ijkb0 = ijkb0 + nh(np)
                   END IF
                END DO
             END DO
             DO jpol = 1, ipol
                IF (noncolin) THEN
                   DO i = 1, npw
                      work1(i) = evc(i     ,ibnd) * gk(jpol,i)
                      work2(i) = evc(i+npwx,ibnd) * gk(jpol,i)
                   END DO
                   sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - &
                              2.D0 * wg(ibnd,ik) * & 
                            ( ddot( 2 * npw, work1, 1, work2_nc(1,1), 1 ) + &
                              ddot( 2 * npw, work2, 1, work2_nc(1,2), 1 ) )
                ELSE
                   DO i = 1, npw
                      work1(i) = evc(i,ibnd) * gk(jpol,i)
                   END DO
                   sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - &
                                      2.D0 * wg(ibnd,ik) * & 
                                      ddot( 2 * npw, work1, 1, work2, 1 )
                END IF
             END DO
          END DO
       END DO
       !
10     CONTINUE
       !
       IF (noncolin) THEN
           DEALLOCATE( work2_nc )
           DEALLOCATE( deff_nc )
       ELSE
           DEALLOCATE( work2 )
           DEALLOCATE( deff )
       ENDIF
       DEALLOCATE( dvkb )
       DEALLOCATE( work1 )
       !
       RETURN
       !
     END SUBROUTINE stres_us_k
     !
END SUBROUTINE stres_us

#ifdef USE_CUDA

SUBROUTINE myDdotAB( n, A, B, res )
  !
  use kinds, ONLY : DP
  use cudafor
  use cublas
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n
  !
!!!!!pgi$ ignore_tkr A, B
  REAL(DP), DEVICE, INTENT(IN) :: A(n), B(n)
  !
  REAL(DP), INTENT(OUT) :: res
  !
  res = cublasDdot( n, A, 1, B, 1 )
  !
  return
  !
END SUBROUTINE myDdotAB


!----------------------------------------------------------------------------
SUBROUTINE stres_us_gpu( ik, gk_d, sigmanlc )
  !----------------------------------------------------------------------------
  !
  ! nonlocal (separable pseudopotential) contribution to the stress
  ! NOTICE: sum of partial results over procs is performed in calling routine
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, ityp_d
  USE constants,            ONLY : eps8
  USE klist,                ONLY : nks, xk, ngk, igk_k
  USE lsda_mod,             ONLY : current_spin, lsda, isk
  USE wvfct,                ONLY : npwx, nbnd, wg, wg_d, et, et_d
  USE control_flags,        ONLY : gamma_only
  USE uspp_param,           ONLY : upf, lmaxkb, nh, nh_d, newpseudo, nhm
  USE uspp,                 ONLY : nkb, vkb, qq, deeq, qq_d, deeq_d, deeq_nc, qq_so
  USE wavefunctions_module, ONLY : evc
  USE spin_orb,             ONLY : lspinorb
  USE lsda_mod,             ONLY : nspin
  USE noncollin_module,     ONLY : noncolin, npol
  USE mp_pools,             ONLY : me_pool, root_pool
  USE mp_bands,             ONLY : intra_bgrp_comm, me_bgrp, root_bgrp
  USE becmod,               ONLY : allocate_bec_type, deallocate_bec_type, &
                                   bec_type, becp, calbec
  USE mp,                   ONLY : mp_sum, mp_get_comm_null, mp_circular_shift_left 
  !
  USE cudafor
  USE cublas
  USE wavefunctions_module, ONLY : evc_d
  USE uspp,                 ONLY : vkb_d, qq_d, deeq_d
  USE klist,                ONLY : nks, xk, ngk, igk_k_d
  USE ep_debug,             ONLY : compare

  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)    :: ik
  REAL(DP), DEVICE, INTENT(IN)   :: gk_d(3,npwx)
  REAL(DP), INTENT(INOUT):: sigmanlc(3,3)
  !
  REAL(DP), ALLOCATABLE  :: qm1(:)
  REAL(DP)               :: q
  INTEGER                :: npw, i
  !
  REAL(DP), ALLOCATABLE   :: gk(:,:)
  REAL(DP), DEVICE, ALLOCATABLE  :: qm1_d(:)
  !
  IF ( nkb == 0 ) RETURN
  !
  IF ( lsda ) current_spin = isk(ik)
  npw = ngk(ik)

!#ifdef USE_CUDA
!  vkb = vkb_d
!#endif

  IF ( nks > 1 ) THEN
    !CALL init_us_2    ( npw, igk_k(1,ik), xk(1,ik), vkb )
    CALL init_us_2_gpu( npw, igk_k_d(1,ik), xk(1,ik), vkb_d )
    !CALL compare( vkb, vkb_d, "vkb" )
  ENDIF
  !
  CALL allocate_bec_type ( nkb, nbnd, becp, intra_bgrp_comm ) 
  !CALL calbec( npw, vkb, evc, becp )
  CALL calbec( npw, vkb_d, evc_d, becp )
  !CALL compare( becp%k, becp%k_d, "becp_k")
  !
  !ALLOCATE( gk, source=gk_d )
  ALLOCATE( qm1_d( npwx ) )
!$cuf kernel do(1) <<<*,*>>>
  DO i = 1, npw 
     q = SQRT( gk_d(1,i)*gk_d(1,i) + gk_d(2,i)*gk_d(2,i) + gk_d(3,i)*gk_d(3,i) )
     IF ( q > eps8 ) THEN
        qm1_d(i) = 1.D0 / q   
     ELSE
        qm1_d(i) = 0.D0
     END IF
  END DO

!  ALLOCATE( qm1( npwx ) )
!  DO i = 1, npw
!     q = SQRT( gk(1,i)**2 + gk(2,i)**2 + gk(3,i)**2 )
!     IF ( q > eps8 ) THEN
!        qm1(i) = 1.D0 / q
!     ELSE
!        qm1(i) = 0.D0
!     END IF
!  END DO

  !CALL compare( qm1(1:npw), qm1_d(1:npw),"qm1" )
  !
  IF ( gamma_only ) THEN
     !
     CALL stres_us_gamma_gpu()
     !
  ELSE
     !
     CALL stres_us_k_gpu(gk_d)
     !
  END IF
  !
  DEALLOCATE( qm1_d )
  !CALL deallocate_bec_type ( becp ) 
  !
  RETURN
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     SUBROUTINE stres_us_gamma_gpu()
       !-----------------------------------------------------------------------
       ! 
       ! ... gamma version
       !
       IMPLICIT NONE
       !
       ! ... local variables
       !
       INTEGER                       :: na, np, ibnd, ipol, jpol, l, i, &
                                        ikb, jkb, ih, jh, ijkb0, ibnd_loc, &
                                        nproc, nbnd_loc, nbnd_begin, icyc
       INTEGER, EXTERNAL :: ldim_block, lind_block, gind_block
       REAL(DP)                 :: fac, xyz(3,3), evps, ddot
       REAL(DP), ALLOCATABLE    :: deff(:,:,:)
       COMPLEX(DP), ALLOCATABLE :: work1(:), work2(:), dvkb(:,:)
       ! dvkb contains the derivatives of the kb potential
       COMPLEX(DP)              :: ps
       ! xyz are the three unit vectors in the x,y,z directions
       DATA xyz / 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0 /
       !
       !
       IF( becp%comm /= mp_get_comm_null() ) THEN
          nproc   = becp%nproc
          nbnd_loc   = becp%nbnd_loc
          nbnd_begin = becp%ibnd_begin
          IF( ( nbnd_begin + nbnd_loc - 1 ) > nbnd ) nbnd_loc = nbnd - nbnd_begin + 1
       ELSE
          nproc = 1
          nbnd_loc = nbnd
          nbnd_begin = 1
       END IF

       ALLOCATE( work1( npwx ), work2( npwx ) ) 
       ALLOCATE( deff(nhm,nhm,nat) )
       !
       ! ... diagonal contribution - if the result from "calbec" are not 
       ! ... distributed, must be calculated on a single processor
       !
       evps = 0.D0
       IF ( nproc == 1 .AND. me_pool /= root_pool ) GO TO 100
       !
       DO ibnd_loc = 1, nbnd_loc
          ibnd = ibnd_loc + becp%ibnd_begin - 1 
          CALL compute_deff ( deff, et(ibnd,ik) )
          fac = wg(ibnd,ik)
          ijkb0 = 0
          DO np = 1, ntyp
             DO na = 1, nat
                IF ( ityp(na) == np ) THEN
                   DO ih = 1, nh(np)
                      ikb = ijkb0 + ih
                      evps = evps + fac * deff(ih,ih,na) * &
                                    ABS( becp%r(ikb,ibnd_loc) )**2
                      !
                      IF ( upf(np)%tvanp .OR. newpseudo(np) ) THEN
                         !
                         ! ... only in the US case there is a contribution 
                         ! ... for jh<>ih
                         ! ... we use here the symmetry in the interchange of 
                         ! ... ih and jh
                         !
                         DO jh = ( ih + 1 ), nh(np)
                            jkb = ijkb0 + jh
                            evps = evps + deff(ih,jh,na) * fac * 2.D0 * &
                                 becp%r(ikb,ibnd_loc) * becp%r(jkb,ibnd_loc)
                         END DO
                      END IF
                   END DO
                   ijkb0 = ijkb0 + nh(np)
                END IF
             END DO
          END DO
       END DO
       !
100    CONTINUE
       !
       ! ... non diagonal contribution - derivative of the bessel function
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ALLOCATE( dvkb( npwx, nkb ) )
       !
       CALL gen_us_dj( ik, dvkb )
       !
       DO icyc = 0, nproc -1
          !
          DO ibnd_loc = 1, nbnd_loc
             !  
             ibnd = ibnd_loc + becp%ibnd_begin - 1 
             CALL compute_deff ( deff, et(ibnd,ik) )
             work2(:) = (0.D0,0.D0)
             ijkb0 = 0
             DO np = 1, ntyp
                DO na = 1, nat
                   IF ( ityp(na) == np ) THEN
                      DO ih = 1, nh(np)
                         ikb = ijkb0 + ih
                         IF ( .NOT. ( upf(np)%tvanp .OR. newpseudo(np) ) ) THEN
                            ps = becp%r(ikb,ibnd_loc) * deff(ih,ih,na)
                         ELSE
                            !
                            ! ... in the US case there is a contribution 
                            ! ... also for jh<>ih
                            !
                            ps = (0.D0,0.D0)
                            DO jh = 1, nh(np)
                               jkb = ijkb0 + jh
                               ps = ps + becp%r(jkb,ibnd_loc) * deff(ih,jh,na)
                            END DO
                         END IF
                         CALL zaxpy( npw, ps, dvkb(1,ikb), 1, work2, 1 )
                      END DO
                      ijkb0 = ijkb0 + nh(np)
                   END IF
                END DO
             END DO
             !
             ! ... a factor 2 accounts for the other half of the G-vector sphere
             !
             DO ipol = 1, 3
                DO jpol = 1, ipol
                   DO i = 1, npw
                      work1(i) = evc(i,ibnd) * gk(ipol,i) * gk(jpol,i) * qm1(i)
                   END DO
                   sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - &
                        4.D0 * wg(ibnd,ik) * &
                        ddot( 2 * npw, work1, 1, work2, 1 )
                END DO
             END DO
          END DO
          IF ( nproc > 1 ) THEN
             CALL mp_circular_shift_left(becp%r, icyc, becp%comm)
             CALL mp_circular_shift_left(becp%ibnd_begin, icyc, becp%comm)
             CALL mp_circular_shift_left(nbnd_loc, icyc, becp%comm)
          END IF
       END DO
       !
       ! ... non diagonal contribution - derivative of the spherical harmonics
       ! ... (no contribution from l=0)
       !
       IF ( lmaxkb == 0 ) GO TO 10
       !
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO ipol = 1, 3
          CALL gen_us_dy( ik, xyz(1,ipol), dvkb )
                
          DO icyc = 0, nproc -1
                
             DO ibnd_loc = 1, nbnd_loc
                ibnd = ibnd_loc + becp%ibnd_begin - 1 
                CALL compute_deff ( deff, et(ibnd,ik) )
                work2(:) = (0.D0,0.D0)
                ijkb0 = 0
                DO np = 1, ntyp
                   DO na = 1, nat
                      IF ( ityp(na) == np ) THEN
                         DO ih = 1, nh(np)
                            ikb = ijkb0 + ih
                            IF ( .NOT. ( upf(np)%tvanp .OR. newpseudo(np) ) ) THEN
                               ps = becp%r(ikb,ibnd_loc) * deff(ih,ih,na)
                            ELSE 
                               !
                               ! ... in the US case there is a contribution 
                               ! ... also for jh<>ih
                               !
                               ps = (0.D0,0.D0)
                               DO jh = 1, nh(np)
                                  jkb = ijkb0 + jh
                                  ps = ps + becp%r(jkb,ibnd_loc)*deff(ih,jh,na)
                               END DO
                            END IF
                            CALL zaxpy( npw, ps, dvkb(1,ikb), 1, work2, 1 )
                         END DO
                         ijkb0 = ijkb0 + nh(np)
                      END IF
                   END DO
                END DO
                !
                ! ... a factor 2 accounts for the other half of the G-vector sphere
                !
                DO jpol = 1, ipol
                   DO i = 1, npw
                      work1(i) = evc(i,ibnd) * gk(jpol,i)
                   END DO
                   sigmanlc(ipol,jpol) = sigmanlc(ipol,jpol) - &
                        4.D0 * wg(ibnd,ik) * &
                        ddot( 2 * npw, work1, 1, work2, 1 )
                END DO
             END DO

             IF ( nproc > 1 ) THEN
                CALL mp_circular_shift_left(becp%r, icyc, becp%comm)
                CALL mp_circular_shift_left(becp%ibnd_begin, icyc, becp%comm)
                CALL mp_circular_shift_left(nbnd_loc, icyc, becp%comm)
             END IF

          ENDDO
       END DO

10     CONTINUE
       !
       DO l = 1, 3
          sigmanlc(l,l) = sigmanlc(l,l) - evps
       END DO
       !
       DEALLOCATE( dvkb )
       DEALLOCATE( deff, work2, work1 )
       !
       RETURN
       !
     END SUBROUTINE stres_us_gamma_gpu     
     !
     !
     !----------------------------------------------------------------------
     SUBROUTINE stres_us_k_gpu(gk_d_)
       !----------------------------------------------------------------------  
       !
       ! ... k-points version
       !
       IMPLICIT NONE
       !
       REAL(DP), DEVICE, INTENT(IN)   :: gk_d_(3,npwx) ! workaround for PGI bug
       ! ... local variables
       !
       INTEGER                       :: na, np, ibnd, ipol, jpol, l, i, &
                                        ikb, jkb, ih, jh, ijkb0, is, js, ijs
       REAL(DP)                 :: fac, xyz (3, 3), evps, ddot
       COMPLEX(DP), ALLOCATABLE :: work1(:), work2(:), dvkb(:,:)
       COMPLEX(DP), ALLOCATABLE :: work2_nc(:,:)
       COMPLEX(DP), ALLOCATABLE :: deff_nc(:,:,:,:)
       REAL(DP), ALLOCATABLE :: deff(:,:,:)
       ! dvkb contains the derivatives of the kb potential
       COMPLEX(DP)              :: ps, ps_nc(2), cv, cv2
       ! xyz are the three unit vectors in the x,y,z directions
       DATA xyz / 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0 /
       !
       COMPLEX(DP), DEVICE, ALLOCATABLE :: dvkb_d(:,:,:), ps_d(:,:)
       REAL(DP), DEVICE, ALLOCATABLE :: deff_d(:,:,:)
       REAL(DP) :: workdot, dot11, dot21, dot22, dot31, dot32, dot33
       REAL(DP) :: gk1, gk2, gk3, wgi, temp
       COMPLEX(DP) :: work_sum, wsum1,wsum2,wsum3
       COMPLEX(DP) :: evci, qm1i
       COMPLEX(DP), DEVICE, POINTER :: becp_k(:,:)
       integer, allocatable :: shift(:)
       integer, allocatable, device :: shift_d(:)
       integer :: nhnp, nhmax
       !
       becp_k => becp%k_d


       !
       if (noncolin) then
          ALLOCATE( work2_nc(npwx,npol) )
          ALLOCATE( deff_nc(nhm,nhm,nat,nspin) )
       else
          ALLOCATE( deff(nhm,nhm,nat) )
          ALLOCATE( deff_d(nhm,nhm,nat) )
       endif
       !
       ALLOCATE( work1(npwx), work2(npwx) )
       !
       evps = 0.D0
       ! ... diagonal contribution
       !
       IF ( me_bgrp /= root_bgrp ) GO TO 100
       !
       ! ... the contribution is calculated only on one processor because
       ! ... partial results are later summed over all processors
       !
       DO ibnd = 1, nbnd
          fac = wg(ibnd,ik)
          IF (ABS(fac) < 1.d-9) CYCLE
          IF (noncolin) THEN
             CALL compute_deff_nc(deff_nc,et(ibnd,ik))
          ELSE
             CALL compute_deff(deff,et(ibnd,ik))
          ENDIF
          ijkb0 = 0
          DO np = 1, ntyp
             DO na = 1, nat
                IF ( ityp(na) == np ) THEN
                   DO ih = 1, nh(np)
                      ikb = ijkb0 + ih
                      IF (noncolin) THEN
                         ijs=0
                         DO is=1,npol
                            DO js=1,npol
                               ijs=ijs+1
                               evps=evps+fac*deff_nc(ih,ih,na,ijs)*   &
                                         CONJG(becp%nc(ikb,is,ibnd))* &
                                               becp%nc(ikb,js,ibnd)
                            END DO
                         END DO
                      ELSE
                         evps = evps+fac*deff(ih,ih,na)*ABS(becp%k(ikb,ibnd) )**2
                      END IF
                      !
                      IF ( upf(np)%tvanp .OR. newpseudo(np) ) THEN
                         !
                         ! ... only in the US case there is a contribution 
                         ! ... for jh<>ih
                         ! ... we use here the symmetry in the interchange of 
                         ! ... ih and jh
                         !
                         DO jh = ( ih + 1 ), nh(np)
                            jkb = ijkb0 + jh
                            IF (noncolin) THEN
                               ijs=0
                               DO is=1,npol
                                  DO js=1,npol
                                     ijs=ijs+1
                                     evps = evps+2.d0*fac&
                                            *DBLE(deff_nc(ih,jh,na,ijs)*      &
                                            (CONJG( becp%nc(ikb,is,ibnd) ) * &
                                                    becp%nc(jkb,js,ibnd))  )
                                  END DO
                               END DO
                            ELSE
                               evps = evps + deff(ih,jh,na) * fac * 2.D0 * &
                                     DBLE( CONJG( becp%k(ikb,ibnd) ) * &
                                                  becp%k(jkb,ibnd) )
                            END IF
                         END DO
                      END IF
                   END DO
                   ijkb0 = ijkb0 + nh(np)
                END IF
             END DO
          END DO
       END DO
       DO l = 1, 3
          sigmanlc(l,l) = sigmanlc(l,l) - evps
       END DO
       !
100    CONTINUE
       !
       ! ... non diagonal contribution - derivative of the bessel function
       !
       ALLOCATE( dvkb_d( npwx, nkb, 4 ) )
       !
       CALL gen_us_dj_gpu( ik, dvkb_d(:,:,4) )
       !
       IF ( lmaxkb .ne. 0 ) THEN
         DO ipol = 1, 3
            CALL gen_us_dy_gpu( ik, xyz(1,ipol), dvkb_d(:,:,ipol) )
         END DO
       ENDIF
   
       ! Precompute ijkb0 shifts (can be avoided if atoms are sorted by type)
       ALLOCATE(shift(nat))
       ijkb0 = 0
       DO np = 1, ntyp
         DO na = 1, nat
           IF ( ityp(na) == np ) THEN
             shift(na) = ijkb0
             ijkb0 = ijkb0 + nh(np)
           ENDIF
         END DO
       END DO
       ALLOCATE(shift_d, source = shift)

       nhmax = 0
       do np = 1, ntyp
         nhmax = max(nhmax, nh(np))
       end do

       ALLOCATE(ps_d(nhmax*nat, ntyp))

       DO ibnd = 1, nbnd
         !IF (noncolin) THEN
           !CALL compute_deff_nc(deff_nc,et(ibnd,ik))
         !ELSE
           CALL compute_deff_gpu(deff_d,et(ibnd,ik))
         !ENDIF

         ! Compute PS
         ps_d = 0.d0
         !$cuf kernel do (1) <<<*,*>>>
         DO na = 1, nat
           np = ityp_d(na)
           ijkb0 = shift_d(na)
           nhnp = nh_d(np)

           DO ih = 1, nhnp
             ikb = ijkb0 + ih
             ! TODO: Alternative case not currently supported. Needs to be updated.
             !IF ( .NOT. ( upf(np)%tvanp .OR. newpseudo(np) ) ) THEN
             !  IF (noncolin) THEN
             !    if (lspinorb) call errore('stres_us','wrong case',1)
             !    ijs=0
             !    ps_nc=(0.D0, 0.D0)
             !    DO is=1,npol
             !      DO js=1,npol
             !         ijs=ijs+1
             !         ps_nc_d(is,ikb,np)=ps_nc_d(is,ikb,np)+becp%nc(ikb,js,ibnd)* &
             !         deff_nc(ih,ih,na,ijs)
             !      END DO
             !    END DO
             !  ELSE
             !    ps_d(ikb, np) = becp%k(ikb, ibnd) * deeq(ih,ih,na,current_spin)
             !  ENDIF
             !ELSE
               DO jh = 1, nhnp
                 jkb = ijkb0 + jh
                 ps_d(ikb, np) = ps_d(ikb, np) + becp_k(jkb,ibnd) * deff_d(ih,jh,na)
               END DO
             !ENDIF
           END DO

         END DO

         dot11 = 0.d0 
         dot21 = 0.d0 
         dot22 = 0.d0 
         dot31 = 0.d0 
         dot32 = 0.d0 
         dot33 = 0.d0 

         !$cuf kernel do(2) <<<*,*>>>
         DO na = 1, nat
           DO i = 1, npw
             work_sum = (0.D0,0.D0)

             np = ityp_d(na)
             ijkb0 = shift_d(na)
             nhnp = nh_d(np)

             DO ih = 1, nhnp
               ikb = ijkb0 + ih
               work_sum = work_sum + ps_d(ikb, np)*dvkb_d(i,ikb,4)
             END DO

             evci = evc_d(i, ibnd)
             qm1i = qm1_d(i)
             gk1 = gk_d_(1, i)
             gk2 = gk_d_(2, i)
             gk3 = gk_d_(3, i)


             !TODO: add updated logic for noncolinear case here
             cv = evci * gk1 * gk1 * qm1i
             dot11 = dot11 + work_sum%re * cv%re + work_sum%im * cv%im

             cv = evci * gk2 * gk1 * qm1i
             dot21 = dot21 + work_sum%re * cv%re + work_sum%im * cv%im

             cv = evci * gk3 * gk1 * qm1i
             dot31 = dot31 + work_sum%re * cv%re + work_sum%im * cv%im

             cv = evci * gk2 * gk2 * qm1i
             dot22 = dot22 + work_sum%re * cv%re + work_sum%im * cv%im

             dot32 = dot32 + work_sum%re * cv%re + work_sum%im * cv%im

             cv = evci * gk3 * gk3 * qm1i
             dot33 = dot33 + work_sum%re * cv%re + work_sum%im * cv%im
           END DO
         END DO

         !TODO: add updated logic for noncolinear case here
         sigmanlc(1,1) = sigmanlc(1,1) - 2.D0 * wg(ibnd,ik) * dot11
         sigmanlc(2,1) = sigmanlc(2,1) - 2.D0 * wg(ibnd,ik) * dot21
         sigmanlc(2,2) = sigmanlc(2,2) - 2.D0 * wg(ibnd,ik) * dot22
         sigmanlc(3,1) = sigmanlc(3,1) - 2.D0 * wg(ibnd,ik) * dot31
         sigmanlc(3,2) = sigmanlc(3,2) - 2.D0 * wg(ibnd,ik) * dot32
         sigmanlc(3,3) = sigmanlc(3,3) - 2.D0 * wg(ibnd,ik) * dot33

         !
         ! ... non diagonal contribution - derivative of the spherical harmonics
         ! ... (no contribution from l=0)
         !
         IF ( lmaxkb == 0 ) CYCLE
         !
       
         dot11 = 0.d0 
         dot21 = 0.d0 
         dot22 = 0.d0 
         dot31 = 0.d0 
         dot32 = 0.d0 
         dot33 = 0.d0 

         !$cuf kernel do(2) <<<*,*>>>
         DO na = 1, nat
           DO i = 1, npw
             np = ityp_d(na)
             ijkb0 = shift_d(na)
             nhnp = nh_d(np)

             wsum1 = (0.D0,0.D0)
             wsum2 = (0.D0,0.D0)
             wsum3 = (0.D0,0.D0)

             DO ih = 1, nhnp
               ikb = ijkb0 + ih
               gk1 = ps_d(ikb, np)

               wsum1 = wsum1 + gk1*dvkb_d(i,ikb,1)
               wsum2 = wsum2 + gk1*dvkb_d(i,ikb,2)
               wsum3 = wsum3 + gk1*dvkb_d(i,ikb,3)
             END DO
             
             evci = evc_d(i, ibnd)
             gk1 = gk_d_(1, i)
             gk2 = gk_d_(2, i)
             gk3 = gk_d_(3, i)

             !TODO: add updated logic for noncolinear case here
             cv = evci * gk1
             dot11 = dot11 + wsum1%re * cv%re + wsum1%im * cv%im
             dot21 = dot21 + wsum2%re * cv%re + wsum2%im * cv%im
             dot31 = dot31 + wsum3%re * cv%re + wsum3%im * cv%im

             cv = evci * gk2
             dot22 = dot22 + wsum2%re * cv%re + wsum2%im * cv%im
             dot32 = dot32 + wsum3%re * cv%re + wsum3%im * cv%im

             cv = evci * gk3
             dot33 = dot33 + wsum3%re * cv%re + wsum3%im * cv%im

           END DO
         END DO

         !TODO: add updated logic for noncolinear case here
         sigmanlc(1,1) = sigmanlc(1,1) - 2.D0 * wg(ibnd,ik) * dot11
         sigmanlc(2,1) = sigmanlc(2,1) - 2.D0 * wg(ibnd,ik) * dot21
         sigmanlc(2,2) = sigmanlc(2,2) - 2.D0 * wg(ibnd,ik) * dot22
         sigmanlc(3,1) = sigmanlc(3,1) - 2.D0 * wg(ibnd,ik) * dot31
         sigmanlc(3,2) = sigmanlc(3,2) - 2.D0 * wg(ibnd,ik) * dot32
         sigmanlc(3,3) = sigmanlc(3,3) - 2.D0 * wg(ibnd,ik) * dot33

       END DO
       !
10     CONTINUE
       !
       IF (noncolin) THEN
           DEALLOCATE( work2_nc )
           DEALLOCATE( deff_nc )
       ELSE
           DEALLOCATE( work2 )
           DEALLOCATE( deff )
           DEALLOCATE( deff_d )
       ENDIF
       !DEALLOCATE( dvkb )
       DEALLOCATE( work1 )
       DEALLOCATE( dvkb_d )
       !
       DEALLOCATE(shift)
       DEALLOCATE(shift_d)
       DEALLOCATE(ps_d)
       RETURN
       !
     END SUBROUTINE stres_us_k_gpu
     !
END SUBROUTINE stres_us_gpu

#endif
