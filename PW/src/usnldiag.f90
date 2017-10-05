!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE usnldiag (npw, h_diag, s_diag)
  !-----------------------------------------------------------------------
  !
  !    add nonlocal pseudopotential term to diagonal part of Hamiltonian
  !    compute the diagonal part of the S matrix
  !
  USE kinds, ONLY: DP
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp
  USE wvfct, ONLY: npwx
  USE lsda_mod, ONLY: current_spin
  USE uspp,  ONLY: deeq, vkb, qq, qq_so, deeq_nc, indv_ijkb0
  USE uspp_param, ONLY: upf, nh, newpseudo
  USE spin_orb, ONLY: lspinorb
  USE noncollin_module, ONLY: noncolin, npol
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(in) :: npw
  ! number of plane waves
  !
  REAL(dp), INTENT(inout) :: h_diag (npwx,npol)
  ! the diagonal part of the hamiltonian
  REAL(dp), INTENT(out)   :: s_diag (npwx,npol)
  ! the diagonal part of the S matrix
  !
  INTEGER :: ikb, jkb, ih, jh, na, nt, ig, ipol
  COMPLEX(DP) :: ps1(2), ps2(2), ar
  !
  ! initialise s_diag
  !
  s_diag = 1.d0
  !
  !    multiply on projectors
  !
  DO nt = 1, ntyp
     DO na = 1, nat
        IF (ityp (na) == nt) THEN
           DO ih = 1, nh(nt)
              ikb = indv_ijkb0(na) + ih
              IF (lspinorb) THEN
                 ps1(1) = deeq_nc (ih, ih, na, 1)
                 ps1(2) = deeq_nc (ih, ih, na, 4)
                 ps2(1) = qq_so(ih, ih, 1, nt)
                 ps2(2) = qq_so(ih, ih, 4, nt)
              ELSEIF (noncolin) THEN
                 ps1(1) = deeq_nc (ih, ih, na, 1)
                 ps1(2) = deeq_nc (ih, ih, na, 4)
                 ps2(1) = qq (ih, ih, nt)
                 ps2(2) = qq (ih, ih, nt)
              ELSE
                 ps1(1) = deeq (ih, ih, na, current_spin)
                 ps2(1) = qq (ih, ih, nt)
              ENDIF
              DO ipol =1, npol
                 DO ig = 1, npw
                    ar = vkb (ig, ikb)*conjg(vkb (ig, ikb))
                    h_diag (ig,ipol) = h_diag (ig,ipol) + ps1(ipol) * ar
                    s_diag (ig,ipol) = s_diag (ig,ipol) + ps2(ipol) * ar
                 ENDDO
              ENDDO
              IF ( upf(nt)%tvanp .or.newpseudo (nt) ) THEN
                 DO jh = 1, nh (nt)
                    IF (jh/=ih) THEN
                       jkb = indv_ijkb0(na) + jh
                       IF (lspinorb) THEN
                          ps1(1) = deeq_nc (ih, jh, na, 1)
                          ps1(2) = deeq_nc (ih, jh, na, 4)
                          ps2(1) = qq_so(ih, jh, 1, nt)
                          ps2(2) = qq_so(ih, jh, 4, nt)
                       ELSEIF (noncolin) THEN
                          ps1(1) = deeq_nc (ih, jh, na, 1)
                          ps1(2) = deeq_nc (ih, jh, na, 4)
                          ps2(1) = qq (ih, jh, nt)
                          ps2(2) = qq (ih, jh, nt)
                       ELSE
                          ps1(1) = deeq (ih, jh, na, current_spin)
                          ps2(1) = qq (ih, jh, nt)
                       ENDIF
                       DO ipol = 1, npol
                          DO ig = 1, npw
                             ar = vkb (ig, ikb) *conjg( vkb (ig, jkb))
                             h_diag (ig,ipol) = h_diag (ig,ipol) + &
                                  ps1(ipol) * ar
                             s_diag (ig,ipol) = s_diag (ig,ipol) + &
                                  ps2(ipol) * ar
                          ENDDO
                       ENDDO
                    ENDIF
                 ENDDO
              ENDIF
           ENDDO
        ENDIF
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE usnldiag

#ifdef USE_CUDA
!-----------------------------------------------------------------------
SUBROUTINE usnldiag_gpu (npw, h_diag, h_diag_d, s_diag, s_diag_d)
  !-----------------------------------------------------------------------
  !
  !    add nonlocal pseudopotential term to diagonal part of Hamiltonian
  !    compute the diagonal part of the S matrix
  !
  USE kinds, ONLY: DP
  USE ions_base,  ONLY : nat, ityp, ntyp => nsp
  USE wvfct, ONLY: npwx
  USE lsda_mod, ONLY: current_spin
  USE uspp,  ONLY: deeq, vkb, qq, qq_so, deeq_nc, indv_ijkb0
  USE uspp_param, ONLY: upf, nh, newpseudo
  USE spin_orb, ONLY: lspinorb
  USE noncollin_module, ONLY: noncolin, npol
  !
  USE uspp, ONLY: deeq_d, vkb_d, qq_d, indv_ijkb0_d
  USE cudafor
  !
  IMPLICIT NONE
  INTEGER, INTENT(in) :: npw
  ! number of plane waves
  !
  REAL(dp), INTENT(inout) :: h_diag (npwx,npol)
  ! the diagonal part of the hamiltonian
  REAL(dp), INTENT(out)   :: s_diag (npwx,npol)
  ! the diagonal part of the S matrix
  REAL(dp), device, INTENT(inout) :: h_diag_d (npwx,npol)
  ! the diagonal part of the hamiltonian
  REAL(dp), device, INTENT(out)   :: s_diag_d (npwx,npol)
  ! the diagonal part of the S matrix

  !
  INTEGER :: ikb, jkb, ih, jh, na, nt, ig, ipol
  COMPLEX(DP) :: ps1(2), ps2(2), ar

  INTEGER :: nhnt, i_ijkb0, istat
  COMPLEX(DP) :: ps11, ps21, cv
  REAL(DP) :: timer, sum1, sum2
  CALL start_clock( 'usnldiag' )
  !
  ! initialise s_diag
  !
  !s_diag = 1.d0
  s_diag_d = 1.d0
  h_diag_d = h_diag

  !deeq_d = deeq
  !qq_d = qq

  if ( lspinorb .or. noncolin ) call errore('usnldiag_gpu','lspinorb or noncolin not yet supported on GPU',1)
  !
  !    multiply on projectors
  !

  DO nt = 1, ntyp
  IF ( upf(nt)%tvanp .or.newpseudo (nt) ) THEN
     DO na = 1, nat
        IF (ityp (na) == nt) THEN

           nhnt = nh(nt)
           i_ijkb0 = indv_ijkb0(na)
!$cuf kernel do(1) <<<*,*>>>
              DO ig = 1, npw
                 sum1 = 0.d0
                 sum2 = 0.d0

                 DO ih = 1, nhnt
                    ikb = i_ijkb0 + ih
                    cv = vkb_d(ig, ikb)
                    DO jh = 1, nhnt
                       jkb = i_ijkb0 + jh
                       ps11 = deeq_d (ih, jh, na, current_spin)
                       ps21 =   qq_d (ih, jh, nt)
                       ar   =  cv *conjg( vkb_d (ig, jkb))
                       sum1 = sum1 + dble(ps11 * ar)
                       sum2 = sum2 + dble(ps21 * ar)
                    ENDDO
                 ENDDO
                 h_diag_d (ig,1) = h_diag_d (ig,1) + sum1
                 s_diag_d (ig,1) = s_diag_d (ig,1) + sum2

              ENDDO

        ENDIF
     ENDDO
  ELSE
     DO na = 1, nat
        IF (ityp (na) == nt) THEN

           nhnt = nh(nt)
           i_ijkb0 = indv_ijkb0(na)
!$cuf kernel do(1) <<<*,*>>>
              DO ig = 1, npw
                 sum1 = 0.d0
                 sum2 = 0.d0

                 DO ih = 1, nhnt
                    ikb = i_ijkb0 + ih
                    ps11 = deeq_d (ih, ih, na, current_spin)
                    ps21 =   qq_d (ih, ih, nt)
                    cv = vkb_d(ig, ikb)
                    ar = cv * conjg(cv)
                    sum1 = sum1 + dble(ps11 * ar)
                    sum2 = sum2 + dble(ps21 * ar)
                 ENDDO
                 h_diag_d (ig,1) = h_diag_d (ig,1) + sum1
                 s_diag_d (ig,1) = s_diag_d (ig,1) + sum2

              ENDDO

        ENDIF
     ENDDO
  ENDIF
  ENDDO

  h_diag = h_diag_d
  s_diag = s_diag_d

  CALL stop_clock( 'usnldiag' )

  RETURN
END SUBROUTINE usnldiag_gpu
#endif
