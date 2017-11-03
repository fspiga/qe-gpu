!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine force_lc (forcelc)
  !----------------------------------------------------------------------
  !
  USE kinds
  USE constants, ONLY : tpi
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE esm,       ONLY : esm_force_lc, do_comp_esm, esm_bc

  USE ions_base,     ONLY : nat, ityp, tau
  USE cell_base,     ONLY : alat, omega  
  USE gvect,         ONLY : ngm, gstart, ngl
  USE gvect,         ONLY : nl, igtongl, g
  USE lsda_mod,      ONLY : nspin
  USE vlocal,        ONLY : vloc
  USE scf,           ONLY : rho
  USE control_flags, ONLY : gamma_only
  implicit none

  real(DP), intent(out) :: forcelc (3, nat)
  ! the local-potential contribution to forces on atoms

  integer :: ipol, ig, na
  ! counter on polarizations
  ! counter on G vectors
  ! counter on atoms

  complex(DP), allocatable :: aux (:)
  ! auxiliary space for FFT
  real(DP) :: arg, fact
  !
  ! contribution to the force from the local part of the bare potential
  ! F_loc = Omega \Sum_G n*(G) d V_loc(G)/d R_i
  !
  allocate (aux(dfftp%nnr))
  if ( nspin == 2) then
      aux(:) = CMPLX( rho%of_r(:,1)+rho%of_r(:,2), 0.0_dp, kind=dp )
  else
      aux(:) = CMPLX( rho%of_r(:,1), 0.0_dp, kind=dp )
  end if
  CALL fwfft ('Dense', aux, dfftp)
  !
  !    aux contains now  n(G)
  !
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  do na = 1, nat
     do ipol = 1, 3
        forcelc (ipol, na) = 0.d0
     enddo
     ! contribution from G=0 is zero
     do ig = gstart, ngm
        arg = (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) + &
               g (3, ig) * tau (3, na) ) * tpi
        do ipol = 1, 3
           forcelc (ipol, na) = forcelc (ipol, na) + &
                g (ipol, ig) * vloc (igtongl (ig), ityp (na) ) * &
                (sin(arg)*DBLE(aux(nl(ig))) + cos(arg)*AIMAG(aux(nl(ig))) )
        enddo
     enddo
     do ipol = 1, 3
        forcelc (ipol, na) = fact * forcelc (ipol, na) * omega * tpi / alat
     enddo
  enddo
  IF ( do_comp_esm .and. ( esm_bc .ne. 'pbc' ) ) THEN
     !
     ! ... Perform corrections for ESM method (add long-range part)
     !
     CALL esm_force_lc ( aux, forcelc )
  ENDIF
  !
  call mp_sum(  forcelc, intra_bgrp_comm )
  !
  deallocate (aux)
  return
end subroutine force_lc

#ifdef USE_CUDA
subroutine force_lc_gpu ( forcelc )
  !----------------------------------------------------------------------
  !
  USE kinds
  USE constants, ONLY : tpi
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces, ONLY : fwfft
  USE esm,       ONLY : esm_force_lc, do_comp_esm, esm_bc

  USE ions_base,     ONLY : nat, ityp, tau
  USE cell_base,     ONLY : alat, omega  
  USE gvect,         ONLY : ngm, gstart, ngl
  USE gvect,         ONLY : nl_d, igtongl_d, g_d
  USE lsda_mod,      ONLY : nspin
  USE vlocal,        ONLY : vloc_d
  USE scf,           ONLY : rho
  USE control_flags, ONLY : gamma_only

  implicit none

  real(DP), intent(out) :: forcelc (3, nat)
  ! the local-potential contribution to forces on atoms

  integer :: ipol, ig, na, i, it
  ! counter on polarizations
  ! counter on G vectors
  ! counter on atoms

  complex(DP), allocatable :: aux (:)
  complex(DP), device, allocatable :: aux_d (:)
  ! auxiliary space for FFT
  real(DP) :: arg, fact
  real(DP), device, pointer :: rho_of_r_d(:,:)
  real(DP) :: tau1, tau2, tau3
  real(DP) :: flc1, flc2, flc3
  real(DP) :: prod
  !
  ! contribution to the force from the local part of the bare potential
  ! F_loc = Omega \Sum_G n*(G) d V_loc(G)/d R_i
  !
  forcelc(:,:) = 0.d0

  allocate (aux(dfftp%nnr))
  allocate (aux_d(dfftp%nnr))
  rho%of_r_d = rho%of_r
  rho_of_r_d => rho%of_r_d

  if ( nspin == 2) then
      !$cuf kernel do(1) <<<*, *>>>
      do i = lbound(aux_d, 1), ubound(aux_d, 1)
        aux_d(i) = CMPLX( rho_of_r_d(i,1)+rho_of_r_d(i,2), 0.0_dp, kind=dp )
      end do
  else
      !$cuf kernel do(1) <<<*, *>>>
      do i = lbound(aux_d, 1), ubound(aux_d, 1)
        aux_d(i) = CMPLX( rho_of_r_d(i,1), 0.0_dp, kind=dp )
      end do
  end if
  CALL fwfft ('Dense', aux_d, dfftp)

  !
  !    aux contains now  n(G)
  !
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  do na = 1, nat
     tau1 = tau(1, na)
     tau2 = tau(2, na)
     tau3 = tau(3, na)
     
     flc1 = 0.d0
     flc2 = 0.d0
     flc3 = 0.d0

     it = ityp(na)

     ! contribution from G=0 is zero
     !$cuf kernel do(1) <<<*, *>>>
     do ig = gstart, ngm
        arg = (g_d (1, ig) * tau1 + g_d (2, ig) * tau2 + &
               g_d (3, ig) * tau3) * tpi

        prod = vloc_d (igtongl_d (ig), it ) * DBLE(cmplx(sin(arg), cos(arg), DP) * conjg(aux_d(nl_d(ig))))

        flc1 = flc1 + g_d (1, ig) * prod
        flc2 = flc2 + g_d (2, ig) * prod
        flc3 = flc3 + g_d (3, ig) * prod
     enddo
     forcelc (1, na) = fact * flc1 * omega * tpi / alat
     forcelc (2, na) = fact * flc2 * omega * tpi / alat
     forcelc (3, na) = fact * flc3 * omega * tpi / alat
  enddo
  IF ( do_comp_esm .and. ( esm_bc .ne. 'pbc' ) ) THEN
     !
     ! ... Perform corrections for ESM method (add long-range part)
     !
     aux = aux_d
     CALL esm_force_lc ( aux, forcelc )
  ENDIF
  !
  call mp_sum(  forcelc, intra_bgrp_comm )
  !
  deallocate (aux)
  deallocate (aux_d)
  return
end subroutine force_lc_gpu
#endif
