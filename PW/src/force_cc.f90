!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine force_cc (forcecc)
  !----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE atom,                 ONLY : rgrid
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,            ONLY : alat, omega, tpiba, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, gstart, nl, g, gg, ngl, gl, igtongl
  USE ener,                 ONLY : etxc, vtxc
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho, rho_core, rhog_core
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin
  USE wavefunctions_module, ONLY : psic
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  implicit none
  !
  !   first the dummy variable
  !
  real(DP) :: forcecc (3, nat)
  ! output: the local forces on atoms

  integer :: ipol, ig, ir, nt, na
  ! counter on polarizations
  ! counter on G vectors
  ! counter on FFT grid points
  ! counter on types of atoms
  ! counter on atoms


  real(DP), allocatable :: vxc (:,:), rhocg (:)
  ! exchange-correlation potential
  ! radial fourier trasform of rho core
  real(DP)  ::  arg, fact

  !
  forcecc(:,:) = 0.d0
  if ( ANY ( upf(1:ntyp)%nlcc ) ) go to 15
  return
  !
15 continue
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  !
  ! recalculate the exchange-correlation potential
  !
  allocate ( vxc(dfftp%nnr,nspin) )
  !
  call v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxc)
  !
  psic=(0.0_DP,0.0_DP)
  if (nspin == 1 .or. nspin == 4) then
     do ir = 1, dfftp%nnr
        psic (ir) = vxc (ir, 1)
     enddo
  else
     do ir = 1, dfftp%nnr
        psic (ir) = 0.5d0 * (vxc (ir, 1) + vxc (ir, 2) )
     enddo
  endif
  deallocate (vxc)
  CALL fwfft ('Dense', psic, dfftp)
  !
  ! psic contains now Vxc(G)
  !
  allocate ( rhocg(ngl) )
  !
  ! core correction term: sum on g of omega*ig*exp(-i*r_i*g)*n_core(g)*vxc
  ! g = 0 term gives no contribution
  !
  do nt = 1, ntyp
     if ( upf(nt)%nlcc ) then

        call drhoc (ngl, gl, omega, tpiba2, rgrid(nt)%mesh, rgrid(nt)%r,&
             rgrid(nt)%rab, upf(nt)%rho_atc, rhocg)
        do na = 1, nat
           if (nt.eq.ityp (na) ) then
              do ig = gstart, ngm
                 arg = (g (1, ig) * tau (1, na) + g (2, ig) * tau (2, na) &
                      + g (3, ig) * tau (3, na) ) * tpi
                 do ipol = 1, 3
                    forcecc (ipol, na) = forcecc (ipol, na) + tpiba * omega * &
                         rhocg (igtongl (ig) ) * CONJG(psic (nl (ig) ) ) * &
                         CMPLX( sin (arg), cos (arg) ,kind=DP) * g (ipol, ig) * fact
                 enddo
              enddo
           endif
        enddo
     endif
  enddo
  !
  call mp_sum(  forcecc, intra_bgrp_comm )
  !
  deallocate (rhocg)
  !
  return
end subroutine force_cc

#ifdef USE_CUDA
subroutine force_cc_gpu (forcecc)
  !----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE atom,                 ONLY : rgrid
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,            ONLY : alat, omega, tpiba, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, gstart, nl_d, g_d, ngl, gl_d, igtongl_d
  USE ener,                 ONLY : etxc, vtxc
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho, rho_core, rhog_core, rho_core_d, rhog_core_d
  USE control_flags,        ONLY : gamma_only
  USE noncollin_module,     ONLY : noncolin
  USE wavefunctions_module, ONLY : psic, psic_d
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE funct,                ONLY : get_iexch, get_icorr, get_igcx, get_igcc
  USE cudafor
  !
  implicit none
  !
  !   first the dummy variable
  !
  real(DP) :: forcecc (3, nat)
  ! output: the local forces on atoms

  integer :: ipol, ig, ir, nt, na, istat
  ! counter on polarizations
  ! counter on G vectors
  ! counter on FFT grid points
  ! counter on types of atoms
  ! counter on atoms

  real(DP), allocatable, device :: vxc_d (:,:), rhocg_d (:)
  real(DP), allocatable :: vxc (:,:)
  ! exchange-correlation potential
  ! radial fourier trasform of rho core
  real(DP)  ::  arg, fact, prod
  real(DP)  ::  tau1, tau2, tau3
  real(DP)  ::  fcc1, fcc2, fcc3
  integer   :: iexch, icorr, igcx, igcc


  !
  forcecc(:,:) = 0.d0
!  if ( ANY ( upf(1:ntyp)%nlcc ) ) go to 15
!  return
  !
!15 continue

  if (COUNT (upf(1:ntyp)%nlcc) .eq. 0) then
    return
  endif

  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  !
  ! recalculate the exchange-correlation potential
  !
  allocate ( vxc_d(dfftp%nnr,nspin) )
  !
  iexch = get_iexch
  icorr = get_icorr
  igcx = get_igcx
  igcc = get_igcc

  ! If calling PBE functional configuration, use GPU path
  if (iexch .eq. 1 .and. icorr .eq. 4 .and. (igcx .eq. 2 .or. igcx .eq. 3) .and. (igcc .eq. 2 .or. igcc .eq. 4)) then
    !TODO Might be able to remove some of these copies
    rho_core_d = rho_core
    rhog_core_d = rhog_core
    rho%of_r_d = rho%of_r
    rho%of_g_d = rho%of_g
    CALL v_xc_gpu( rho, rho_core_d, rhog_core_d, etxc, vtxc, vxc_d)

  ! Otherwise fall back to CPU path
  else
    allocate ( vxc(dfftp%nnr,nspin) )
    CALL v_xc( rho, rho_core, rhog_core, etxc, vtxc, vxc)
    vxc_d = vxc
    deallocate(vxc)
  endif
  !
  psic_d=cmplx(0, 0, DP)
  if (nspin == 1 .or. nspin == 4) then
     !$cuf kernel do (1) <<<*, *>>>
     do ir = 1, dfftp%nnr
        psic_d (ir) = vxc_d (ir, 1)
     enddo
  else
     !$cuf kernel do (1) <<<*, *>>>
     do ir = 1, dfftp%nnr
        psic_d (ir) = 0.5d0 * (vxc_d (ir, 1) + vxc_d (ir, 2) )
     enddo
  endif
  deallocate (vxc_d)
  CALL fwfft ('Dense', psic_d, dfftp)
  !
  ! psic contains now Vxc(G)
  !
  allocate ( rhocg_d(ngl) )
  !
  ! core correction term: sum on g of omega*ig*exp(-i*r_i*g)*n_core(g)*vxc
  ! g = 0 term gives no contribution
  !
  do nt = 1, ntyp
     if ( upf(nt)%nlcc ) then
        call drhoc_gpu (ngl, gl_d, omega, tpiba2, rgrid(nt)%mesh,&
            rgrid(nt)%r_d, rgrid(nt)%rab_d, upf(nt)%rho_atc_d, rhocg_d)

        do na = 1, nat
           if (nt.eq.ityp (na) ) then

              tau1 = tau(1, na)
              tau2 = tau(2, na)
              tau3 = tau(3, na)

              fcc1 = 0.d0
              fcc2 = 0.d0
              fcc3 = 0.d0

              !$cuf kernel do (1) <<<*, *>>>
              do ig = gstart, ngm
                 arg = (g_d (1, ig) * tau1 + g_d (2, ig) * tau2 &
                      + g_d (3, ig) * tau3 ) * tpi
                 prod = tpiba * omega * &
                      rhocg_d (igtongl_d (ig) ) * dble(CONJG(psic_d (nl_d (ig) ) ) * &
                      CMPLX( sin (arg), cos (arg), kind=DP))* fact

                 fcc1 = fcc1 + g_d (1, ig) * prod
                 fcc2 = fcc2 + g_d (2, ig) * prod
                 fcc3 = fcc3 + g_d (3, ig) * prod
              enddo

              forcecc(1, na) = forcecc(1, na) + fcc1
              forcecc(2, na) = forcecc(2, na) + fcc2
              forcecc(3, na) = forcecc(3, na) + fcc3
           endif
        enddo
     endif
  enddo
  !
  call mp_sum(  forcecc, intra_bgrp_comm )
  !
  deallocate (rhocg_d)
  !
  return
end subroutine force_cc_gpu
#endif
