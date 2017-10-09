!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine stres_cc (sigmaxcc)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE atom,                 ONLY : rgrid
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : ntyp => nsp
  USE cell_base,            ONLY : alat, omega, tpiba, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, gstart, nl, g, gg, ngl, gl,igtongl
  USE ener,                 ONLY : etxc, vtxc
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho, rho_core, rhog_core
  USE vlocal,               ONLY : strf
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  implicit none
  ! output
  real(DP) :: sigmaxcc (3, 3)
  ! local variables

  integer :: nt, ng, l, m, ir
  ! counters
  real(DP) :: fact, sigmadiag
  real(DP) , allocatable:: rhocg (:), vxc (:,:)

  sigmaxcc(:,:) = 0.d0
  if ( ANY (upf(1:ntyp)%nlcc) ) goto 15

  return

15 continue
  !
  ! recalculate the exchange-correlation potential
  !
  allocate ( vxc(dfftp%nnr,nspin) )
  call v_xc (rho, rho_core, rhog_core, etxc, vtxc, vxc)
  if (nspin.eq.1.or.nspin.eq.4) then
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
  allocate(rhocg(ngl))
  sigmadiag = 0.0d0
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if
  do nt = 1, ntyp
     if ( upf(nt)%nlcc ) then
        call drhoc (ngl, gl, omega, tpiba2, rgrid(nt)%mesh, rgrid(nt)%r, &
              rgrid(nt)%rab, upf(nt)%rho_atc, rhocg)
        ! diagonal term
        if (gstart==2) sigmadiag = sigmadiag + &
             CONJG(psic (nl(1) ) ) * strf (1,nt) * rhocg (igtongl (1) )
        do ng = gstart, ngm
           sigmadiag = sigmadiag + CONJG(psic (nl (ng) ) ) * &
                strf (ng,nt) * rhocg (igtongl (ng) ) * fact
        enddo

        call deriv_drhoc (ngl, gl, omega, tpiba2, rgrid(nt)%mesh, &
             rgrid(nt)%r, rgrid(nt)%rab, upf(nt)%rho_atc, rhocg)
        ! non diagonal term (g=0 contribution missing)
        do ng = gstart, ngm
           do l = 1, 3
              do m = 1, 3
                 sigmaxcc (l, m) = sigmaxcc (l, m) + CONJG(psic (nl (ng) ) ) &
                      * strf (ng, nt) * rhocg (igtongl (ng) ) * tpiba * &
                      g (l, ng) * g (m, ng) / sqrt (gg (ng) ) * fact
              enddo
           enddo
        enddo
     endif
  enddo

  do l = 1, 3
     sigmaxcc (l, l) = sigmaxcc (l, l) + sigmadiag
  enddo
  call mp_sum(  sigmaxcc, intra_bgrp_comm )
  deallocate (rhocg)
  return
end subroutine stres_cc

#ifdef USE_CUDA
subroutine stres_cc_gpu (sigmaxcc)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE atom,                 ONLY : rgrid
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : ntyp => nsp
  USE cell_base,            ONLY : alat, omega, tpiba, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE gvect,                ONLY : ngm, gstart, nl=>nl_d, g=>g_d, gg=>gg_d, ngl, gl=>gl_d,igtongl=>igtongl_d
  USE ener,                 ONLY : etxc, vtxc
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho, rho_core, rhog_core, rho_core_d, rhog_core_d
  USE vlocal,               ONLY : strf=>strf_d
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic=>psic_d
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE funct,                ONLY : get_iexch, get_icorr, get_igcx, get_igcc
  USE cudafor

  !
  implicit none
  ! output
  real(DP) :: sigmaxcc (3, 3)
  ! local variables

  integer :: nt, ng, l, m, ir
  ! counters
  real(DP) :: fact, sigmadiag
  real(DP) , allocatable,device:: rhocg(:), vxc_d(:,:)
  real(DP) , allocatable:: vxc(:,:)
  real(DP) :: s11,s12,s13,s21,s22,s23,s31,s32,s33
  real(DP) :: tmpf

  integer   :: iexch, icorr, igcx, igcc, blocks
  type(dim3) :: threads

  sigmaxcc(:,:) = 0.d0
  if ( ANY (upf(1:ntyp)%nlcc) ) goto 15

  return

15 continue
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
  if (iexch .eq. 1 .and. icorr .eq. 4 .and. (igcx .eq. 2 .or. igcx .eq. 3) .and.  (igcc .eq. 2 .or. igcc .eq. 4)) then
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

  psic=cmplx(0.,0.,DP)

  if (nspin.eq.1.or.nspin.eq.4) then
     !$cuf kernel do (1) <<<*, *>>>
     do ir = 1, dfftp%nnr
        psic (ir) = vxc_d (ir, 1)
     enddo
  else
     !$cuf kernel do (1) <<<*, *>>>
     do ir = 1, dfftp%nnr
        psic (ir) = 0.5d0 * (vxc_d (ir, 1) + vxc_d (ir, 2) )
     enddo
  endif
  deallocate (vxc_d)

  CALL fwfft ('Dense', psic, dfftp)
  !
  ! psic contains now Vxc(G)
  !
  allocate(rhocg(ngl))

  sigmadiag = 0.0d0
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if

  do nt = 1, ntyp
     if ( upf(nt)%nlcc ) then
        call drhoc_gpu (ngl, gl, omega, tpiba2, rgrid(nt)%mesh, rgrid(nt)%r_d, &
              rgrid(nt)%rab_d, upf(nt)%rho_atc_d, rhocg)
        ! diagonal term
        !$cuf kernel do (1) <<<*, *>>>
        do ng = gstart, ngm
           if (ng==2 .and. gstart==2) sigmadiag = sigmadiag + &
             dble(CONJG(psic(nl(1) ) ) * strf(1,nt)) * rhocg(igtongl (1) )
           sigmadiag = sigmadiag + dble(CONJG(psic(nl(ng) ) ) * &
                strf(ng,nt)) * rhocg(igtongl(ng) ) * fact
        enddo

        call deriv_drhoc_gpu (ngl, gl, omega, tpiba2, rgrid(nt)%mesh, &
             rgrid(nt)%r_d, rgrid(nt)%rab_d, upf(nt)%rho_atc_d, rhocg)

        ! non diagonal term (g=0 contribution missing)
        s11=0.d0
        s12=0.d0
        s13=0.d0
        s21=0.d0
        s22=0.d0
        s23=0.d0
        s31=0.d0
        s32=0.d0
        s33=0.d0
        !$cuf kernel do (1) <<<*, *>>>
        do ng = gstart, ngm
           tmpf=dble( CONJG( psic(nl(ng)) ) &
                * strf(ng, nt)) * rhocg(igtongl(ng) ) * tpiba  &
                / sqrt(gg(ng)) * fact 
           s11 = s11 + tmpf* g(1,ng) * g(1,ng)
           s12 = s12 + tmpf* g(1,ng) * g(2,ng)
           s13 = s13 + tmpf* g(1,ng) * g(3,ng)
           s21 = s21 + tmpf* g(2,ng) * g(1,ng)
           s22 = s22 + tmpf* g(2,ng) * g(2,ng)
           s23 = s23 + tmpf* g(2,ng) * g(3,ng)
           s31 = s31 + tmpf* g(3,ng) * g(1,ng)
           s32 = s32 + tmpf* g(3,ng) * g(2,ng)
           s33 = s33 + tmpf* g(3,ng) * g(3,ng)
        enddo

        sigmaxcc(1,1) = s11
        sigmaxcc(1,2) = s12
        sigmaxcc(1,3) = s13
        sigmaxcc(2,1) = s21
        sigmaxcc(2,2) = s22
        sigmaxcc(2,3) = s23
        sigmaxcc(3,1) = s31
        sigmaxcc(3,2) = s32
        sigmaxcc(3,3) = s33
     endif
  enddo

  do l = 1, 3
     sigmaxcc (l, l) = sigmaxcc (l, l) + sigmadiag
  enddo
  call mp_sum(  sigmaxcc, intra_bgrp_comm )
  deallocate (rhocg)
  return
end subroutine stres_cc_gpu
#endif
