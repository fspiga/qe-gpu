!
! Copyright (C) 2004-2016 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-------------------------------------------------------------------
#ifdef USE_CUDA
module funct_gpu
USE kinds, ONLY : DP
use cudafor
contains
!-----------------------------------------------------------------------
! subroutines ported from  funct.f90
!-----------------------------------------------------------------------
attributes(device) subroutine xc_dev (iexch, icorr, rho, ex, ec, vx, vc)
  !-----------------------------------------------------------------------
  !     lda exchange and correlation functionals - Hartree a.u.
  !
  !     exchange   :  Slater, relativistic Slater
  !     correlation:  Ceperley-Alder (Perdew-Zunger parameters)
  !                   Vosko-Wilk-Nusair
  !                   Lee-Yang-Parr
  !                   Perdew-Wang
  !                   Wigner
  !                   Hedin-Lundqvist
  !                   Ortiz-Ballone (Perdew-Zunger formula)
  !                   Ortiz-Ballone (Perdew-Wang formula)
  !                   Gunnarsson-Lundqvist
  !
  !     input : iexch, icorr, rho=rho(r)
  !     definitions: E_x = \int E_x(rho) dr, E_x(rho) = rho\epsilon_c(rho)
  !                  same for correlation
  !     output: ex = \epsilon_x(rho) ( NOT E_x(rho) )
  !             vx = dE_x(rho)/drho  ( NOT d\epsilon_x(rho)/drho )
  !             ec, vc as above for correlation
  !
  USE kinds, ONLY : DP
  implicit none

  integer, value  :: iexch, icorr
  real(DP), value :: rho
  real(DP), device, intent(out) :: ec, vc, ex, vx
  real(DP) :: ec__, vc__
  !
  real(DP), parameter :: small = 1.E-10_DP,  third = 1.0_DP / 3.0_DP, &
       pi34 = 0.6203504908994_DP  ! pi34=(3/4pi)^(1/3)
  real(DP) :: rs
  !
  if (rho <= small) then
     ec = 0.0_DP
     vc = 0.0_DP
     ex = 0.0_DP
     vx = 0.0_DP
     return
  else
     rs = pi34 / rho**third
     ! rs as in the theory of metals: rs=(3/(4pi rho))^(1/3)
  endif
  !..exchange
  if (iexch == 1) THEN             !  'sla'
     call slater_dev (rs, ex, vx)
  !ELSEIF (iexch == 2) THEN         !  'sl1'
  !   call slater1(rs, ex, vx)
  !ELSEIF (iexch == 3) THEN         !  'rxc'
  !   CALL slater_rxc(rs, ex, vx)
  !ELSEIF ((iexch == 4).or.(iexch==5)) THEN  ! 'oep','hf'
  !   IF (exx_started) then
  !      ex = 0.0_DP
  !      vx = 0.0_DP
  !   else
  !      call slater (rs, ex, vx)
  !   endif
  !ELSEIF (iexch == 6) THEN         !  'pb0x'
  !   CALL slater(rs, ex, vx)
  !   if (exx_started) then
  !      ex = (1.0_DP - exx_fraction) * ex 
  !      vx = (1.0_DP - exx_fraction) * vx 
  !   end if
  !ELSEIF (iexch == 7) THEN         !  'B3LYP'
  !   CALL slater(rs, ex, vx)
  !   if (exx_started) then
  !      ex = 0.8_DP * ex 
  !      vx = 0.8_DP * vx 
  !   end if
  !ELSEIF (iexch == 8) THEN         !  'sla+kzk'
  !   if (.NOT. finite_size_cell_volume_set) call errore ('XC',&
  !        'finite size corrected exchange used w/o initialization',1)
  !   call slaterKZK (rs, ex, vx, finite_size_cell_volume)
  !   !
  !ELSEIF (iexch == 9) THEN         !  'X3LYP'
  !   CALL slater(rs, ex, vx)
  !   if (exx_started) then
  !      ex = 0.782_DP * ex 
  !      vx = 0.782_DP * vx 
  !   end if
  else
     ex = 0.0_DP
     vx = 0.0_DP
  endif
  !!..correlation
  !if (icorr == 1) then
  !   call pz_dev (rs, 1, ec, vc)
  !elseif (icorr == 2) then
  !   call vwn (rs, ec, vc)
  !elseif (icorr == 3) then
  !   call lyp (rs, ec, vc)
  if (icorr == 4) then
  !elseif (icorr == 4) then
     call pw_dev (rs, 1, ec, vc)
  !elseif (icorr == 5) then
  !   call wigner (rs, ec, vc)
  !elseif (icorr == 6) then
  !   call hl (rs, ec, vc)
  !elseif (icorr == 7) then
  !   call pz (rs, 2, ec, vc)
  !elseif (icorr == 8) then
  !   call pw (rs, 2, ec, vc)
  !elseif (icorr == 9) then
  !   call gl (rs, ec, vc)
  !elseif (icorr ==10) then
  !   if (.NOT. finite_size_cell_volume_set) call errore ('XC',&
  !        'finite size corrected correlation used w/o initialization',1)
  !   call pzKZK (rs, ec, vc, finite_size_cell_volume)
  !elseif (icorr ==11) then
  !   call vwn1_rpa (rs, ec, vc)
  !elseif (icorr ==12) then  ! 'B3LYP'
  !   call vwn (rs, ec, vc)
  !   ec = 0.19_DP * ec
  !   vc = 0.19_DP * vc
  !   call lyp( rs, ec__, vc__ )
  !   ec = ec + 0.81_DP * ec__
  !   vc = vc + 0.81_DP * vc__
  !elseif (icorr ==13) then  ! 'B3LYP-V1R'
  !   call vwn1_rpa (rs, ec, vc)
  !   ec = 0.19_DP * ec
  !   vc = 0.19_DP * vc
  !   call lyp( rs, ec__, vc__ )
  !   ec = ec + 0.81_DP * ec__
  !   vc = vc + 0.81_DP * vc__
  !elseif (icorr ==14) then  ! 'X3LYP'
  !   call vwn1_rpa (rs, ec, vc)
  !   ec = 0.129_DP * ec
  !   vc = 0.129_DP * vc
  !   call lyp( rs, ec__, vc__ )
  !   ec = ec + 0.871_DP * ec__
  !   vc = vc + 0.871_DP * vc__
  else
     ec = 0.0_DP
     vc = 0.0_DP
  endif
  !
  return
end subroutine xc_dev
!-----------------------------------------------------------------------
attributes(device) subroutine gcxc_dev (igcx, igcc, rho, grho, sx, sc, v1x, v2x, v1c, v2c)
  !-----------------------------------------------------------------------
  !     gradient corrections for exchange and correlation - Hartree a.u.
  !     See comments at the beginning of module for implemented cases
  !
  !     input:  igcx, igcc, rho, grho=|\nabla rho|^2
  !     definition:  E_x = \int E_x(rho,grho) dr
  !     output: sx = E_x(rho,grho)
  !             v1x= D(E_x)/D(rho)
  !             v2x= D(E_x)/D( D rho/D r_alpha ) / |\nabla rho|
  !             sc, v1c, v2c as above for correlation
  !
  implicit none

  integer, value :: igcx, igcc
  real(DP), value :: rho, grho
  real(DP), device, intent(out) :: sx, sc, v1x, v2x, v1c, v2c
  real(DP) :: sx__,v1x__, v2x__
  real(DP) :: sxsr, v1xsr, v2xsr
  real(DP), parameter:: small = 1.E-10_DP

  ! exchange
  if (rho <= small) then
     sx = 0.0_DP
     v1x = 0.0_DP
     v2x = 0.0_DP
  !elseif (igcx == 1) then
  !   call becke88 (rho, grho, sx, v1x, v2x)
  elseif (igcx == 2) then
     call ggax_dev (rho, grho, sx, v1x, v2x)
  elseif (igcx == 3) then
     call pbex_dev (rho, grho, 1, sx, v1x, v2x)
  !elseif (igcx == 4) then
  !   call pbex (rho, grho, 2, sx, v1x, v2x)
  !elseif (igcx == 5 .and. igcc == 5) then
  !   call hcth(rho, grho, sx, v1x, v2x)
  !elseif (igcx == 6) then
  !   call optx (rho, grho, sx, v1x, v2x)
  !! case igcx == 7 (meta-GGA) must be treated in a separate call to another
  !! routine: needs kinetic energy density in addition to rho and grad rho
  !elseif (igcx == 8) then ! 'PBE0'
  !   call pbex (rho, grho, 1, sx, v1x, v2x)
  !   if (exx_started) then
  !      sx  = (1.0_DP - exx_fraction) * sx
  !      v1x = (1.0_DP - exx_fraction) * v1x
  !      v2x = (1.0_DP - exx_fraction) * v2x
  !   end if
  !elseif (igcx == 9) then ! 'B3LYP'
  !   call becke88 (rho, grho, sx, v1x, v2x)
  !   if (exx_started) then
  !      sx  = 0.72_DP * sx
  !      v1x = 0.72_DP * v1x
  !      v2x = 0.72_DP * v2x
  !   end if
  !elseif (igcx ==10) then ! 'pbesol'
  !   call pbex (rho, grho, 3, sx, v1x, v2x)
  !elseif (igcx ==11) then ! 'wc'
  !   call wcx (rho, grho, sx, v1x, v2x)
  !elseif (igcx ==12) then ! 'pbexsr'
  !   call pbex (rho, grho, 1, sx, v1x, v2x)
  !   if(exx_started) then
  !     call pbexsr (rho, grho, sxsr, v1xsr, v2xsr, screening_parameter)
  !     sx = sx - exx_fraction * sxsr
  !     v1x = v1x - exx_fraction * v1xsr
  !     v2x = v2x - exx_fraction * v2xsr
  !   endif 
  !elseif (igcx ==13) then ! 'rPW86'
  !   call rPW86 (rho, grho, sx, v1x, v2x)
  !elseif (igcx ==16) then ! 'C09x'
  !   call c09x (rho, grho, sx, v1x, v2x)
  !elseif (igcx ==17) then ! 'sogga'
  !   call sogga(rho, grho, sx, v1x, v2x)
  !elseif (igcx ==19) then ! 'pbeq2d'
  !   call pbex (rho, grho, 4, sx, v1x, v2x)
  !elseif (igcx ==20) then ! 'gau-pbe'
  !   call pbex (rho, grho, 1, sx, v1x, v2x)
  !   if(exx_started) then
  !     call pbexgau (rho, grho, sxsr, v1xsr, v2xsr, gau_parameter)
  !     sx = sx - exx_fraction * sxsr
  !     v1x = v1x - exx_fraction * v1xsr
  !     v2x = v2x - exx_fraction * v2xsr
  !   endif
  !elseif (igcx == 21) then ! 'pw86'
  !   call pw86 (rho, grho, sx, v1x, v2x)
  !elseif (igcx == 22) then ! 'b86b'
  !   call becke86b (rho, grho, sx, v1x, v2x)
  !   ! call b86b (rho, grho, 1, sx, v1x, v2x)
  !elseif (igcx == 23) then ! 'optB88'
  !   call pbex (rho, grho, 5, sx, v1x, v2x)
  !elseif (igcx == 24) then ! 'optB86b'
  !   call pbex (rho, grho, 6, sx, v1x, v2x)
  !   ! call b86b (rho, grho, 2, sx, v1x, v2x)
  !elseif (igcx == 25) then ! 'ev93'
  !   call pbex (rho, grho, 7, sx, v1x, v2x)
  !elseif (igcx == 26) then ! 'b86r'
  !   call b86b (rho, grho, 3, sx, v1x, v2x)
  !elseif (igcx == 27) then ! 'cx13'
  !   call cx13 (rho, grho, sx, v1x, v2x)
  !elseif (igcx == 28) then ! 'X3LYP'
  !   call becke88 (rho, grho, sx, v1x, v2x)
  !   call pbex (rho, grho, 1, sx__, v1x__, v2x__)
  !   if (exx_started) then
  !      sx  = real(0.765*0.709,DP) * sx
  !      v1x = real(0.765*0.709,DP) * v1x
  !      v2x = real(0.765*0.709,DP) * v2x
  !      sx  = sx  + real(0.235*0.709) * sx__
  !      v1x = v1x + real(0.235*0.709) * v1x__
  !      v2x = v2x + real(0.235*0.709) * v2x__
  !   end if
  else
     sx = 0.0_DP
     v1x = 0.0_DP
     v2x = 0.0_DP
  endif
  ! correlation
  if (rho.le.small) then
     sc = 0.0_DP
     v1c = 0.0_DP
     v2c = 0.0_DP
  !elseif (igcc == 1) then
  !   call perdew86 (rho, grho, sc, v1c, v2c)
  elseif (igcc == 2) then
     call ggac_dev (rho, grho, sc, v1c, v2c)
  !elseif (igcc == 3) then
  !   call glyp (rho, grho, sc, v1c, v2c)
  elseif (igcc == 4) then
     call pbec_dev (rho, grho, 1, sc, v1c, v2c)
  ! igcc == 5 (HCTH) is calculated together with case igcx=5
  ! igcc == 6 (meta-GGA) is treated in a different routine
  !elseif (igcc == 7) then !'B3LYP'
  !   call glyp (rho, grho, sc, v1c, v2c)
  !   if (exx_started) then
  !      sc  = 0.81_DP * sc
  !      v1c = 0.81_DP * v1c
  !      v2c = 0.81_DP * v2c
  !   end if
  !elseif (igcc == 8) then ! 'PBEsol'
  !   call pbec (rho, grho, 2, sc, v1c, v2c)
  !! igcc == 9 set to 5, back-compatibility
  !! igcc ==10 set to 6, back-compatibility
  !! igcc ==11 M06L calculated in another routine
  !else if (igcc == 12) then ! 'Q2D'
  !   call pbec (rho, grho, 3, sc, v1c, v2c)
  !elseif (igcc == 13) then !'X3LYP'
  !   call glyp (rho, grho, sc, v1c, v2c)
  !   if (exx_started) then
  !      sc  = 0.871_DP * sc
  !      v1c = 0.871_DP * v1c
  !      v2c = 0.871_DP * v2c
  !   end if
  else
     sc = 0.0_DP
     v1c = 0.0_DP
     v2c = 0.0_DP
  endif
  !
  return
end subroutine gcxc_dev

attributes(device) subroutine xc_spin_dev (iexch, icorr, rho, zeta, ex, ec, vxup, vxdw, vcup, vcdw)
  !-----------------------------------------------------------------------
  !     lsd exchange and correlation functionals - Hartree a.u.
  !
  !     exchange  :  Slater (alpha=2/3)
  !     correlation: Ceperley & Alder (Perdew-Zunger parameters)
  !                  Perdew & Wang
  !
  !     input : iexch, icorr
  !             rho = rhoup(r)+rhodw(r)
  !             zeta=(rhoup(r)-rhodw(r))/rho
  !
  implicit none

  integer, value :: iexch, icorr
  real(DP), value :: rho, zeta
  real(DP), device, intent(out)  :: ex, ec, vxup, vxdw, vcup, vcdw
  real(DP) :: ec__, vcup__, vcdw__
  !
  real(DP), parameter :: small= 1.E-10_DP, third = 1.0_DP/3.0_DP, &
       pi34= 0.6203504908994_DP ! pi34=(3/4pi)^(1/3)
  real(DP) :: rs
  !
  if (rho <= small) then
     ec = 0.0_DP
     vcup = 0.0_DP
     vcdw = 0.0_DP
     ex = 0.0_DP
     vxup = 0.0_DP
     vxdw = 0.0_DP
     return
  else
     rs = pi34 / rho**third
  endif
  !..exchange
  IF (iexch == 1) THEN      ! 'sla'
     call slater_spin_dev (rho, zeta, ex, vxup, vxdw)
  !ELSEIF (iexch == 2) THEN  ! 'sl1'
  !   call slater1_spin (rho, zeta, ex, vxup, vxdw)
  !ELSEIF (iexch == 3) THEN  ! 'rxc'
  !   call slater_rxc_spin ( rho, zeta, ex, vxup, vxdw )
  !ELSEIF ((iexch == 4).or.(iexch==5)) THEN  ! 'oep','hf'
  !   IF (exx_started) then
  !      ex   = 0.0_DP
  !      vxup = 0.0_DP 
  !      vxdw = 0.0_DP 
  !   else
  !      call slater_spin (rho, zeta, ex, vxup, vxdw)
  !   endif
  !ELSEIF (iexch == 6) THEN  ! 'pb0x'
  !   call slater_spin (rho, zeta, ex, vxup, vxdw)
  !   if (exx_started) then
  !      ex   = (1.0_DP - exx_fraction) * ex
  !      vxup = (1.0_DP - exx_fraction) * vxup 
  !      vxdw = (1.0_DP - exx_fraction) * vxdw 
  !   end if
  !ELSEIF (iexch == 7) THEN  ! 'B3LYP'
  !   call slater_spin (rho, zeta, ex, vxup, vxdw)
  !   if (exx_started) then
  !      ex   = 0.8_DP * ex
  !      vxup = 0.8_DP * vxup 
  !      vxdw = 0.8_DP * vxdw 
  !   end if
  ELSE
     ex = 0.0_DP
     vxup = 0.0_DP
     vxdw = 0.0_DP
  ENDIF
  !..correlation
  if (icorr == 0) then
     ec = 0.0_DP
     vcup = 0.0_DP
     vcdw = 0.0_DP
  !elseif (icorr == 1) then
  !   call pz_spin (rs, zeta, ec, vcup, vcdw)
  !elseif (icorr == 2) then
  !   call vwn_spin (rs, zeta, ec, vcup, vcdw)
  !elseif (icorr == 3) then
  !   call lsd_lyp (rho, zeta, ec, vcup, vcdw) ! from CP/FPMD (more_functionals)
  elseif (icorr == 4) then
     call pw_spin_dev (rs, zeta, ec, vcup, vcdw)
  !elseif (icorr == 12) then ! 'B3LYP'
  !   call vwn_spin (rs, zeta, ec, vcup, vcdw)
  !   ec = 0.19_DP * ec
  !   vcup = 0.19_DP * vcup
  !   vcdw = 0.19_DP * vcdw
  !   call lsd_lyp (rho, zeta, ec__, vcup__, vcdw__) ! from CP/FPMD (more_functionals)
  !   ec = ec + 0.81_DP * ec__
  !   vcup = vcup + 0.81_DP * vcup__
  !   vcdw = vcdw + 0.81_DP * vcdw__
  !elseif (icorr == 13) then   ! 'B3LYP-V1R'
  !   call vwn1_rpa_spin (rs, zeta, ec, vcup, vcdw)
  !   ec = 0.19_DP * ec
  !   vcup = 0.19_DP * vcup
  !   vcdw = 0.19_DP * vcdw
  !   call lsd_lyp (rho, zeta, ec__, vcup__, vcdw__) ! from CP/FPMD (more_functionals)
  !   ec = ec + 0.81_DP * ec__
  !   vcup = vcup + 0.81_DP * vcup__
  !   vcdw = vcdw + 0.81_DP * vcdw__
  !else
  !   call errore ('lsda_functional (xc_spin)', 'not implemented', icorr)
  endif
  !
  return
end subroutine xc_spin_dev

attributes(device) subroutine gcx_spin_dev (igcx, rhoup, rhodw, grhoup2, grhodw2, &
                     sx, v1xup, v1xdw, v2xup, v2xdw)
  !-----------------------------------------------------------------------
  !     gradient corrections for exchange - Hartree a.u.
  !
  implicit none
  !
  !     dummy arguments
  !
  integer, value :: igcx
  real(DP), value :: rhoup, rhodw, grhoup2, grhodw2
  real(DP), device, intent(out) :: sx, v1xup, v1xdw, &
       v2xup, v2xdw
  ! up and down charge
  ! up and down gradient of the charge
  ! exchange and correlation energies
  ! derivatives of exchange wr. rho
  ! derivatives of exchange wr. grho
  !
  real(DP) :: sxsr, v1xupsr, v2xupsr, v1xdwsr, v2xdwsr
  real(DP), parameter :: small = 1.E-10_DP
  real(DP) :: rho, sxup, sxdw
  integer :: iflag
  !
  !
  ! exchange
  rho = rhoup + rhodw
  if (rho <= small .or. igcx == 0) then
     sx = 0.0_DP
     v1xup = 0.0_DP
     v2xup = 0.0_DP
     v1xdw = 0.0_DP
     v2xdw = 0.0_DP
  !elseif (igcx == 1) then
  !   if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
  !      call becke88_spin (rhoup, grhoup2, sxup, v1xup, v2xup)
  !   else
  !      sxup = 0.0_DP
  !      v1xup = 0.0_DP
  !      v2xup = 0.0_DP
  !   endif
  !   if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
  !      call becke88_spin (rhodw, grhodw2, sxdw, v1xdw, v2xdw)
  !   else
  !      sxdw = 0.0_DP
  !      v1xdw = 0.0_DP
  !      v2xdw = 0.0_DP
  !   endif
  !   sx = sxup + sxdw
  elseif (igcx == 2) then
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call ggax_dev (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call ggax_dev (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw
  elseif (igcx == 3 .or. igcx == 4 .or. igcx == 8 .or. &
          igcx == 10 .or. igcx == 12 .or. igcx == 20 .or. igcx == 25) then
     ! igcx=3: PBE, igcx=4: revised PBE, igcx=8: PBE0, igcx=10: PBEsol
     ! igcx=12: HSE,  igcx=20: gau-pbe, igcx=25: ev93
     if (igcx == 4) then
        iflag = 2
     elseif (igcx == 10) then
        iflag = 3
     elseif (igcx == 25) then
        iflag = 7
     else
        iflag = 1
     endif
     if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
        call pbex_dev (2.0_DP * rhoup, 4.0_DP * grhoup2, iflag, sxup, v1xup, v2xup)
     else
        sxup = 0.0_DP
        v1xup = 0.0_DP
        v2xup = 0.0_DP
     endif
     if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
        call pbex_dev (2.0_DP * rhodw, 4.0_DP * grhodw2, iflag, sxdw, v1xdw, v2xdw)
     else
        sxdw = 0.0_DP
        v1xdw = 0.0_DP
        v2xdw = 0.0_DP
     endif
     sx = 0.5_DP * (sxup + sxdw)
     v2xup = 2.0_DP * v2xup
     v2xdw = 2.0_DP * v2xdw
     !if (igcx == 8 .and. exx_started ) then
     !  sx = (1.0_DP - exx_fraction) * sx
     !  v1xup = (1.0_DP - exx_fraction) * v1xup
     !  v1xdw = (1.0_DP - exx_fraction) * v1xdw
     !  v2xup = (1.0_DP - exx_fraction) * v2xup
     !  v2xdw = (1.0_DP - exx_fraction) * v2xdw
     !end if
     !if (igcx == 12 .and. exx_started ) then

     !   call pbexsr_lsd (rhoup, rhodw, grhoup2, grhodw2, sxsr,  &
     !                    v1xupsr, v2xupsr, v1xdwsr, v2xdwsr, &
     !                    screening_parameter)
     !   sx  = sx - exx_fraction*sxsr
     !   v1xup = v1xup - exx_fraction*v1xupsr
     !   v2xup = v2xup - exx_fraction*v2xupsr
     !   v1xdw = v1xdw - exx_fraction*v1xdwsr
     !   v2xdw = v2xdw - exx_fraction*v2xdwsr
     !end if

     !if (igcx == 20 .and. exx_started ) then
     !   ! gau-pbe
     !   call pbexgau_lsd (rhoup, rhodw, grhoup2, grhodw2, sxsr,  &
     !                    v1xupsr, v2xupsr, v1xdwsr, v2xdwsr, &
     !                    gau_parameter)
     !   sx  = sx - exx_fraction*sxsr
     !   v1xup = v1xup - exx_fraction*v1xupsr
     !   v2xup = v2xup - exx_fraction*v2xupsr
     !   v1xdw = v1xdw - exx_fraction*v1xdwsr
     !   v2xdw = v2xdw - exx_fraction*v2xdwsr
     !end if

  !elseif (igcx == 9) then
  !   if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
  !      call becke88_spin (rhoup, grhoup2, sxup, v1xup, v2xup)
  !   else
  !      sxup = 0.0_DP
  !      v1xup = 0.0_DP
  !      v2xup = 0.0_DP
  !   endif
  !   if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
  !      call becke88_spin (rhodw, grhodw2, sxdw, v1xdw, v2xdw)
  !   else
  !      sxdw = 0.0_DP
  !      v1xdw = 0.0_DP
  !      v2xdw = 0.0_DP
  !   endif
  !   sx = sxup + sxdw

  !   if (exx_started ) then
  !     sx = 0.72_DP * sx
  !     v1xup = 0.72_DP * v1xup
  !     v1xdw = 0.72_DP * v1xdw
  !     v2xup = 0.72_DP * v2xup
  !     v2xdw = 0.72_DP * v2xdw
  !   end if

  !elseif (igcx == 11) then ! 'Wu-Cohen'
  !   if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
  !      call wcx (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
  !   else
  !      sxup = 0.0_DP
  !      v1xup = 0.0_DP
  !      v2xup = 0.0_DP
  !   endif
  !   if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
  !      call wcx (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
  !   else
  !      sxdw = 0.0_DP
  !      v1xdw = 0.0_DP
  !      v2xdw = 0.0_DP
  !   endif
  !   sx = 0.5_DP * (sxup + sxdw)
  !   v2xup = 2.0_DP * v2xup
  !   v2xdw = 2.0_DP * v2xdw

  !elseif (igcx == 13) then ! 'revised PW86 for vdw-df2'
  !   if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
  !      call rPW86 (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
  !   else
  !      sxup = 0.0_DP
  !      v1xup = 0.0_DP
  !      v2xup = 0.0_DP
  !   endif
  !   if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
  !      call rPW86 (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
  !   else
  !      sxdw = 0.0_DP
  !      v1xdw = 0.0_DP
  !      v2xdw = 0.0_DP
  !   endif
  !   sx = 0.5_DP * (sxup + sxdw)
  !   v2xup = 2.0_DP * v2xup
  !   v2xdw = 2.0_DP * v2xdw

  !elseif (igcx == 16) then ! 'c09x for vdw-df-c09.'
  !   if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
  !      call c09x (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
  !   else
  !      sxup = 0.0_DP
  !      v1xup = 0.0_DP
  !      v2xup = 0.0_DP
  !   endif
  !   if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
  !      call c09x (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
  !   else
  !      sxdw = 0.0_DP
  !      v1xdw = 0.0_DP
  !      v2xdw = 0.0_DP
  !   endif
  !   sx = 0.5_DP * (sxup + sxdw)
  !   v2xup = 2.0_DP * v2xup
  !   v2xdw = 2.0_DP * v2xdw

  !elseif (igcx == 21) then ! 'PW86'
  !   if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
  !      call pw86 (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
  !   else
  !      sxup = 0.0_DP
  !      v1xup = 0.0_DP
  !      v2xup = 0.0_DP
  !   endif
  !   if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
  !      call pw86 (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
  !   else
  !      sxdw = 0.0_DP
  !      v1xdw = 0.0_DP
  !      v2xdw = 0.0_DP
  !   endif
  !   sx = 0.5_DP * (sxup + sxdw)
  !   v2xup = 2.0_DP * v2xup
  !   v2xdw = 2.0_DP * v2xdw

  !elseif (igcx == 22) then ! 'B86B'
  !   if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
  !      call becke86b (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
  !   else
  !      sxup = 0.0_DP
  !      v1xup = 0.0_DP
  !      v2xup = 0.0_DP
  !   endif
  !   if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
  !      call becke86b (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
  !   else
  !      sxdw = 0.0_DP
  !      v1xdw = 0.0_DP
  !      v2xdw = 0.0_DP
  !   endif
  !   sx = 0.5_DP * (sxup + sxdw)
  !   v2xup = 2.0_DP * v2xup
  !   v2xdw = 2.0_DP * v2xdw

  ! elseif (igcx == 26) then ! 'B86R for rev-vdW-DF2'
  !   if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
  !      call b86b (2.0_DP * rhoup, 4.0_DP * grhoup2, 3, sxup, v1xup, v2xup)
  !   else
  !      sxup = 0.0_DP
  !      v1xup = 0.0_DP
  !      v2xup = 0.0_DP
  !   endif
  !   if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
  !      call b86b (2.0_DP * rhodw, 4.0_DP * grhodw2, 3, sxdw, v1xdw, v2xdw)
  !   else
  !      sxdw = 0.0_DP
  !      v1xdw = 0.0_DP
  !      v2xdw = 0.0_DP
  !   endif
  !   sx = 0.5_DP * (sxup + sxdw)
  !   v2xup = 2.0_DP * v2xup
  !   v2xdw = 2.0_DP * v2xdw

  !elseif (igcx == 27) then ! 'cx13 for vdw-df-cx'
  !   if (rhoup > small .and. sqrt (abs (grhoup2) ) > small) then
  !      call cx13 (2.0_DP * rhoup, 4.0_DP * grhoup2, sxup, v1xup, v2xup)
  !   else
  !      sxup = 0.0_DP
  !      v1xup = 0.0_DP
  !      v2xup = 0.0_DP
  !   endif
  !   if (rhodw > small .and. sqrt (abs (grhodw2) ) > small) then
  !      call cx13 (2.0_DP * rhodw, 4.0_DP * grhodw2, sxdw, v1xdw, v2xdw)
  !   else
  !      sxdw = 0.0_DP
  !      v1xdw = 0.0_DP
  !      v2xdw = 0.0_DP
  !   endif
  !   sx = 0.5_DP * (sxup + sxdw)
  !   v2xup = 2.0_DP * v2xup
  !   v2xdw = 2.0_DP * v2xdw


  !! case igcx == 5 (HCTH) and 6 (OPTX) not implemented
  !! case igcx == 7 (meta-GGA) must be treated in a separate call to another
  !! routine: needs kinetic energy density in addition to rho and grad rho

  !else
  !   call errore ('gcx_spin', 'not implemented', igcx)
  endif
  !
  return
end subroutine gcx_spin_dev

attributes(device) subroutine gcc_spin_dev (igcc, rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  !-----------------------------------------------------------------------
  !     gradient corrections for correlations - Hartree a.u.
  !     Implemented:  Perdew86, GGA (PW91), PBE
  !
  implicit none
  !
  !     dummy arguments
  !
  integer, value :: igcc
  real(DP), value :: rho, zeta, grho
  real(DP), device, intent(out) :: sc, v1cup, v1cdw, v2c
  ! the total charge
  ! the magnetization
  ! the gradient of the charge squared
  ! exchange and correlation energies
  ! derivatives of correlation wr. rho
  ! derivatives of correlation wr. grho

  real(DP), parameter :: small = 1.E-10_DP, epsr=1.E-6_DP
  !
  if ( abs(zeta) > 1.0_DP ) then
     sc = 0.0_DP
     v1cup = 0.0_DP
     v1cdw = 0.0_DP
     v2c = 0.0_DP
     return
  else
     !
     ! ... ( - 1.0 + epsr )  <  zeta  <  ( 1.0 - epsr )
     zeta = SIGN( MIN( ABS( zeta ), ( 1.0_DP - epsr ) ) , zeta )
  endif

  if (igcc == 0 .or. rho <= small .or. sqrt(abs(grho)) <= small) then
     sc = 0.0_DP
     v1cup = 0.0_DP
     v1cdw = 0.0_DP
     v2c = 0.0_DP
  !elseif (igcc == 1) then
  !   call perdew86_spin (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  elseif (igcc == 2) then
     call ggac_spin_dev (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  elseif (igcc == 4) then
     call pbec_spin_dev (rho, zeta, grho, 1, sc, v1cup, v1cdw, v2c)
  !elseif (igcc == 8) then
  !   call pbec_spin (rho, zeta, grho, 2, sc, v1cup, v1cdw, v2c)
  !else
  !   call errore ('lsda_functionals (gcc_spin)', 'not implemented', igcc)
  endif
  !
  return
end subroutine gcc_spin_dev


!-----------------------------------------------------------------------
! subroutines ported from functionals.f90
!-----------------------------------------------------------------------

attributes(device) subroutine slater_dev (rs, ex, vx)
  !-----------------------------------------------------------------------
  !        Slater exchange with alpha=2/3
  !
  USE kinds, ONLY : DP
  implicit none
  real(dp), value :: rs
  real(dp), device, intent(out):: ex, vx
  real(dp), parameter  :: f= -0.687247939924714d0, alpha = 2.0d0/3.0d0
  ! f = -9/8*(3/2pi)^(2/3)
  !
  ex = f * alpha / rs
  vx = 4.d0 / 3.d0 * f * alpha / rs
  !
  return
end subroutine slater_dev

attributes(device) subroutine pw_dev (rs, iflag, ec, vc)
  !-----------------------------------------------------------------------
  !     iflag=1: J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
  !     iflag=2: G. Ortiz and P. Ballone, PRB 50, 1391 (1994)
  !
  implicit none
  real(dp), value :: rs
  real(dp), device, intent(out):: ec, vc
  integer, value :: iflag 
  real(DP) :: a, b1, b2, c0, c1, c2, c3, d0, d1
  parameter (a = 0.031091d0, b1 = 7.5957d0, b2 = 3.5876d0, c0 = a, &
       c1 = 0.046644d0, c2 = 0.00664d0, c3 = 0.01043d0, d0 = 0.4335d0, &
       d1 = 1.4408d0)
  real(DP) :: lnrs, rs12, rs32, rs2, om, dom, olog
  real(DP) :: a1 (2), b3 (2), b4 (2)
  data a1 / 0.21370d0, 0.026481d0 /, b3 / 1.6382d0, -0.46647d0 /, &
       b4 / 0.49294d0, 0.13354d0 /
  !
  ! high- and low-density formulae implemented but not used in PW case
  ! (reason: inconsistencies in PBE/PW91 functionals)
  !
  if (rs.lt.1d0.and.iflag.eq.2) then
     ! high density formula
     lnrs = log (rs)
     ec = c0 * lnrs - c1 + c2 * rs * lnrs - c3 * rs
     vc = c0 * lnrs - (c1 + c0 / 3.d0) + 2.d0 / 3.d0 * c2 * rs * &
          lnrs - (2.d0 * c3 + c2) / 3.d0 * rs
  elseif (rs.gt.100.d0.and.iflag.eq.2) then
     ! low density formula
     ec = - d0 / rs + d1 / rs**1.5d0
     vc = - 4.d0 / 3.d0 * d0 / rs + 1.5d0 * d1 / rs**1.5d0
  else
     ! interpolation formula
     rs12 = sqrt (rs)
     rs32 = rs * rs12
     rs2 = rs**2
     om = 2.d0 * a * (b1 * rs12 + b2 * rs + b3 (iflag) * rs32 + b4 ( &
          iflag) * rs2)
     dom = 2.d0 * a * (0.5d0 * b1 * rs12 + b2 * rs + 1.5d0 * b3 ( &
          iflag) * rs32 + 2.d0 * b4 (iflag) * rs2)
     olog = log (1.d0 + 1.0d0 / om)
     ec = - 2.d0 * a * (1.d0 + a1 (iflag) * rs) * olog
     vc = - 2.d0 * a * (1.d0 + 2.d0 / 3.d0 * a1 (iflag) * rs) &
          * olog - 2.d0 / 3.d0 * a * (1.d0 + a1 (iflag) * rs) * dom / &
          (om * (om + 1.d0) )
  endif
  return
end subroutine pw_dev


attributes(device) subroutine pbex_dev (rho, grho, iflag, sx, v1x, v2x)
  !---------------------------------------------------------------
  !
  ! PBE exchange (without Slater exchange):
  ! iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
  ! iflag=2  "revised' PBE: Y. Zhang et al., PRL 80, 890 (1998)
  ! iflag=3  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
  ! iflag=4  PBEQ2D: L. Chiodo et al., PRL 108, 126402 (2012)
  ! iflag=5  optB88: Klimes et al., J. Phys. Cond. Matter, 22, 022201 (2010)
  ! iflag=6  optB86b: Klimes et al., Phys. Rev. B 83, 195131 (2011)
  ! iflag=7  ev: Engel and Vosko, PRB 47, 13164 (1991)
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : pi
  implicit none
  real(dp), value :: rho, grho
  ! input: charge and squared gradient
  real(dp), device, intent(out):: sx, v1x, v2x
  ! output: energy, potential
  integer, value :: iflag
  ! local variables
  real(DP) :: kf, agrho, s1, s2, ds, dsg, exunif, fx
  ! (3*pi2*|rho|)^(1/3)
  ! |grho|
  ! |grho|/(2*kf*|rho|)
  ! s^2
  ! n*ds/dn
  ! n*ds/d(gn)
  ! exchange energy LDA part
  ! exchange energy gradient part
  real(DP) :: dxunif, dfx, f1, f2, f3, dfx1
  real(DP) :: p,  amu, ab, c, dfxdp, dfxds, upbe, uge, s, ak, aa
  ! numerical coefficients (NB: c2=(3 pi^2)^(1/3) )
  real(DP), parameter :: third = 1._DP / 3._DP, c1 = 0.75_DP / pi , &
       c2 = 3.093667726280136_DP, c5 = 4._DP * third, &
       c6 = c2*2.51984210, c7=5._DP/6._DP, c8=0.8_DP ! (3pi^2)^(1/3)*2^(4/3)
  ! parameters of the functional
  real(DP) :: k (6), mu(6), ev(6)
  !           pbe        rpbe        pbesol   pbeq2d      optB88  optB86b
  data k / 0.804_DP,   1.2450_DP,   0.804_DP , 0.804_DP,  0.0_dp,  0.0_dp/, &
       mu/ 0.21951_DP, 0.21951_DP, 0.12345679012345679_DP,             &
                                   0.12345679012345679_DP,  0.22_dp, 0.1234_dp/, &
       ev / 1.647127_DP, 0.980118_DP, 0.017399_DP, 1.523671_DP, 0.367229_DP, &
                                   0.011282_DP /  ! a and b parameters of Engel and Vosko
  agrho = sqrt (grho)
  kf = c2 * rho**third
  dsg = 0.5_DP / kf
  s1 = agrho * dsg / rho
  s2 = s1 * s1
  ds = - c5 * s1
  !
  !   Energy
  !
  if ( iflag == 4) then
     p=s1*s1
     s=s1
     ak=0.804_DP
     amu=10._DP/81._DP
     ab=0.5217_DP
     c=2._DP
     fx =  ak - ak / (1.0_dp + amu * p / ak)  + p**2 * (1 + p) &
           /(10**c + p**3) * (-1.0_dp - ak + ak / (1.0_dp + amu * p / ak) &
           + ab * p ** (-0.1d1/ 0.4D1))
  elseif ( iflag == 5) then
     ab=mu(iflag)*c7 ! mu/ab=1.2
     p=s1*c6
     c=log(p+sqrt(p*p+1)) ! asinh(p)
     dfx1=1+ab*s1*c
     fx =  mu(iflag)*s1*s1/dfx1
  elseif ( iflag == 6) then
     p=mu(iflag)*s1*s1
     fx =  p / ( 1 + p )**c8
  elseif ( iflag == 7) then
     s=s2*s2
     f1 =  1 + ev(1)*s2 + ev(2)*s + ev(3)*s*s2
     f2 =  1 + ev(4)*s2 + ev(5)*s + ev(6)*s*s2
     fx = f1 / f2 - 1
  else
     f1 = s2 * mu(iflag) / k (iflag)
     f2 = 1._DP + f1
     f3 = k (iflag) / f2
     fx = k (iflag) - f3
  end if
  exunif = - c1 * kf
  sx = exunif * fx
  !
  !   Potential
  !
  dxunif = exunif * third
  if ( iflag == 4) then
      dfxdp = dble(1 / (1 + amu * p / ak) ** 2 * amu) + dble(2 * p * (1 &
     + p) / (10 ** c + p ** 3) * (-1 - ak + ak / (1 + amu * p / ak) + ab &
      * p ** (-0.1d1 / 0.4D1))) + dble(p ** 2 / (10 ** c + p ** 3) * ( &
     -1 - ak + ak / (1 + amu * p / ak) + ab * p ** (-0.1d1 / 0.4D1))) - &
      dble(3 * p ** 4 * (1 + p) / (10 ** c + p ** 3) ** 2 * (-1 - ak + &
     ak / (1 + amu * p / ak) + ab * p ** (-0.1d1 / 0.4D1))) + dble(p ** &
      2) * dble(1 + p) / dble(10 ** c + p ** 3) * (-dble(1 / (1 + amu * &
      p / ak) ** 2 * amu) - dble(ab * p ** (-0.5d1 / 0.4D1)) / 0.4D1)

      dfxds=dfxdp*2._DP*s
      dfx=dfxds
  elseif (iflag == 5) then
     dfx=2*fx/s1-fx/dfx1*(ab*c+ab*s1/sqrt(p*p+1)*c6)
  elseif (iflag == 6) then
     dfx=2*mu(iflag)*s1*fx*(1+(1-c8)*p)/(p*(1+p))
  elseif (iflag == 7) then
    dfx  =  ev(1) + 2*ev(2)*s2 + 3*ev(3)*s  
    dfx1 =  ev(4) + 2*ev(5)*s2 + 3*ev(6)*s 
    dfx  = 2 * s1 * ( dfx - f1*dfx1/f2 ) / f2
  else
     dfx1 = f2 * f2
     dfx = 2._DP * mu(iflag) * s1 / dfx1
  end if
  v1x = sx + dxunif * fx + exunif * dfx * ds
  v2x = exunif * dfx * dsg / agrho
  sx = sx * rho
  return
end subroutine pbex_dev

attributes(device) subroutine pbec_dev (rho, grho, iflag, sc, v1c, v2c)
  !---------------------------------------------------------------
  !
  ! PBE correlation (without LDA part)
  ! iflag=1: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
  ! iflag=2: J.P.Perdew et al., PRL 100, 136406 (2008).
  ! iflag=3: L. Chiodo et al, PRL 108, 126402 (2012)  (PBEQ2D)
  !
  USE kinds, ONLY : DP
  implicit none
  integer, value :: iflag
  real(DP), value :: rho, grho
  real(DP), device, intent(out):: sc, v1c, v2c
  real(DP), parameter :: ga = 0.031091d0
  real(DP) :: be (3)
!             pbe           pbesol   pbeq2d
  data be / 0.066725d0, 0.046d0,     0.066725d0/
  real(DP), parameter :: third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0
  real(DP), parameter :: xkf = 1.919158292677513d0, xks = 1.128379167095513d0
  ! pi34=(3/4pi)^(1/3), xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
  real(DP) :: kf, ks, rs, ec, vc, t, expe, af, bf, y, xy, qy
  real(DP) :: s1, h0, dh0, ddh0, sc2D, v1c2D, v2c2D
  !
  rs = pi34 / rho**third
  call pw_dev (rs, 1, ec, vc)
  kf = xkf / rs
  ks = xks * sqrt (kf)
  t = sqrt (grho) / (2.d0 * ks * rho)
  expe = exp ( - ec / ga)
  af = be(iflag) / ga * (1.d0 / (expe-1.d0) )
  bf = expe * (vc - ec)
  y = af * t * t
  xy = (1.d0 + y) / (1.d0 + y + y * y)
  qy = y * y * (2.d0 + y) / (1.d0 + y + y * y) **2
  s1 = 1.d0 + be(iflag) / ga * t * t * xy
  h0 = ga * log (s1)
  dh0 = be(iflag) * t * t / s1 * ( - 7.d0 / 3.d0 * xy - qy * (af * bf / &
       be(iflag)-7.d0 / 3.d0) )
  ddh0 = be(iflag) / (2.d0 * ks * ks * rho) * (xy - qy) / s1
  sc = rho * h0
  v1c = h0 + dh0
  v2c = ddh0
! q2D
  !if (iflag == 3)then
  !   call cpbe2d(rho,grho,sc2D,v1c2D,v2c2D)
  !   sc=sc+sc2D
  !   v1c=v1c+v1c2D
  !   v2c=v2c+v2c2D
  !endif
  !
  return
end subroutine pbec_dev

attributes(device) subroutine ggax_dev (rho, grho, sx, v1x, v2x)
  !-----------------------------------------------------------------------
  ! Perdew-Wang GGA (PW91), exchange part:
  ! J.P. Perdew et al.,PRB 46, 6671 (1992)
  !
  USE kinds, ONLY : DP
  implicit none
  real(DP), value :: rho, grho
  real(DP), device, intent(out) :: sx, v1x, v2x
  real(DP) :: f1, f2, f3, f4, f5
  parameter (f1 = 0.19645d0, f2 = 7.7956d0, f3 = 0.2743d0, f4 = &
       0.1508d0, f5 = 0.004d0)
  real(DP) :: fp1, fp2
  parameter (fp1 = -0.019292021296426d0, fp2 = 0.161620459673995d0)
  ! fp1 = -3/(16 pi)*(3 pi^2)^(-1/3)
  ! fp2 = (1/2)(3 pi^2)**(-1/3)
  real(DP) :: rhom43, s, s2, s3, s4, exps, as, sa2b8, shm1, bs, das, &
       dbs, dls
  !
  rhom43 = rho** ( - 4.d0 / 3.d0)
  s = fp2 * sqrt (grho) * rhom43
  s2 = s * s
  s3 = s2 * s
  s4 = s2 * s2
  exps = f4 * exp ( - 100.d0 * s2)
  as = f3 - exps - f5 * s2
  sa2b8 = sqrt (1.0d0 + f2 * f2 * s2)
  shm1 = log (f2 * s + sa2b8)
  bs = 1.d0 + f1 * s * shm1 + f5 * s4
  das = (200.d0 * exps - 2.d0 * f5) * s
  dbs = f1 * (shm1 + f2 * s / sa2b8) + 4.d0 * f5 * s3
  dls = (das / as - dbs / bs)
  sx = fp1 * grho * rhom43 * as / bs
  v1x = - 4.d0 / 3.d0 * sx / rho * (1.d0 + s * dls)
  v2x = fp1 * rhom43 * as / bs * (2.d0 + s * dls)
  !
  return
end subroutine ggax_dev

attributes(device) subroutine ggac_dev (rho, grho, sc, v1c, v2c)
  !-----------------------------------------------------------------------
  ! Perdew-Wang GGA (PW91) correlation part
  !
  USE kinds, ONLY : DP
  implicit none
  real(DP), value :: rho, grho
  real(DP), device, intent(out) :: sc, v1c, v2c
  real(DP) :: al, pa, pb, pc, pd, cx, cxc0, cc0
  parameter (al = 0.09d0, pa = 0.023266d0, pb = 7.389d-6, pc = &
       8.723d0, pd = 0.472d0)
  parameter (cx = -0.001667d0, cxc0 = 0.002568d0, cc0 = - cx + cxc0)
  real(DP) :: third, pi34, nu, be, xkf, xks
  parameter (third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0)
  parameter (nu = 15.755920349483144d0, be = nu * cc0)
  parameter (xkf = 1.919158292677513d0, xks = 1.128379167095513d0)
  ! pi34=(3/4pi)^(1/3),  nu=(16/pi)*(3 pi^2)^(1/3)
  ! xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
  real(DP) :: kf, ks, rs, rs2, rs3, ec, vc, t, expe, af, bf, y, xy, &
       qy, s1
  real(DP) :: h0, dh0, ddh0, ee, cn, dcn, cna, dcna, cnb, dcnb, h1, &
       dh1, ddh1
  !
  rs = pi34 / rho**third
  rs2 = rs * rs
  rs3 = rs * rs2
  call pw_dev (rs, 1, ec, vc)
  kf = xkf / rs
  ks = xks * sqrt (kf)
  t = sqrt (grho) / (2.d0 * ks * rho)
  expe = exp ( - 2.d0 * al * ec / (be * be) )
  af = 2.d0 * al / be * (1.d0 / (expe-1.d0) )
  bf = expe * (vc - ec)
  y = af * t * t
  xy = (1.d0 + y) / (1.d0 + y + y * y)
  qy = y * y * (2.d0 + y) / (1.d0 + y + y * y) **2
  s1 = 1.d0 + 2.d0 * al / be * t * t * xy
  h0 = be * be / (2.d0 * al) * log (s1)
  dh0 = be * t * t / s1 * ( - 7.d0 / 3.d0 * xy - qy * (af * bf / &
       be-7.d0 / 3.d0) )
  ddh0 = be / (2.d0 * ks * ks * rho) * (xy - qy) / s1
  ee = - 100.d0 * (ks / kf * t) **2
  cna = cxc0 + pa * rs + pb * rs2
  dcna = pa * rs + 2.d0 * pb * rs2
  cnb = 1.d0 + pc * rs + pd * rs2 + 1.d4 * pb * rs3
  dcnb = pc * rs + 2.d0 * pd * rs2 + 3.d4 * pb * rs3
  cn = cna / cnb - cx
  dcn = dcna / cnb - cna * dcnb / (cnb * cnb)
  h1 = nu * (cn - cc0 - 3.d0 / 7.d0 * cx) * t * t * exp (ee)
  dh1 = - third * (h1 * (7.d0 + 8.d0 * ee) + nu * t * t * exp (ee) &
       * dcn)
  ddh1 = 2.d0 * h1 * (1.d0 + ee) * rho / grho
  sc = rho * (h0 + h1)
  v1c = h0 + h1 + dh0 + dh1
  v2c = ddh0 + ddh1
  !
  return
end subroutine ggac_dev

!-----------------------------------------------------------------------
! subroutines ported from lsda_functionals.f90
!-----------------------------------------------------------------------

attributes(device) subroutine slater_spin_dev (rho, zeta, ex, vxup, vxdw)
  !-----------------------------------------------------------------------
  !     Slater exchange with alpha=2/3, spin-polarized case
  !
  USE kinds, ONLY : DP
  implicit none
  real(DP), value :: rho, zeta
  real(DP), device, intent(out) :: ex, vxup, vxdw
  real(DP) :: f, alpha, third, p43
  parameter (f = - 1.10783814957303361d0, alpha = 2.0d0 / 3.0d0)
  ! f = -9/8*(3/pi)^(1/3)
  parameter (third = 1.d0 / 3.d0, p43 = 4.d0 / 3.d0)
  real(DP) :: exup, exdw, rho13
  !
  rho13 = ( (1.d0 + zeta) * rho) **third
  exup = f * alpha * rho13
  vxup = p43 * f * alpha * rho13
  rho13 = ( (1.d0 - zeta) * rho) **third
  exdw = f * alpha * rho13
  vxdw = p43 * f * alpha * rho13
  ex = 0.5d0 * ( (1.d0 + zeta) * exup + (1.d0 - zeta) * exdw)
  !
  return
end subroutine slater_spin_dev

attributes(device) subroutine pw_spin_dev (rs, zeta, ec, vcup, vcdw)
  !-----------------------------------------------------------------------
  !     J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)
  !
  USE kinds, ONLY : DP
  implicit none
  real(DP), value :: rs, zeta
  real(DP), device, intent(out) :: ec, vcup, vcdw
  ! xc parameters, unpolarised
  real(DP) :: a, a1, b1, b2, b3, b4, c0, c1, c2, c3, d0, d1
  parameter (a = 0.031091d0, a1 = 0.21370d0, b1 = 7.5957d0, b2 = &
       3.5876d0, b3 = 1.6382d0, b4 = 0.49294d0, c0 = a, c1 = 0.046644d0, &
       c2 = 0.00664d0, c3 = 0.01043d0, d0 = 0.4335d0, d1 = 1.4408d0)
  ! xc parameters, polarised
  real(DP) :: ap, a1p, b1p, b2p, b3p, b4p, c0p, c1p, c2p, c3p, d0p, &
       d1p
  parameter (ap = 0.015545d0, a1p = 0.20548d0, b1p = 14.1189d0, b2p &
       = 6.1977d0, b3p = 3.3662d0, b4p = 0.62517d0, c0p = ap, c1p = &
       0.025599d0, c2p = 0.00319d0, c3p = 0.00384d0, d0p = 0.3287d0, d1p &
       = 1.7697d0)
  ! xc parameters, antiferro
  real(DP) :: aa, a1a, b1a, b2a, b3a, b4a, c0a, c1a, c2a, c3a, d0a, &
       d1a
  parameter (aa = 0.016887d0, a1a = 0.11125d0, b1a = 10.357d0, b2a = &
       3.6231d0, b3a = 0.88026d0, b4a = 0.49671d0, c0a = aa, c1a = &
       0.035475d0, c2a = 0.00188d0, c3a = 0.00521d0, d0a = 0.2240d0, d1a &
       = 0.3969d0)
  real(DP) :: fz0
  parameter (fz0 = 1.709921d0)
  real(DP) :: rs12, rs32, rs2, zeta2, zeta3, zeta4, fz, dfz
  real(DP) :: om, dom, olog, epwc, vpwc
  real(DP) :: omp, domp, ologp, epwcp, vpwcp
  real(DP) :: oma, doma, ologa, alpha, vpwca
  !
  !     if(rs.lt.0.5d0) then
  ! high density formula (not implemented)
  !
  !     else if(rs.gt.100.d0) then
  ! low density formula  (not implemented)
  !
  !     else
  ! interpolation formula
  zeta2 = zeta * zeta
  zeta3 = zeta2 * zeta
  zeta4 = zeta3 * zeta
  rs12 = sqrt (rs)
  rs32 = rs * rs12
  rs2 = rs**2
  ! unpolarised
  om = 2.d0 * a * (b1 * rs12 + b2 * rs + b3 * rs32 + b4 * rs2)
  dom = 2.d0 * a * (0.5d0 * b1 * rs12 + b2 * rs + 1.5d0 * b3 * rs32 &
       + 2.d0 * b4 * rs2)
  olog = log (1.d0 + 1.0d0 / om)
  epwc = - 2.d0 * a * (1.d0 + a1 * rs) * olog
  vpwc = - 2.d0 * a * (1.d0 + 2.d0 / 3.d0 * a1 * rs) * olog - 2.d0 / &
       3.d0 * a * (1.d0 + a1 * rs) * dom / (om * (om + 1.d0) )
  ! polarized
  omp = 2.d0 * ap * (b1p * rs12 + b2p * rs + b3p * rs32 + b4p * rs2)
  domp = 2.d0 * ap * (0.5d0 * b1p * rs12 + b2p * rs + 1.5d0 * b3p * &
       rs32 + 2.d0 * b4p * rs2)
  ologp = log (1.d0 + 1.0d0 / omp)
  epwcp = - 2.d0 * ap * (1.d0 + a1p * rs) * ologp
  vpwcp = - 2.d0 * ap * (1.d0 + 2.d0 / 3.d0 * a1p * rs) * ologp - &
       2.d0 / 3.d0 * ap * (1.d0 + a1p * rs) * domp / (omp * (omp + 1.d0) &
       )
  ! antiferro
  oma = 2.d0 * aa * (b1a * rs12 + b2a * rs + b3a * rs32 + b4a * rs2)
  doma = 2.d0 * aa * (0.5d0 * b1a * rs12 + b2a * rs + 1.5d0 * b3a * &
       rs32 + 2.d0 * b4a * rs2)
  ologa = log (1.d0 + 1.0d0 / oma)
  alpha = 2.d0 * aa * (1.d0 + a1a * rs) * ologa
  vpwca = + 2.d0 * aa * (1.d0 + 2.d0 / 3.d0 * a1a * rs) * ologa + &
       2.d0 / 3.d0 * aa * (1.d0 + a1a * rs) * doma / (oma * (oma + 1.d0) &
       )
  !
  fz = ( (1.d0 + zeta) ** (4.d0 / 3.d0) + (1.d0 - zeta) ** (4.d0 / &
       3.d0) - 2.d0) / (2.d0** (4.d0 / 3.d0) - 2.d0)
  dfz = ( (1.d0 + zeta) ** (1.d0 / 3.d0) - (1.d0 - zeta) ** (1.d0 / &
       3.d0) ) * 4.d0 / (3.d0 * (2.d0** (4.d0 / 3.d0) - 2.d0) )
  !
  ec = epwc + alpha * fz * (1.d0 - zeta4) / fz0 + (epwcp - epwc) &
       * fz * zeta4
  !
  vcup = vpwc + vpwca * fz * (1.d0 - zeta4) / fz0 + (vpwcp - vpwc) &
       * fz * zeta4 + (alpha / fz0 * (dfz * (1.d0 - zeta4) - 4.d0 * fz * &
       zeta3) + (epwcp - epwc) * (dfz * zeta4 + 4.d0 * fz * zeta3) ) &
       * (1.d0 - zeta)

  vcdw = vpwc + vpwca * fz * (1.d0 - zeta4) / fz0 + (vpwcp - vpwc) &
       * fz * zeta4 - (alpha / fz0 * (dfz * (1.d0 - zeta4) - 4.d0 * fz * &
       zeta3) + (epwcp - epwc) * (dfz * zeta4 + 4.d0 * fz * zeta3) ) &
       * (1.d0 + zeta)
  !      endif
  !
  return
end subroutine pw_spin_dev


attributes(device) subroutine pbec_spin_dev (rho, zeta, grho, iflag, sc, v1cup, v1cdw, v2c)
  !---------------------------------------------------------------
  !
  ! PBE correlation (without LDA part) - spin-polarized
  ! iflag = 1: J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996).
  ! iflag = 2: J.P.Perdew et al., PRL 100, 136406 (2008)
  !
  USE kinds, ONLY : DP
  implicit none
  integer, value :: iflag
  real(DP), value :: rho, zeta, grho
  real(DP), device, intent(out) :: sc, v1cup, v1cdw, v2c
  real(DP) :: ga, be(2)
  parameter (ga = 0.031091d0)
  data be / 0.066725d0 ,  0.046d0 /
  real(DP) :: third, pi34, xkf, xks
  parameter (third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0)
  parameter (xkf = 1.919158292677513d0, xks = 1.128379167095513d0)
  ! pi34=(3/4pi)^(1/3), xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
  real(DP) :: kf, ks, rs, ec, vcup, vcdw, t, expe, af, y, xy, qy, &
       s1, h0, ddh0
  real(DP) :: fz, fz2, fz3, fz4, dfz, bfup, bfdw, dh0up, dh0dw, &
       dh0zup, dh0zdw
  !
  rs = pi34 / rho**third
  call pw_spin_dev (rs, zeta, ec, vcup, vcdw)
  kf = xkf / rs
  ks = xks * sqrt (kf)
  fz = 0.5d0 * ( (1.d0 + zeta) ** (2.d0 / 3.d0) + (1.d0 - zeta) ** ( &
       2.d0 / 3.d0) )
  fz2 = fz * fz
  fz3 = fz2 * fz
  fz4 = fz3 * fz
  dfz = ( (1.d0 + zeta) ** ( - 1.d0 / 3.d0) - (1.d0 - zeta) ** ( - &
       1.d0 / 3.d0) ) / 3.d0
  t = sqrt (grho) / (2.d0 * fz * ks * rho)
  expe = exp ( - ec / (fz3 * ga) )
  af = be(iflag) / ga * (1.d0 / (expe-1.d0) )
  bfup = expe * (vcup - ec) / fz3
  bfdw = expe * (vcdw - ec) / fz3
  y = af * t * t
  xy = (1.d0 + y) / (1.d0 + y + y * y)
  qy = y * y * (2.d0 + y) / (1.d0 + y + y * y) **2
  s1 = 1.d0 + be(iflag) / ga * t * t * xy
  h0 = fz3 * ga * log (s1)
  dh0up = be(iflag) * t * t * fz3 / s1 * ( - 7.d0 / 3.d0 * xy - qy * &
       (af * bfup / be(iflag)-7.d0 / 3.d0) )
  dh0dw = be(iflag) * t * t * fz3 / s1 * ( - 7.d0 / 3.d0 * xy - qy * &
       (af * bfdw / be(iflag)-7.d0 / 3.d0) )
  dh0zup = (3.d0 * h0 / fz - be(iflag) * t * t * fz2 / s1 * (2.d0 * xy - &
  qy * (3.d0 * af * expe * ec / fz3 / be(iflag)+2.d0) ) ) * dfz * (1.d0 - zeta)
  dh0zdw = - (3.d0 * h0 / fz - be(iflag) * t * t * fz2 / s1 * (2.d0 * xy - &
  qy * (3.d0 * af * expe * ec / fz3 / be(iflag)+2.d0) ) ) * dfz * (1.d0 + zeta)

  ddh0 = be(iflag) * fz / (2.d0 * ks * ks * rho) * (xy - qy) / s1
  sc = rho * h0
  v1cup = h0 + dh0up + dh0zup
  v1cdw = h0 + dh0dw + dh0zdw
  v2c = ddh0
  return
end subroutine pbec_spin_dev
!
attributes(device) subroutine ggac_spin_dev (rho, zeta, grho, sc, v1cup, v1cdw, v2c)
  !-----------------------------------------------------------------------
  ! Perdew-Wang GGA (PW91) correlation part - spin-polarized
  !
  USE kinds, ONLY : DP
  implicit none
  real(DP), value :: rho, zeta, grho
  real(DP), device, intent(out) ::  sc, v1cup, v1cdw, v2c
  real(DP) :: al, pa, pb, pc, pd, cx, cxc0, cc0
  parameter (al = 0.09d0, pa = 0.023266d0, pb = 7.389d-6, pc = &
       8.723d0, pd = 0.472d0)
  parameter (cx = - 0.001667d0, cxc0 = 0.002568d0, cc0 = - cx + &
       cxc0)
  real(DP) :: third, pi34, nu, be, xkf, xks
  parameter (third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0)
  parameter (nu = 15.755920349483144d0, be = nu * cc0)
  parameter (xkf = 1.919158292677513d0, xks = 1.128379167095513d0)
  ! pi34=(3/4pi)^(1/3),  nu=(16/pi)*(3 pi^2)^(1/3)
  ! xkf=(9 pi/4)^(1/3), xks= sqrt(4/pi)
  real(DP) :: kf, ks, rs, rs2, rs3, ec, vcup, vcdw, t, expe, af, y, &
       xy, qy, s1, h0, ddh0, ee, cn, dcn, cna, dcna, cnb, dcnb, h1, dh1, &
       ddh1, fz, fz2, fz3, fz4, dfz, bfup, bfdw, dh0up, dh0dw, dh0zup, &
       dh0zdw, dh1zup, dh1zdw
  !
  rs = pi34 / rho**third
  rs2 = rs * rs
  rs3 = rs * rs2
  call pw_spin_dev (rs, zeta, ec, vcup, vcdw)
  kf = xkf / rs
  ks = xks * sqrt (kf)
  fz = 0.5d0 * ( (1.d0 + zeta) ** (2.d0 / 3.d0) + (1.d0 - zeta) ** ( &
       2.d0 / 3.d0) )
  fz2 = fz * fz
  fz3 = fz2 * fz
  fz4 = fz3 * fz
  dfz = ( (1.d0 + zeta) ** ( - 1.d0 / 3.d0) - (1.d0 - zeta) ** ( - &
       1.d0 / 3.d0) ) / 3.d0
  t = sqrt (grho) / (2.d0 * fz * ks * rho)
  expe = exp ( - 2.d0 * al * ec / (fz3 * be * be) )
  af = 2.d0 * al / be * (1.d0 / (expe-1.d0) )
  bfup = expe * (vcup - ec) / fz3
  bfdw = expe * (vcdw - ec) / fz3
  y = af * t * t
  xy = (1.d0 + y) / (1.d0 + y + y * y)
  qy = y * y * (2.d0 + y) / (1.d0 + y + y * y) **2
  s1 = 1.d0 + 2.d0 * al / be * t * t * xy
  h0 = fz3 * be * be / (2.d0 * al) * log (s1)
  dh0up = be * t * t * fz3 / s1 * ( - 7.d0 / 3.d0 * xy - qy * &
       (af * bfup / be-7.d0 / 3.d0) )
  dh0dw = be * t * t * fz3 / s1 * ( - 7.d0 / 3.d0 * xy - qy * &
       (af * bfdw / be-7.d0 / 3.d0) )
  dh0zup = (3.d0 * h0 / fz - be * t * t * fz2 / s1 * (2.d0 * xy - &
       qy * (3.d0 * af * expe * ec / fz3 / be+2.d0) ) ) * dfz * (1.d0 - &
       zeta)
  dh0zdw = - (3.d0 * h0 / fz - be * t * t * fz3 / s1 * (2.d0 * xy - &
       qy * (3.d0 * af * expe * ec / fz3 / be+2.d0) ) ) * dfz * (1.d0 + &
       zeta)
  ddh0 = be * fz / (2.d0 * ks * ks * rho) * (xy - qy) / s1
  ee = - 100.d0 * fz4 * (ks / kf * t) **2
  cna = cxc0 + pa * rs + pb * rs2
  dcna = pa * rs + 2.d0 * pb * rs2
  cnb = 1.d0 + pc * rs + pd * rs2 + 1.d4 * pb * rs3
  dcnb = pc * rs + 2.d0 * pd * rs2 + 3.d4 * pb * rs3
  cn = cna / cnb - cx
  dcn = dcna / cnb - cna * dcnb / (cnb * cnb)
  h1 = nu * (cn - cc0 - 3.d0 / 7.d0 * cx) * fz3 * t * t * exp (ee)
  dh1 = - third * (h1 * (7.d0 + 8.d0 * ee) + fz3 * nu * t * t * exp &
       (ee) * dcn)
  ddh1 = 2.d0 * h1 * (1.d0 + ee) * rho / grho
  dh1zup = (1.d0 - zeta) * dfz * h1 * (1.d0 + 2.d0 * ee / fz)
  dh1zdw = - (1.d0 + zeta) * dfz * h1 * (1.d0 + 2.d0 * ee / fz)
  sc = rho * (h0 + h1)
  v1cup = h0 + h1 + dh0up + dh1 + dh0zup + dh1zup
  v1cdw = h0 + h1 + dh0up + dh1 + dh0zdw + dh1zdw
  v2c = ddh0 + ddh1
  return
end subroutine ggac_spin_dev
!
!
end module
#endif
