!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine stres_har (sigmahar)
  !----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : e2, fpi
  USE cell_base, ONLY: omega, tpiba2
  USE ener,      ONLY: ehart
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft
  USE gvect,     ONLY: ngm, gstart, nl, g, gg
  USE lsda_mod,  ONLY: nspin
  USE scf,       ONLY: rho
  USE control_flags,        ONLY: gamma_only
  USE wavefunctions_module, ONLY : psic
  USE mp_bands,  ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum

  implicit none
  !
  real(DP) :: sigmahar (3, 3), shart, g2
  real(DP), parameter :: eps = 1.d-8
  integer :: is, ig, l, m, nspin0

  sigmahar(:,:) = 0.d0
  psic (:) = (0.d0, 0.d0)
  nspin0=nspin
  if (nspin==4) nspin0=1
  do is = 1, nspin0
     call daxpy (dfftp%nnr, 1.d0, rho%of_r (1, is), 1, psic, 2)
  enddo

  CALL fwfft ('Dense', psic, dfftp)
  ! psic contains now the charge density in G space
  ! the  G=0 component is not computed
  do ig = gstart, ngm
     g2 = gg (ig) * tpiba2
     shart = psic (nl (ig) ) * CONJG(psic (nl (ig) ) ) / g2
     do l = 1, 3
        do m = 1, l
           sigmahar (l, m) = sigmahar (l, m) + shart * tpiba2 * 2 * &
                g (l, ig) * g (m, ig) / g2
        enddo
     enddo
  enddo
  !
  call mp_sum(  sigmahar, intra_bgrp_comm )
  !
  if (gamma_only) then
     sigmahar(:,:) =       fpi * e2 * sigmahar(:,:)
  else
     sigmahar(:,:) = 0.5d0 * fpi * e2 * sigmahar(:,:)
  end if
  do l = 1, 3
     sigmahar (l, l) = sigmahar (l, l) - ehart / omega
  enddo
  do l = 1, 3
     do m = 1, l - 1
        sigmahar (m, l) = sigmahar (l, m)
     enddo
  enddo

  sigmahar(:,:) = -sigmahar(:,:)

  return
end subroutine stres_har

#ifdef USE_CUDA
subroutine stres_har_gpu (sigmahar)
  !----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : e2, fpi
  USE cell_base, ONLY: omega, tpiba2
  USE ener,      ONLY: ehart
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft
  USE gvect,     ONLY: ngm, gstart, nl=>nl_d, g=>g_d, gg=>gg_d
  USE lsda_mod,  ONLY: nspin
  USE scf,       ONLY: rho
  USE control_flags,        ONLY: gamma_only
  USE wavefunctions_module, ONLY : psic=>psic_d
  USE mp_bands,  ONLY: intra_bgrp_comm
  USE mp,        ONLY: mp_sum
  USE cudafor

  implicit none
  !
  real(DP) :: sigmahar (3, 3), shart, g2
  real(DP), parameter :: eps = 1.d-8
  integer :: is, ig, l, m, nspin0
  real(DP), pointer, device :: rho_of_r_d(:,:)
  real(DP):: s11,s21,s22,s31,s32,s33

  sigmahar(:,:) = 0.d0
  psic (:) = (0.d0, 0.d0)
  nspin0=nspin
  if (nspin==4) nspin0=1
   rho%of_r_d=rho%of_r
   rho_of_r_d =>rho%of_r_d

 !$cuf kernel do(1) <<<*,*>>>
  do l=1,dfftp%nnr
   do is = 1, nspin0
    psic(l) = psic(l) + rho_of_r_d (l, is)
   enddo
  enddo

  CALL fwfft ('Dense', psic, dfftp)
  ! psic contains now the charge density in G space
  ! the  G=0 component is not computed
  s11=0.d0
  s21=0.d0
  s22=0.d0
  s31=0.d0
  s32=0.d0
  s33=0.d0
 !$cuf kernel do(1) <<<*,*>>>
  do ig = gstart, ngm
     g2 = gg (ig) * tpiba2
     shart = 2.d0 * tpiba2 *psic (nl (ig) ) * CONJG(psic (nl (ig) ) ) / (g2*g2)
      s11 = s11+ shart *  g(1,ig) * g(1,ig) 
      s21 = s21+ shart *  g(2,ig) * g(1,ig) 
      s22 = s22+ shart *  g(2,ig) * g(2,ig) 
      s31 = s31+ shart *  g(3,ig) * g(1,ig) 
      s32 = s32+ shart *  g(3,ig) * g(2,ig) 
      s33 = s33+ shart *  g(3,ig) * g(3,ig) 
  enddo
  sigmahar(1,1)=s11
  sigmahar(2,1)=s21
  sigmahar(2,2)=s22
  sigmahar(3,1)=s31
  sigmahar(3,2)=s32
  sigmahar(3,3)=s33
  !
  call mp_sum(  sigmahar, intra_bgrp_comm )
  !
  if (gamma_only) then
     sigmahar(:,:) =       fpi * e2 * sigmahar(:,:)
  else
     sigmahar(:,:) = 0.5d0 * fpi * e2 * sigmahar(:,:)
  end if
  do l = 1, 3
     sigmahar (l, l) = sigmahar (l, l) - ehart / omega
  enddo
  do l = 1, 3
     do m = 1, l - 1
        sigmahar (m, l) = sigmahar (l, m)
     enddo
  enddo

  sigmahar(:,:) = -sigmahar(:,:)

  return
end subroutine stres_har_gpu
#endif
