!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine stres_knl (sigmanlc, sigmakin)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY: DP
  USE constants,            ONLY: pi, e2
  USE cell_base,            ONLY: omega, alat, at, bg, tpiba
  USE gvect,                ONLY: g
  USE gvecw,                ONLY: qcutz, ecfixed, q2sigma
  USE klist,                ONLY: nks, xk, ngk, igk_k
  USE io_files,             ONLY: iunwfc, nwordwfc
  USE buffers,              ONLY: get_buffer
  USE symme,                ONLY: symmatrix
  USE wvfct,                ONLY: npwx, nbnd, wg
  USE control_flags,        ONLY: gamma_only
  USE noncollin_module,     ONLY: noncolin, npol
  USE wavefunctions_module, ONLY: evc
#ifdef USE_CUDA
  USE wavefunctions_module, ONLY: evc_d
#endif
  USE mp_pools,             ONLY: inter_pool_comm
  USE mp_bands,             ONLY: intra_bgrp_comm
  USE mp,                   ONLY: mp_sum
  implicit none
  real(DP) :: sigmanlc (3, 3), sigmakin (3, 3)
  real(DP), allocatable :: gk (:,:), kfac (:)
  real(DP) :: twobysqrtpi, gk2, arg
  integer :: npw, ik, l, m, i, ibnd, is

  allocate (gk(  3, npwx))    
  allocate (kfac(   npwx))    

  sigmanlc(:,:) =0.d0
  sigmakin(:,:) =0.d0
  twobysqrtpi = 2.d0 / sqrt (pi)

  kfac(:) = 1.d0

#ifdef USE_CUDA
  evc = evc_d
#endif

  do ik = 1, nks
     if (nks > 1) &
        call get_buffer (evc, nwordwfc, iunwfc, ik)
     npw = ngk(ik)
     do i = 1, npw
        gk (1, i) = (xk (1, ik) + g (1, igk_k(i,ik) ) ) * tpiba
        gk (2, i) = (xk (2, ik) + g (2, igk_k(i,ik) ) ) * tpiba
        gk (3, i) = (xk (3, ik) + g (3, igk_k(i,ik) ) ) * tpiba
        if (qcutz.gt.0.d0) then
           gk2 = gk (1, i) **2 + gk (2, i) **2 + gk (3, i) **2
           arg = ( (gk2 - ecfixed) / q2sigma) **2
           kfac (i) = 1.d0 + qcutz / q2sigma * twobysqrtpi * exp ( - arg)
        endif
     enddo
     !
     !   kinetic contribution
     !
     do l = 1, 3
        do m = 1, l
           do ibnd = 1, nbnd
              do i = 1, npw
                 if (noncolin) then
                    sigmakin (l, m) = sigmakin (l, m) + wg (ibnd, ik) * &
                     gk (l, i) * gk (m, i) * kfac (i) * &
                     ( DBLE (CONJG(evc(i     ,ibnd))*evc(i     ,ibnd)) + &
                       DBLE (CONJG(evc(i+npwx,ibnd))*evc(i+npwx,ibnd)))
                 else
                    sigmakin (l, m) = sigmakin (l, m) + wg (ibnd, ik) * &
                        gk (l, i) * gk (m, i) * kfac (i) * &
                          DBLE (CONJG(evc (i, ibnd) ) * evc (i, ibnd) )
                 end if
              enddo
           enddo
        enddo
     enddo
     !
     !  contribution from the  nonlocal part
     !
     call stres_us (ik, gk, sigmanlc)

  enddo
  !
  ! add the US term from augmentation charge derivatives
  !
  call addusstres (sigmanlc)
  !
  call mp_sum( sigmakin, intra_bgrp_comm )
  call mp_sum( sigmanlc, intra_bgrp_comm )
  call mp_sum( sigmakin, inter_pool_comm )
  call mp_sum( sigmanlc, inter_pool_comm )
  !
  do l = 1, 3
     do m = 1, l - 1
        sigmanlc (m, l) = sigmanlc (l, m)
        sigmakin (m, l) = sigmakin (l, m)
     enddo
  enddo
  !
  if (gamma_only) then
     sigmakin(:,:) = 2.d0 * e2 / omega * sigmakin(:,:)
  else
     sigmakin(:,:) = e2 / omega * sigmakin(:,:)
  end if
  sigmanlc(:,:) = -1.d0 / omega * sigmanlc(:,:)
  !
  ! symmetrize stress
  !
  call symmatrix ( sigmakin )
  call symmatrix ( sigmanlc )

  deallocate(kfac)
  deallocate(gk)
  return
end subroutine stres_knl

#ifdef USE_CUDA
subroutine stres_knl_gpu (sigmanlc, sigmakin)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY: DP
  USE constants,            ONLY: pi, e2
  USE cell_base,            ONLY: omega, alat, at, bg, tpiba
  USE gvect,                ONLY: g=>g_d
  USE gvecw,                ONLY: qcutz, ecfixed, q2sigma
  USE klist,                ONLY: nks, xk, ngk, igk_k=>igk_k_d
  USE io_files,             ONLY: iunwfc, nwordwfc
  USE buffers,              ONLY: get_buffer
  USE symme,                ONLY: symmatrix
  USE wvfct,                ONLY: npwx, nbnd, wg=>wg_d
  USE control_flags,        ONLY: gamma_only
  USE noncollin_module,     ONLY: noncolin, npol
  USE wavefunctions_module, ONLY: evc
#ifdef USE_CUDA
  USE wavefunctions_module, ONLY: evc_d
#endif
  USE mp_pools,             ONLY: inter_pool_comm
  USE mp_bands,             ONLY: intra_bgrp_comm
  USE mp,                   ONLY: mp_sum
  implicit none
  real(DP) :: sigmanlc (3, 3), sigmakin (3, 3)
  real(DP), allocatable,device :: gk (:,:), kfac (:)
  real(DP) :: twobysqrtpi, gk2, arg
  real(DP) :: xk1,xk2,xk3,tmpf
  real(DP) :: s11,s21,s22,s31,s32,s33
  integer :: npw, ik, l, m, i, ibnd, is

  allocate (gk(  3, npwx))    
  allocate (kfac(   npwx))    

  sigmanlc(:,:) =0.d0
  sigmakin(:,:) =0.d0
  twobysqrtpi = 2.d0 / sqrt (pi)

  kfac(:) = 1.d0

#ifdef USE_CUDA
  evc = evc_d
#endif

  do ik = 1, nks
     if (nks > 1) then
        call get_buffer (evc, nwordwfc, iunwfc, ik)
       evc_d=evc ! this copy can be async
     endif
     
     npw = ngk(ik)
     xk1=xk (1, ik)
     xk2=xk (2, ik)
     xk3=xk (3, ik)
   !$cuf kernel do (1) <<<*, *>>>
     do i = 1, npw
        gk (1, i) = (xk1 + g (1, igk_k(i,ik) ) ) * tpiba
        gk (2, i) = (xk2 + g (2, igk_k(i,ik) ) ) * tpiba
        gk (3, i) = (xk3 + g (3, igk_k(i,ik) ) ) * tpiba
        if (qcutz.gt.0.d0) then
           gk2 = gk(1,i)*gk(1,i) + gk(2,i)*gk(2,i) + gk(3,i) *gk(3,i)
           arg = ( (gk2 - ecfixed) / q2sigma) **2
           kfac (i) = 1.d0 + qcutz / q2sigma * twobysqrtpi * exp ( - arg)
        endif
     enddo
     !
     !   kinetic contribution
     !
          s11=0.d0
          s21=0.d0
          s22=0.d0
          s31=0.d0
          s32=0.d0
          s33=0.d0
           do ibnd = 1, nbnd
             !$cuf kernel do (1) <<<*, *>>>
              do i = 1, npw
                 if (noncolin) then
                    tmpf = wg (ibnd, ik) * kfac (i) * &
                     ( DBLE (CONJG(evc_d(i     ,ibnd))*evc_d(i     ,ibnd)) + &
                       DBLE (CONJG(evc_d(i+npwx,ibnd))*evc_d(i+npwx,ibnd)))
                 else
                    tmpf =  wg (ibnd, ik) * kfac (i) * &
                          DBLE (CONJG(evc_d (i, ibnd) ) * evc_d (i, ibnd) )
                 end if
                    s11= s11+ tmpf * gk(1, i) * gk(1,i) 
                    s21= s21+ tmpf * gk(2, i) * gk(1,i) 
                    s22= s22+ tmpf * gk(2, i) * gk(2,i) 
                    s31= s31+ tmpf * gk(3, i) * gk(1,i) 
                    s32= s32+ tmpf * gk(3, i) * gk(2,i) 
                    s33= s33+ tmpf * gk(3, i) * gk(3,i) 
              enddo
            enddo
            sigmakin (1, 1) = sigmakin (1, 1) + s11
            sigmakin (2, 1) = sigmakin (2, 1) + s21
            sigmakin (2, 2) = sigmakin (2, 2) + s22
            sigmakin (3, 1) = sigmakin (3, 1) + s31
            sigmakin (3, 2) = sigmakin (3, 2) + s32
            sigmakin (3, 3) = sigmakin (3, 3) + s33

     !
     !  contribution from the  nonlocal part
     !
     !call stres_us (ik, gk, sigmanlc)

  enddo
  !
  ! add the US term from augmentation charge derivatives
  !
  call addusstres_gpu (sigmanlc)
  !
  call mp_sum( sigmakin, intra_bgrp_comm )
  call mp_sum( sigmanlc, intra_bgrp_comm )
  call mp_sum( sigmakin, inter_pool_comm )
  call mp_sum( sigmanlc, inter_pool_comm )
  !
  do l = 1, 3
     do m = 1, l - 1
        sigmanlc (m, l) = sigmanlc (l, m)
        sigmakin (m, l) = sigmakin (l, m)
     enddo
  enddo
  !
  if (gamma_only) then
     sigmakin(:,:) = 2.d0 * e2 / omega * sigmakin(:,:)
  else
     sigmakin(:,:) = e2 / omega * sigmakin(:,:)
  end if
  sigmanlc(:,:) = -1.d0 / omega * sigmanlc(:,:)
  !
  ! symmetrize stress
  !
  call symmatrix ( sigmakin )
  call symmatrix ( sigmanlc )

  deallocate(kfac)
  deallocate(gk)
  return
end subroutine stres_knl_gpu

#endif
