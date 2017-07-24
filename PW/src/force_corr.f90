!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#ifdef USE_CUDA
module compute_rhocgnt_gpu_m
contains
  attributes(global) subroutine  compute_rhocgnt_gpu(n, mesh, tpiba, gl, r, rho_at, rab, rhocgnt)
  use cudafor
  use kinds, ONLY: DP
  implicit none
 
  integer, value :: n, mesh
  real(DP), value :: tpiba
  real(DP), device, intent(in) :: gl(n), r(mesh), rho_at(mesh), rab(mesh)
  real(DP), device, intent(out) :: rhocgnt(n)
 
  integer :: tx, ty, igl, ir
  real(DP):: mysum, val, gx, x
 
  tx = threadIdx%x
  ty = threadIdx%y
 
  igl = (blockIdx%x - 1) * blockDim%y + ty
 
  if (igl > n) return
 
  gx = sqrt(gl(igl)) * tpiba
  mysum = 0.d0
 
  do ir = tx, mesh, blockDim%x 
    val = rho_at(ir) * rab(ir)

    if (r(ir)  .ge. 1.0d-8) then
      val = val * sin (gx*r(ir)) / (gx *r(ir)) 
    end if

    if (ir == 1 .or. ir == mesh) then
      mysum = mysum + val
    else if (mod(ir,2)) then
      mysum = mysum + 4.d0*val
    else
      mysum = mysum + 2.d0*val
    endif
  end do 
 
 ! Reduce by warp
       val = __shfl_down(mysum,1)
       mysum = mysum + val
       val = __shfl_down(mysum,2)
       mysum = mysum + val
       val = __shfl_down(mysum,4)
       mysum = mysum + val
       val = __shfl_down(mysum,8)
       mysum = mysum + val
       val = __shfl_down(mysum,16)
       mysum = mysum + val
 
       if (tx == 1) then
         rhocgnt(igl) =  mysum / 3.d0 
       endif 
 
  end subroutine compute_rhocgnt_gpu
end module compute_rhocgnt_gpu_m

#endif
!-----------------------------------------------------------------------
subroutine force_corr (forcescc)
  !-----------------------------------------------------------------------
  !   This routine calculates the force term vanishing at full
  !     self-consistency. It follows the suggestion of Chan-Bohnen-Ho
  !     (PRB 47, 4771 (1993)). The true charge density is approximated
  !     by means of a free atom superposition.
  !     (alessio f.)
  ! Uses superposition of atomic charges contained in the array rho_at
  ! and read from pseudopotential files
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : tpi
  USE atom,                 ONLY : msh, rgrid
  USE uspp_param,           ONLY : upf
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE cell_base,            ONLY : tpiba
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : vnew
  USE control_flags,        ONLY : gamma_only
#ifdef USE_CUDA
  USE wavefunctions_module, ONLY : psic=>psic_d
  USE compute_rhocgnt_gpu_m
  USE cudafor
  USE gvect,                ONLY : ngm, gstart, nl=>nl_d, g=>g_d, ngl, gl=>gl_d, igtongl=>igtongl_d
#else
  USE gvect,                ONLY : ngm, gstart, nl, g, ngl, gl, igtongl
  USE wavefunctions_module, ONLY : psic
#endif 
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  implicit none
  !
  real(DP) :: forcescc (3, nat)
  !
  real(DP), allocatable :: rhocgnt (:), aux (:)
#ifdef USE_CUDA
  attributes(device):: rhocgnt
  real(DP), pointer, device :: vnew_of_r_d(:,:)
  integer    :: blocks
  type(dim3) :: threads
#endif
  ! work space
  real(DP) ::  gx, arg, fact, tmpf, tau1,tau2,tau3
  real(DP) ::  fscc1, fscc2, fscc3
  ! temp factors
  integer :: ir, isup, isdw, ig, nt, na, ipol, ndm, i
  ! counters
  !
  ! vnew is V_out - V_in, psic is the temp space

#ifdef USE_CUDA
! We need to use a pointer and lbound/ubound for the CUF kernels 
   vnew%of_r_d = vnew%of_r
   vnew_of_r_d => vnew%of_r_d
  if (nspin == 1 .or. nspin == 4) then
   !$cuf kernel do(1) <<<*,*>>>
    do i = lbound(psic,1), ubound(psic, 1)
     psic(i) = vnew_of_r_d (i, 1)
    end do
  else
   !$cuf kernel do(1) <<<*,*>>>
    do i = lbound(psic,1), ubound(psic, 1)
     psic(i) = (vnew_of_r_d (i, 1) + vnew_of_r_d (i, 2)) * 0.5d0
    end do
  end if
#else

  if (nspin == 1 .or. nspin == 4) then
     psic(:) = vnew%of_r (:, 1)
  else
     isup = 1
     isdw = 2
     psic(:) = (vnew%of_r (:, isup) + vnew%of_r (:, isdw)) * 0.5d0
  end if
  ndm = MAXVAL ( msh(1:ntyp) )
  allocate ( aux(ndm) )

#endif

  allocate ( rhocgnt(ngl) )

  forcescc(:,:) = 0.d0

  CALL fwfft ('Dense', psic, dfftp)

  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if

  do nt = 1, ntyp
     !
     ! Here we compute the G.ne.0 term
     !
#ifdef USE_CUDA
     threads = dim3(32, 8, 1)
     blocks = ceiling(real(ngl-gstart+1)/8)  ! The loop goes from gstart to ngl, so there are only ngl-gstart+1 elements
     call compute_rhocgnt_gpu <<<blocks, threads>>> (ngl-gstart+1, msh(nt), tpiba, gl(gstart),rgrid(nt)%r_d,&
       upf(nt)%rho_at_d, rgrid(nt)%rab_d, rhocgnt(gstart))
#else
     do ig = gstart, ngl
        gx = sqrt (gl (ig) ) * tpiba
        do ir = 1, msh (nt)
           if (rgrid(nt)%r(ir) .lt.1.0d-8) then
              aux (ir) = upf(nt)%rho_at (ir)
           else
              aux (ir) = upf(nt)%rho_at (ir) * &
                         sin(gx*rgrid(nt)%r(ir)) / (rgrid(nt)%r(ir)*gx)
           endif
        enddo
        call simpson (msh (nt), aux, rgrid(nt)%rab, rhocgnt (ig) )
     enddo
#endif

     do na = 1, nat
        if (nt.eq.ityp (na) ) then
           tau1 = tau(1,na)
           tau2 = tau(2,na)
           tau3 = tau(3,na)
           fscc1 = 0.d0
           fscc2 = 0.d0
           fscc3 = 0.d0
#ifdef USE_CUDA
          !$cuf kernel do(1) <<<*,*>>>
#else
          !$omp parallel do default(shared),private(arg,tmpf),reduction(+:fscc1,fscc2,fscc3)
#endif
           do ig = gstart, ngm
              arg = (g (1, ig) * tau1 + &
                     g (2, ig) * tau2 + &
                     g (3, ig) * tau3 ) * tpi
              tmpf=  fact * rhocgnt (igtongl(ig) ) * tpiba * &
                     DBLE( CMPLX(sin(arg),cos(arg),kind=DP) * CONJG(psic(nl(ig))) )
              fscc1 = fscc1 + tmpf *  g(1,ig)
              fscc2 = fscc2 + tmpf *  g(2,ig)
              fscc3 = fscc3 + tmpf *  g(3,ig)
           enddo
#ifndef USE_CUDA
          !$omp end parallel do 
#endif
           forcescc(1,na) = forcescc(1,na) + fscc1
           forcescc(2,na) = forcescc(2,na) + fscc2
           forcescc(3,na) = forcescc(3,na) + fscc3
        endif
     enddo
  enddo

  call mp_sum(  forcescc, intra_bgrp_comm )
  !
#ifdef USE_CUDA
  deallocate ( rhocgnt )
#else
  deallocate ( aux, rhocgnt )
#endif

  return
end subroutine force_corr

