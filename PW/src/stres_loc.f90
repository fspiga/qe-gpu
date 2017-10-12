
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!


#ifdef USE_CUDA
 module compute_dvloc_gpu_m
 contains
  attributes(global) subroutine dvloc_of_g_gpu (mesh, msh, rab, r, vloc_at, zp, &
                     tpiba2, ngl, gl, omega, dvloc )

 !----------------------------------------------------------------------
  !
  ! dvloc = D Vloc (g^2) / D g^2 = (1/2g) * D Vloc(g) / D g
  !
  USE cudafor
  USE kinds
  USE constants , ONLY : pi, fpi, e2, eps8
  implicit none
  !
  !    first the dummy variables
  !
  integer, value,  intent(in) :: ngl, mesh, msh
  ! the number of shell of G vectors
  ! max number of mesh points
  ! number of mesh points for radial integration

  real(DP), intent(in), value :: zp, tpiba2, omega
  ! valence pseudocharge
  ! 2 pi / alat
  ! the volume of the unit cell
  real(DP), intent(in), device :: rab(mesh), r(mesh), vloc_at(mesh), gl(ngl)
  ! the derivative of the radial grid
  ! the radial grid
  ! the pseudo on the radial grid
  ! the moduli of g vectors for each s
  !
  real(DP), intent(out), device ::  dvloc (ngl)
  ! the fourier transform dVloc/dG 

  real(DP) :: vlcp, g2a

  integer :: tx, ty, igl, ir
  real(DP):: mysum, val, gx

  ! counter on erf functions or gaussians
  ! counter on g shells vectors
  ! first shell with g != 0

  tx = threadIdx%x
  ty = threadIdx%y
 
  igl = (blockIdx%x - 1) * blockDim%y + ty

  if (igl > ngl) return

  gx = sqrt(gl(igl) * tpiba2)
  mysum = 0.d0

  do ir = tx, mesh, blockDim%x 
    val = ( r(ir) * vloc_at(ir) + zp * e2 * erf(r(ir) ) ) * &
           (r(ir) * cos (gx * r(ir) ) / gx - sin (gx  * r(ir) ) / gx**2) * &
          rab(ir) 

    if (ir == 1 .or. ir == mesh) then
      mysum = mysum + val
    else if (mod(ir,2)) then
      mysum = mysum + 2.d0*val
    else
      mysum = mysum + 4.d0*val
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
         mysum  =  ( mysum / 3.d0 ) * fpi / omega / 2.0d0 / gx 
         g2a = gl (igl) * tpiba2 / 4.d0
         mysum = mysum + fpi / omega * zp * e2 * exp ( - g2a) * (g2a + &
          1.d0) / (gl (igl) * tpiba2) **2
         
         dvloc(igl) =  mysum 
       endif  

  ! the  G=0 component is not computed
  if (igl ==1 .and. gl(1) <eps8)   dvloc(1)=0.d0

  end subroutine dvloc_of_g_gpu

 end module compute_dvloc_gpu_m
#endif
!----------------------------------------------------------------------
subroutine stres_loc (sigmaloc)
  !----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE atom,                 ONLY : msh, rgrid
  USE m_gth,                ONLY : dvloc_gth
  USE ions_base,            ONLY : ntyp => nsp
  USE cell_base,            ONLY : omega, tpiba2
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
#ifdef USE_CUDA
  USE cudafor
  USE compute_dvloc_gpu_m
  USE wavefunctions_module, ONLY : psic=>psic_d
  USE gvect,                ONLY : ngm, gstart, nl=>nl_d, g=>g_d, ngl, &
                                   gl=>gl_d, igtongl=>igtongl_d
  USE vlocal,               ONLY : strf=>strf_d, vloc=>vloc_d
#else
  USE wavefunctions_module, ONLY : psic
  USE gvect,                ONLY : ngm, gstart, nl, g, ngl, gl, igtongl
  USE vlocal,               ONLY : strf, vloc
#endif
  USE lsda_mod,             ONLY : nspin
  USE scf,                  ONLY : rho
  USE control_flags,        ONLY : gamma_only
  USE uspp_param,           ONLY : upf
  USE noncollin_module,     ONLY : nspin_lsda
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  !
  implicit none
  !
  real(DP) :: sigmaloc (3, 3)
  real(DP) , allocatable :: dvloc(:)
  real(DP) :: evloc, fact
  integer :: ng, nt, l, m, is
  ! counter on g vectors
  ! counter on atomic type
  ! counter on angular momentum
  ! counter on spin components

#ifdef USE_CUDA
  attributes(device):: dvloc
  real(DP), pointer, device :: rho_of_r_d(:,:)
  integer    :: blocks
  type(dim3) :: threads
  real(DP):: s11,s21,s22,s31,s32,s33,tmpf
#endif
  allocate(dvloc(ngl))
  sigmaloc(:,:) = 0.d0

#ifdef USE_CUDA
! We need to use a pointer and lbound/ubound for the CUF kernels 
 !$cuf kernel do(1) <<<*,*>>>
    do l = lbound(psic,1), ubound(psic, 1)
      psic(l) = (0.d0,0.d0)
    end do
 
   rho%of_r_d=rho%of_r
   rho_of_r_d =>rho%of_r_d

! call daxpy (dfftp%nnr, 1.d0, rho%of_r (1, is), 1, psic, 2)
 
 !$cuf kernel do(1) <<<*,*>>>
  do l=1,dfftp%nnr
   do is = 1, nspin_lsda
    psic(l) = psic(l) + rho_of_r_d (l, is)
   enddo
  enddo

#else

  psic(:)=(0.d0,0.d0)
  do is = 1, nspin_lsda
     call daxpy (dfftp%nnr, 1.d0, rho%of_r (1, is), 1, psic, 2)
  enddo

#endif

  CALL fwfft ('Dense', psic, dfftp)
  ! psic contains now the charge density in G space
  if (gamma_only) then
     fact = 2.d0
  else
     fact = 1.d0
  end if

  evloc = 0.0d0
  do nt = 1, ntyp
#ifdef USE_CUDA
 !$cuf kernel do(1) <<<*,*>>>
#endif
     do ng = gstart, ngm
        if (ng==2 .and. gstart==2) evloc = evloc + &
                 DBLE( psic (nl (1) ) * strf (1, nt) ) * vloc (igtongl (1), nt)
        evloc = evloc +  DBLE (CONJG(psic (nl (ng) ) ) * strf (ng, nt) ) &
             * vloc (igtongl (ng), nt) * fact
     enddo
  enddo
  
!        WRITE( 6,*) ' evloc ', evloc, evloc*omega, gstart   ! DEBUG
  
  do nt = 1, ntyp
     IF ( .NOT. ASSOCIATED ( upf(nt)%vloc ) ) THEN
        IF ( upf(nt)%is_gth ) THEN
           !
           ! special case: GTH pseudopotential
           !
#ifdef USE_CUDA
           print *,"dvloc_gth not yet implemented on GPU"
#else
           call dvloc_gth (nt, upf(nt)%zp, tpiba2, ngl, gl, omega, dvloc)
#endif
           !
        ELSE
           !
           ! special case: pseudopotential is coulomb 1/r potential
           !
#ifdef USE_CUDA
           print *,"dvloc_coul not yet implemented on GPU"
#else
           call dvloc_coul (upf(nt)%zp, tpiba2, ngl, gl, omega, dvloc)
#endif
           !
        END IF
     ELSE
        !
        ! normal case: dvloc contains dV_loc(G)/dG
        !
#ifdef USE_CUDA

        threads = dim3(32, 8, 1)
        blocks = ceiling(real(ngl)/8)  ! The loop goes from 1 to ngl,
        call dvloc_of_g_gpu <<<blocks, threads>>>  (rgrid(nt)%mesh, msh (nt), rgrid(nt)%rab_d, rgrid(nt)%r_d,&
          upf(nt)%vloc_d, upf(nt)%zp, tpiba2, ngl, gl, omega, dvloc)
#else
        call dvloc_of_g (rgrid(nt)%mesh, msh (nt), rgrid(nt)%rab, rgrid(nt)%r,&
          upf(nt)%vloc(:), upf(nt)%zp, tpiba2, ngl, gl, omega, dvloc)
#endif
        !
     END IF
     ! no G=0 contribution
#ifdef USE_CUDA
      s11=0.d0
      s21=0.d0
      s22=0.d0
      s31=0.d0
      s32=0.d0
      s33=0.d0
 !$cuf kernel do(1) <<<*,*>>>
     do ng = 1, ngm
              tmpf =  DBLE( CONJG( psic(nl(ng) ) )  * strf(ng, nt) ) * &
                      2.0d0  * dvloc(igtongl(ng) ) * tpiba2 * fact
              s11= s11+ tmpf* g(1,ng) *g(1,ng)
              s21= s21+ tmpf* g(2,ng) *g(1,ng)
              s22= s22+ tmpf* g(2,ng) *g(2,ng)
              s31= s31+ tmpf* g(3,ng) *g(1,ng)
              s32= s32+ tmpf* g(3,ng) *g(2,ng)
              s33= s33+ tmpf* g(3,ng) *g(3,ng)
     enddo
     sigmaloc(1,1) = sigmaloc(1,1) +s11 
     sigmaloc(2,1) = sigmaloc(2,1) +s21 
     sigmaloc(2,2) = sigmaloc(2,2) +s22 
     sigmaloc(3,1) = sigmaloc(3,1) +s31 
     sigmaloc(3,2) = sigmaloc(3,2) +s32 
     sigmaloc(3,3) = sigmaloc(3,3) +s33 
#else
     do ng = 1, ngm
        do l = 1, 3
           do m = 1, l
              sigmaloc(l, m) = sigmaloc(l, m) +  DBLE( CONJG( psic(nl(ng) ) ) &
                    * strf (ng, nt) ) * 2.0d0 * dvloc (igtongl (ng) ) &
                    * tpiba2 * g (l, ng) * g (m, ng) * fact
           enddo
        enddo
     enddo
#endif
  enddo
  !
  do l = 1, 3
     sigmaloc (l, l) = sigmaloc (l, l) + evloc
     do m = 1, l - 1
        sigmaloc (m, l) = sigmaloc (l, m)
     enddo
  enddo
  !
  call mp_sum(  sigmaloc, intra_bgrp_comm )
  !
  deallocate(dvloc)
  return
end subroutine stres_loc
