!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine set_rhoc
  !-----------------------------------------------------------------------
  !
  !    This routine computes the core charge on the real space 3D mesh
  !
  !
  USE io_global, ONLY : stdout
  USE kinds,     ONLY : DP
  USE atom,      ONLY : msh, rgrid
  USE uspp_param,ONLY : upf
  USE ions_base, ONLY : ntyp => nsp
  USE cell_base, ONLY : omega, tpiba2
  USE ener,      ONLY : etxcc
  USE fft_base,  ONLY : dfftp
  USE fft_interfaces,ONLY : invfft
#ifdef USE_CUDA
  USE gvect,     ONLY : ngm, nl=>nl_d, nlm, ngl, gl=>gl_d, igtongl=>igtongl_d
  USE scf,       ONLY : rho_core_h => rho_core, rhog_core_h => rhog_core, &
                        rho_core=>rho_core_d, rhog_core=>rhog_core_d, scf_type
  USE vlocal,    ONLY : strf=>strf_d
#else
  USE gvect,     ONLY : ngm, nl, nlm, ngl, gl, igtongl
  USE scf,       ONLY : rho_core, rhog_core, scf_type
  USE vlocal,    ONLY : strf
#endif
  USE lsda_mod,  ONLY : nspin
  USE control_flags, ONLY : gamma_only
  USE mp_bands,  ONLY : intra_bgrp_comm
  USE mp,        ONLY : mp_sum
  !
  implicit none
  !
  real(DP), parameter :: eps = 1.d-10

  complex(DP) , allocatable :: aux (:)
  ! used for the fft of the core charge

  real(DP) , allocatable ::  rhocg(:)
#ifdef USE_CUDA
  attributes(device) :: aux, rhocg
#endif
  ! the radial fourier trasform
  real(DP) ::  rhoima, rhoneg, rhorea
  ! used to check the core charge
  real(DP) ::  vtxcc
  ! dummy xc energy term
  type(scf_type) :: dum

  integer :: ir, nt, ng
  ! counter on mesh points
  ! counter on atomic types
  ! counter on g vectors
  CALL start_clock( 'set_rhoc' )

  etxcc = 0.0_DP
  if ( ANY( upf(1:ntyp)%nlcc ) ) goto 10
  
  rhog_core(:) = 0.0_DP
  rho_core(:)  = 0.0_DP

  CALL stop_clock( 'set_rhoc' )
  return

10 continue
  allocate (aux( dfftp%nnr))    
  allocate (rhocg( ngl))    
  aux (:) = (0.0_DP, 0.0_DP)
  !
  !    the sum is on atom types
  !
  do nt = 1, ntyp
     if ( upf(nt)%nlcc ) then
        !
        !     drhoc compute the radial fourier transform for each shell of g vec
        !
#ifdef USE_CUDA
        call drhoc_gpu (ngl, gl, omega, tpiba2, msh(nt), rgrid(nt)%r_d, rgrid(nt)%rab_d, upf(nt)%rho_atc_d, rhocg)
#else
        call drhoc (ngl, gl, omega, tpiba2, msh (nt), rgrid(nt)%r, &
             rgrid(nt)%rab, upf(nt)%rho_atc, rhocg)
#endif
        !
        !     multiply by the structure factor and sum
        !
#ifdef USE_CUDA
        !$cuf kernel do (1) <<<*, *>>>
#endif
        do ng = 1, ngm
           aux(nl(ng)) = aux(nl(ng)) + strf(ng,nt) * rhocg(igtongl(ng))
        enddo
     endif
  enddo
  if (gamma_only) then
#ifdef USE_CUDA
     print*, "set_rhoc: GAMMA ONLY NOT IMPLEMENTED!"
     flush(6); stop
#else
     do ng = 1, ngm
        aux(nlm(ng)) = CONJG(aux(nl (ng)))
     end do
#endif
  end if
  !
#ifdef USE_CUDA
  !$cuf kernel do (1) <<<*, *>>>
#endif
  do ng = lbound(rhog_core, 1), ubound(rhog_core, 1)
    rhog_core(ng) = aux(nl(ng))
  end do
  !
  !   the core charge in real space
  !
  CALL invfft ('Dense', aux, dfftp)
  !
  !    test on the charge and computation of the core energy
  !
  rhoneg = 0.d0
  rhoima = 0.d0
#ifdef USE_CUDA
  !$cuf kernel do (1) <<<*, *>>>
#endif
  do ir = 1, dfftp%nnr
     rhoneg = rhoneg + min (0.d0,  DBLE (aux (ir) ) )
     rhoima = rhoima + abs (AIMAG (aux (ir) ) )
     rho_core(ir) =  DBLE (aux(ir))
     !
     ! NOTE: Core charge is computed in reciprocal space and brought to real
     ! space by FFT. For non smooth core charges (or insufficient cut-off)
     ! this may result in negative values in some grid points.
     ! Up to October 1999 the core charge was forced to be positive definite.
     ! This induces an error in the force, and probably stress, calculation if
     ! the number of grid points where the core charge would be otherwise neg
     ! is large. The error disappears for sufficiently high cut-off, but may be
     ! rather large and it is better to leave the core charge as it is.
     ! If you insist to have it positive definite (with the possible problems
     ! mentioned above) uncomment the following lines.  SdG, Oct 15 1999
     !
     !         rhorea = max ( DBLE (aux (ir) ), eps)
     !         rho_core(ir) = rhorea
     !
  enddo
  rhoneg = rhoneg / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
  rhoima = rhoima / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
  !
  call mp_sum(  rhoneg, intra_bgrp_comm )
  call mp_sum(  rhoima, intra_bgrp_comm )
  !
  IF (rhoneg < -1.0d-6 .OR. rhoima > 1.0d-6) &
       WRITE( stdout, '(/5x,"Check: negative/imaginary core charge=",2f12.6)')&
       rhoneg, rhoima
  !
  ! calculate core_only exch-corr energy etxcc=E_xc[rho_core] if required
  ! The term was present in previous versions of the code but it shouldn't
  !
  !   call create_scf_type(dum)
  !   dum%of_r(:,:) = 0.0_DP
  !   dum%of_g(:,:) = (0.0_DP, 0.0_DP)
  !   
  !   call v_xc( dum, rho_core, rhog_core, etxcc, vtxcc, aux )
  ! 
  !   call destroy_scf_type(dum)
  !   WRITE( stdout, 9000) etxcc
  !   WRITE( stdout,  * ) 'BEWARE it will be subtracted from total energy !'
  !
  deallocate (rhocg)
  deallocate (aux)

#ifdef USE_CUDA
  rho_core_h = rho_core
  rhog_core_h = rhog_core
#endif
  !
  CALL stop_clock( 'set_rhoc' )
  return

  ! 9000 format (5x,'core-only xc energy         = ',f15.8,' Ry')

end subroutine set_rhoc

