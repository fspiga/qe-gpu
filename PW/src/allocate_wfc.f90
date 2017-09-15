!
! Copyright (C) 2001-2008 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE allocate_wfc()
  !----------------------------------------------------------------------------
  !
  ! ... dynamical allocation of arrays: wavefunctions
  ! ... must be called after allocate_nlpot 
  !
  USE io_global, ONLY : stdout
  USE wvfct,     ONLY : npwx, nbnd
  USE basis,     ONLY : natomwfc, swfcatom
  USE fixed_occ, ONLY : one_atom_occupations
  USE ldaU,      ONLY : wfcU, nwfcU, lda_plus_u, U_projection
  USE noncollin_module,     ONLY : noncolin, npol
  USE wavefunctions_module, ONLY : evc
#ifdef USE_CUDA                                               
  USE wavefunctions_module, ONLY : evc_d                      
  USE wvfct,                ONLY : nbnd, nbndx, psi_d, hpsi_d, spsi_d, comm_h_c, comm_s_c
  USE wvfct,  ONLY : hc_d, sc_d, vc_d, vc_temp_d, ew_d, conv_d, conv_idx_d
  USE uspp,                 ONLY : okvan
#endif
  USE wannier_new, ONLY : use_wannier
  !
  IMPLICIT NONE
  !
  !
  ALLOCATE( evc( npwx*npol, nbnd ) )    
#ifdef USE_CUDA                                               
  ALLOCATE( evc_d( npwx*npol, nbnd ) )  
  ALLOCATE( psi_d( npwx, npol, nbndx ) )
  ALLOCATE( hpsi_d( npwx, npol, nbndx ) )
  IF(okvan) &
     ALLOCATE( spsi_d( npwx, npol, nbndx ) )
  ALLOCATE( comm_h_c( nbndx, nbndx ) )
  ALLOCATE( comm_s_c( nbndx, nbndx ) )
  ALLOCATE(hc_d(nbndx, nbndx))
  ALLOCATE(sc_d(nbndx, nbndx))
  ALLOCATE(vc_d(nbndx, nbndx))
  ALLOCATE(vc_temp_d(nbndx, nbndx))
  ALLOCATE(ew_d(nbndx))
  ALLOCATE(conv_d(nbndx))
  ALLOCATE(conv_idx_d(nbnd))
#endif
  IF ( one_atom_occupations .OR. use_wannier ) &
     ALLOCATE( swfcatom( npwx*npol, natomwfc) )
  IF ( lda_plus_u .AND. (U_projection.NE.'pseudo') ) &
     ALLOCATE( wfcU(npwx*npol, nwfcU) )
  !
  RETURN
  !
END SUBROUTINE allocate_wfc
