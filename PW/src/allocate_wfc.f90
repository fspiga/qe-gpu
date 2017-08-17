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
  USE wvfct,                ONLY : nbndx, psi_d, hpsi_d, spsi_d, comm_h_c
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
#endif
  IF ( one_atom_occupations .OR. use_wannier ) &
     ALLOCATE( swfcatom( npwx*npol, natomwfc) )
  IF ( lda_plus_u .AND. (U_projection.NE.'pseudo') ) &
     ALLOCATE( wfcU(npwx*npol, nwfcU) )
  !
  RETURN
  !
END SUBROUTINE allocate_wfc
