
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! TB
! included deallocation of forcefield of monopole 'forcemono'
!
!----------------------------------------------------------------------
SUBROUTINE clean_pw( lflag )
  !----------------------------------------------------------------------
  !    
  ! ... This routine deallocates dynamically allocated arrays
  ! ... if lflag=.TRUE.  all arrays are deallocated (end of calculation)
  ! ... if lflag=.FALSE. ion-related variables and arrays allocated
  ! ... at the very beginning of the calculation (routines iosys, read_file,
  ! ... setup, read_pseudo) are not deallocated; all others arrays are.
  ! ... This is used when a new calculation has to be performed (e.g. in neb,
  ! ... phonon, vc-relax). Beware: the new calculation should not call any
  ! ... of the routines mentioned above
  !
  USE basis,                ONLY : swfcatom
  USE cellmd,               ONLY : lmovecell
  USE ions_base,            ONLY : deallocate_ions_base
  USE gvect,                ONLY : g, gg, gl, nl, nlm, igtongl, mill, &
                                   eigts1, eigts2, eigts3
  USE gvecs,                ONLY : nls, nlsm
  USE fixed_occ,            ONLY : f_inp
  USE ktetra,               ONLY : deallocate_tetra
  USE klist,                ONLY : deallocate_igk
  USE gvect,                ONLY : ig_l2g
  USE vlocal,               ONLY : strf, vloc
  USE wvfct,                ONLY : g2kin, et, wg, btype
  USE force_mod,            ONLY : force
  USE scf,                  ONLY : rho, v, vltot, rho_core, rhog_core, &
                                   vrs, kedtau, destroy_scf_type, vnew
  USE symm_base,            ONLY : irt
  USE symme,                ONLY : sym_rho_deallocate
  USE wavefunctions_module, ONLY : evc, psic, psic_batch, psic_nc
  USE us,                   ONLY : qrad, tab, tab_at, tab_d2y, spline_ps
  USE uspp,                 ONLY : deallocate_uspp
  USE uspp_param,           ONLY : upf
  USE m_gth,                ONLY : deallocate_gth
  USE ldaU,                 ONLY : deallocate_ldaU
  USE extfield,             ONLY : forcefield, forcemono
  USE fft_base,             ONLY : dfftp, dffts  
  USE fft_base,             ONLY : pstickdealloc
  USE fft_types,            ONLY : fft_type_deallocate
  USE spin_orb,             ONLY : lspinorb, fcoef
  USE noncollin_module,     ONLY : deallocate_noncol
  USE dynamics_module,      ONLY : deallocate_dyn_vars
  USE paw_init,             ONLY : deallocate_paw_internals
  USE atom,                 ONLY : msh, rgrid
  USE radial_grids,         ONLY : deallocate_radial_grid
  USE wannier_new,          ONLY : use_wannier
  !
  USE london_module,        ONLY : dealloca_london
  USE xdm_module,           ONLY : cleanup_xdm
  USE constraints_module,   ONLY : deallocate_constraint
  USE realus,               ONLY : deallocate_realsp
  USE pseudo_types,         ONLY : deallocate_pseudo_upf
  USE bp,                   ONLY : deallocate_bp_efield
  USE exx,                  ONLY : deallocate_exx
  !
  USE control_flags,        ONLY : ts_vdw
  USE tsvdw_module,         ONLY : tsvdw_finalize
  !
#ifdef USE_CUDA
  USE vlocal,               ONLY : strf_d,vloc_d
  USE wavefunctions_module, ONLY : evc_d, psic_d, psic_batch_d
  USE gvect,                ONLY : g_d, gg_d, nl_d, mill_d, &
                                   eigts1_d, eigts2_d, eigts3_d, igtongl_d, gl_d
  USE gvecs,                ONLY : nls_d
  USE wvfct,                ONLY : g2kin_d, et_d, wg_d
  USE us,                   ONLY : qrad_d, tab_d, tab_d2y_d
  USE fft_parallel,         ONLY : f_h, aux_h, aux_d, aux2_h, aux2_d
  USE scf,                  ONLY : vltot_d, vrs_d, rho_core_d, rhog_core_d
  USE wvfct,                ONLY : psi_d, hpsi_d, spsi_d, comm_h_c, comm_s_c
  USE wvfct,                ONLY : hc_d, sc_d, vc_d, vc_temp_d, ew_d, conv_d, conv_idx_d
#endif
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: lflag
  !
  INTEGER :: nt, nr1,nr2,nr3
  !
  IF ( lflag ) THEN
     !
     ! ... arrays allocated at the very beginning of the calculation
     !
     IF( ALLOCATED( upf ) ) THEN
        DO nt = 1, SIZE( upf )
           CALL deallocate_pseudo_upf( upf( nt ) )
        END DO
        DEALLOCATE( upf )
     END IF
     IF (ALLOCATED(msh)) DEALLOCATE (msh)
     CALL deallocate_radial_grid(rgrid)
     !
     CALL deallocate_ions_base()
     !
     IF ( ALLOCATED( force ) )      DEALLOCATE( force )
     IF ( ALLOCATED( forcefield ) ) DEALLOCATE( forcefield )
     IF ( ALLOCATED( forcemono ) )  DEALLOCATE( forcemono )
     IF ( ALLOCATED( irt ) )        DEALLOCATE( irt )
     !
     CALL dealloca_london()
     CALL cleanup_xdm()
     CALL deallocate_constraint()
     !
  END IF
  !
  CALL deallocate_bp_efield()
  !
  CALL deallocate_ldaU ( lflag )
  !
  CALL deallocate_tetra ( )
  !
  IF ( ALLOCATED( f_inp ) .and. lflag )      DEALLOCATE( f_inp )
  !
  ! ... arrays allocated in ggen.f90
  !
  IF ( ALLOCATED( ig_l2g ) )     DEALLOCATE( ig_l2g )
  IF ( .NOT. lmovecell ) THEN
     IF ( ASSOCIATED( gl ) )     DEALLOCATE ( gl )
#ifdef USE_CUDA
     IF ( ASSOCIATED( gl_d ) )   DEALLOCATE ( gl_d )
#endif
  END IF
  !
  CALL sym_rho_deallocate ( )
  !
  ! ... arrays allocated in allocate_fft.f90 ( and never deallocated )
  !
  IF ( ALLOCATED( g ) )          DEALLOCATE( g )
  IF ( ALLOCATED( gg ) )         DEALLOCATE( gg )
  IF ( ALLOCATED( nl ) )         DEALLOCATE( nl )  
  IF ( ALLOCATED( nlm ) )        DEALLOCATE( nlm )
  IF ( ALLOCATED( igtongl ) )    DEALLOCATE( igtongl )  
  IF ( ALLOCATED( mill ) )       DEALLOCATE( mill )
#ifdef USE_CUDA
  IF ( ALLOCATED( g_d ) )          DEALLOCATE( g_d )
  IF ( ALLOCATED( gg_d ) )         DEALLOCATE( gg_d )
  IF ( ALLOCATED( nl_d ) )         DEALLOCATE( nl_d )
  IF ( ALLOCATED( igtongl_d ) )    DEALLOCATE( igtongl_d )  
  IF ( ALLOCATED( mill_d ) )       DEALLOCATE( mill_d )
#endif
  call destroy_scf_type(rho)
  call destroy_scf_type(v)
  call destroy_scf_type(vnew)
  IF ( ALLOCATED( kedtau ) )     DEALLOCATE( kedtau )
  IF ( ALLOCATED( vltot ) )      DEALLOCATE( vltot )
  IF ( ALLOCATED( rho_core ) )   DEALLOCATE( rho_core )
  IF ( ALLOCATED( rhog_core ) )  DEALLOCATE( rhog_core )
  IF ( ALLOCATED( psic ) )       DEALLOCATE( psic )
  IF ( ALLOCATED (psic_batch ) ) DEALLOCATE( psic_batch )
  IF ( ALLOCATED( psic_nc ) )    DEALLOCATE( psic_nc )
  IF ( ALLOCATED( vrs ) )        DEALLOCATE( vrs )
  if (spline_ps) then
    IF ( ALLOCATED( tab_d2y) )     DEALLOCATE( tab_d2y )
#ifdef USE_CUDA
    IF ( ALLOCATED( tab_d2y_d) )     DEALLOCATE( tab_d2y_d )
#endif
  endif
  IF ( ALLOCATED( nls ) )     DEALLOCATE( nls )
  IF ( ALLOCATED( nlsm ) )   DEALLOCATE( nlsm )

#ifdef USE_CUDA
  IF ( ALLOCATED( psic_d ) )        DEALLOCATE( psic_d )
  IF ( ALLOCATED( nls_d ) )         DEALLOCATE( nls_d )
  IF ( ALLOCATED( psic_batch_d ) )  DEALLOCATE( psic_batch_d )
  IF ( ALLOCATED( f_h ) )           DEALLOCATE( f_h )
  IF ( ALLOCATED( aux_h ) )         DEALLOCATE( aux_h )
  IF ( ALLOCATED( aux_d ) )         DEALLOCATE( aux_d )
  IF ( ALLOCATED( aux2_h ) )        DEALLOCATE( aux2_h )
  IF ( ALLOCATED( aux2_d ) )        DEALLOCATE( aux2_d )
  IF ( ALLOCATED( vltot_d ) )       DEALLOCATE( vltot_d )
  IF ( ALLOCATED( rho_core_d ) )    DEALLOCATE( rho_core_d )
  IF ( ALLOCATED( rhog_core_d ) )   DEALLOCATE( rhog_core_d )
  IF ( ALLOCATED( vrs_d ) )         DEALLOCATE( vrs_d )
#endif
  !
  ! ... arrays allocated in allocate_locpot.f90 ( and never deallocated )
  !
  IF ( ALLOCATED( vloc ) )       DEALLOCATE( vloc )
  IF ( ALLOCATED( strf ) )       DEALLOCATE( strf )
  IF ( ALLOCATED( eigts1 ) )     DEALLOCATE( eigts1 )
  IF ( ALLOCATED( eigts2 ) )     DEALLOCATE( eigts2 )
  IF ( ALLOCATED( eigts3 ) )     DEALLOCATE( eigts3 )
#ifdef USE_CUDA
  IF ( ALLOCATED( vloc_d ) )       DEALLOCATE( vloc_d )
  IF ( ALLOCATED( strf_d ) )       DEALLOCATE( strf_d )
  IF ( ALLOCATED( eigts1_d ) )     DEALLOCATE( eigts1_d )
  IF ( ALLOCATED( eigts2_d ) )     DEALLOCATE( eigts2_d )
  IF ( ALLOCATED( eigts3_d ) )     DEALLOCATE( eigts3_d )
#endif
  !
  ! ... arrays allocated in allocate_nlpot.f90 ( and never deallocated )
  !
  IF ( ALLOCATED( g2kin ) )      DEALLOCATE( g2kin )
  IF ( ALLOCATED( qrad ) )       DEALLOCATE( qrad )
  IF ( ALLOCATED( tab ) )        DEALLOCATE( tab )
  IF ( ALLOCATED( tab_at ) )     DEALLOCATE( tab_at )
  IF ( lspinorb ) THEN
     IF ( ALLOCATED( fcoef ) )   DEALLOCATE( fcoef )
  END IF
  !
#ifdef USE_CUDA
  IF ( ALLOCATED( g2kin_d ) )      DEALLOCATE( g2kin_d )
  IF ( ALLOCATED( qrad_d ) )       DEALLOCATE( qrad_d )
  IF ( ALLOCATED( tab_d ) )        DEALLOCATE( tab_d )
#endif
  CALL deallocate_igk ( )
  CALL deallocate_uspp() 
  CALL deallocate_gth( lflag ) 
  CALL deallocate_noncol() 
  !
  ! ... arrays allocated in init_run.f90 ( and never deallocated )
  !
  IF ( ALLOCATED( et ) )         DEALLOCATE( et )
  IF ( ALLOCATED( wg ) )         DEALLOCATE( wg )
  IF ( ALLOCATED( btype ) )      DEALLOCATE( btype )

#ifdef USE_CUDA
  IF ( ALLOCATED( et_d ) )         DEALLOCATE( et_d )
  IF ( ALLOCATED( wg_d ) )         DEALLOCATE( wg_d )
#endif

  !
  ! ... arrays allocated in allocate_wfc.f90 ( and never deallocated )
  !
  IF ( ALLOCATED( evc ) )        DEALLOCATE( evc )
  IF ( ALLOCATED( swfcatom ) )   DEALLOCATE( swfcatom )
#ifdef USE_CUDA
  IF ( ALLOCATED( evc_d ) )      DEALLOCATE( evc_d )
  IF ( ALLOCATED( psi_d ) )      DEALLOCATE( psi_d )
  IF ( ALLOCATED( hpsi_d ) )     DEALLOCATE( hpsi_d )
  IF ( ALLOCATED( spsi_d ) )     DEALLOCATE( spsi_d )
  IF ( ALLOCATED( comm_h_c ) )   DEALLOCATE( comm_h_c )
  IF ( ALLOCATED( comm_s_c ) )   DEALLOCATE( comm_s_c )
  IF ( ALLOCATED( hc_d ) )       DEALLOCATE( hc_d )
  IF ( ALLOCATED( sc_d ) )       DEALLOCATE( sc_d )
  IF ( ALLOCATED( vc_d ) )       DEALLOCATE( vc_d )
  IF ( ALLOCATED( vc_temp_d ) )  DEALLOCATE( vc_temp_d )
  IF ( ALLOCATED( ew_d ) )       DEALLOCATE( ew_d ) 
  IF ( ALLOCATED( conv_d ) )     DEALLOCATE( conv_d )
  IF ( ALLOCATED( conv_idx_d ) ) DEALLOCATE( conv_idx_d )

#endif
  !
  ! ... fft structures allocated in data_structure.f90  
  !
  ! UGLY HACK WARNING: unlike previous versions, fft_type_deallocate
  ! removes all information about FFT grids, including FFT dimensions.
  ! If however FFT dimensions were set from input data, one may end
  ! up with a different grid if FFT grids are re-initialized later.
  ! The following workaround restores the previous functionality.
  ! TODO: replace clean_pw with more fine-grained cleaning routines.
  nr1 = dfftp%nr1; nr2 = dfftp%nr2; nr3 = dfftp%nr3
  CALL fft_type_deallocate( dfftp )
  dfftp%nr1 = nr1; dfftp%nr2 = nr2; dfftp%nr3 = nr3
  !
  nr1 = dffts%nr1; nr2 = dffts%nr2; nr3 = dffts%nr3
  CALL fft_type_deallocate( dffts )
  dffts%nr1 = nr1; dffts%nr2 = nr2; dffts%nr3 = nr3
  !
  ! ... stick-owner matrix allocated in sticks_base
  !
  CALL pstickdealloc()
  !
  ! ... arrays allocated for dynamics
  !
  CALL deallocate_dyn_vars()
  !
  ! ... additional arrays for PAW
  !
  CALL deallocate_paw_internals()
  !
  ! ... arrays for real-space algorithm
  !
  CALL  deallocate_realsp()
  !
  ! for Wannier_ac
  if (use_wannier) CALL wannier_clean()
  !
  CALL deallocate_exx ( ) 
  !
  IF (ts_vdw) CALL tsvdw_finalize()
  !
  CALL plugin_clean( lflag )
  !
  RETURN
  !
END SUBROUTINE clean_pw
