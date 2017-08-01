!
! Copyright (C) 2003-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
#ifdef USE_CUDA
SUBROUTINE myDdot( n, A, res )
  !
  use kinds, ONLY : DP
  use cudafor
  use cublas
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n
  !
!!!!!pgi$ ignore_tkr A
  REAL(DP), DEVICE, INTENT(IN) :: A(n)
  !
  REAL(DP), INTENT(OUT) :: res
  !
  res = cublasDdot( n, A, 1, A, 1 )
  !
  return
  !
END SUBROUTINE myDdot
#endif
!
!----------------------------------------------------------------------------
#ifdef USE_CUDA
SUBROUTINE regterg( npw, npwx, nvec, nvecx, evc, evc_d, ethr, &
                    uspp, gstart, e, e_d, btype, notcnv, lrot, dav_iter )
#else
SUBROUTINE regterg( npw, npwx, nvec, nvecx, evc, ethr, &
                    uspp, gstart, e, btype, notcnv, lrot, dav_iter )
#endif
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an uspp matrix, evc is a complex vector
  ! ... (real wavefunctions with only half plane waves stored)
  !
  USE kinds,         ONLY : DP
  USE io_global,     ONLY : stdout
  USE mp_bands,      ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id
  USE mp,            ONLY : mp_sum, mp_bcast
  USE cpu_gpu_interface
#ifdef USE_CUDA
  USE cudafor
  USE cublas!,        ONLY : cublasZgemm, cublasDdot
!  USE ep_debug, ONLY : compare, MPI_Wtime
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx, gstart
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set
    !    (the basis set is refreshed when its dimension would exceed nvecx)
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx,nvec)
    !  evc   contains the  refined estimates of the eigenvectors
#ifdef USE_CUDA
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: evc_d(npwx,npol,nvec)
  REAL(DP), DEVICE, INTENT(OUT) :: e_d(nvec)
#endif
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence: root improvement is stopped,
    ! when two consecutive estimates of the root differ by less than ethr.
  LOGICAL, INTENT(IN) :: uspp
    ! if .FALSE. : S|psi> not needed
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  LOGICAL, INTENT(IN) :: lrot
    ! .TRUE. if the wfc have already been rotated
  REAL(DP), INTENT(OUT) :: e(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer  number of iterations performed
    ! number of unconverged roots
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, n, m, nb1, ibnd
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
    ! counter on the bands
  INTEGER :: ierr
  REAL(DP), ALLOCATABLE :: hr(:,:), sr(:,:), vr(:,:), ew(:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! eigenvectors of the Hamiltonian
    ! eigenvalues of the reduced hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  LOGICAL, ALLOCATABLE :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr 
    ! threshold for empty bands
  INTEGER :: npw2, npwx2
  !
#ifndef USE_CUDA
  REAL(DP), EXTERNAL :: ddot
#endif
#ifdef USE_CUDA
  attributes(pinned) :: vr
  COMPLEX(DP), DEVICE, ALLOCATABLE :: hr_d(:,:), sr_d(:,:), vr_d(:,:), vr_temp_d(:,:)
  COMPLEX(DP), PINNED, ALLOCATABLE :: comm_h_c(:,:)
  COMPLEX(DP), DEVICE, ALLOCATABLE :: psi_d(:,:,:), hpsi_d(:,:,:), spsi_d(:,:,:)
  REAL(DP), DEVICE, ALLOCATABLE :: ew_d(:)
  REAL(DP), PINNED, ALLOCATABLE :: comm_h_r(:)
  LOGICAL, DEVICE, ALLOCATABLE :: conv_d(:)
  INTEGER, DEVICE :: conv_idx_d(nvec)
  INTEGER :: istat, i, j, k
  INTEGER(KIND=8) :: freeMem,totalMem
  REAL(DP) :: rFreeMem,rUsedMem,rTotalMem,minUsedMem,maxUsedMem,minFreeMem,maxFreeMem
  REAL(DP) :: ew_temp1, ew_temp2
!  REAL(DP) :: timer !, it_timer
#endif
  !
  ! EXTERNAL  h_psi, s_psi, g_psi
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi> 
    ! s_psi(npwx,npw,nvec,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,nvec)
    ! g_psi(npwx,npw,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  CALL start_clock( 'regterg' )
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'regter', 'nvecx is too small', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  ALLOCATE( psi(  npwx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate psi ', ABS(ierr) )
  ALLOCATE( hpsi( npwx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate hpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     ALLOCATE( spsi( npwx, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' regterg ',' cannot allocate spsi ', ABS(ierr) )
#ifdef USE_CUDA
     ALLOCATE( spsi_d( npwx, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' cegterg ',' cannot allocate spsi_d ', ABS(ierr) )
#endif
  END IF
  !
  ALLOCATE( sr( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate sr ', ABS(ierr) )
  ALLOCATE( hr( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate hr ', ABS(ierr) )
  ALLOCATE( vr( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate vr ', ABS(ierr) )
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate ew ', ABS(ierr) )
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'regterg ',' cannot allocate conv ', ABS(ierr) )
  !
  npw2  = 2*npw
  npwx2  = 2*npwx
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
#ifdef USE_CUDA

  ALLOCATE( ew_d( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate ew_d ', ABS(ierr) )

  ALLOCATE( comm_h_r( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate comm_h_r ', ABS(ierr) )

  ALLOCATE( conv_d( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate conv_d ', ABS(ierr) )

  ALLOCATE( sr_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate sc_d ', ABS(ierr) )
  ALLOCATE( hr_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate hc_d ', ABS(ierr) )
  ALLOCATE( vr_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate vc_d ', ABS(ierr) )

  ALLOCATE( vr_temp_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate vc_temp_d ', ABS(ierr) )

  ALLOCATE( comm_h_c( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate comm_h_c ', ABS(ierr) )

!#if 0
!      istat=CudaMemGetInfo(freeMem,totalMem)
!      rTotalMem = totalMem/(10.**6)
!      rFreeMem = freeMem/(10.**6);            MinFreeMem = rFreeMem; MaxFreeMem = rFreeMem
!      rUsedMem = (totalMem-freeMem)/(10.**6); MaxUsedMem = rUsedMem; MinUsedMem = rUsedMem
!write(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " GPU memory used: ",minUsedMem," - ",maxUsedMem," / ",rTotalMem," MBytes"
!write(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " GPU memory free: ",minFreeMem," - ",maxFreeMem," / ",rTotalMem," MBytes"
!print *," "
!call flush(6)
!print *,"psi_d: ",npwx,npol,nvecx
!call flush(6)
!#endif

  ALLOCATE(  psi_d( npwx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate psi_d ', ABS(ierr) )

!#if 0
!      istat=CudaMemGetInfo(freeMem,totalMem)
!      rTotalMem = totalMem/(10.**6)
!      rFreeMem = freeMem/(10.**6);            MinFreeMem = rFreeMem; MaxFreeMem = rFreeMem
!      rUsedMem = (totalMem-freeMem)/(10.**6); MaxUsedMem = rUsedMem; MinUsedMem = rUsedMem
!write(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " GPU memory used: ",minUsedMem," - ",maxUsedMem," / ",rTotalMem," MBytes"
!write(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " GPU memory free: ",minFreeMem," - ",maxFreeMem," / ",rTotalMem," MBytes"
!print *," "
!call flush(6)
!#endif

  ALLOCATE( hpsi_d( npwx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate hpsi_d ', ABS(ierr) )
  !
      istat=CudaMemGetInfo(freeMem,totalMem)
      rTotalMem = totalMem/(10.**6)
      rFreeMem = freeMem/(10.**6);            MinFreeMem = rFreeMem; MaxFreeMem = rFreeMem
      rUsedMem = (totalMem-freeMem)/(10.**6); MaxUsedMem = rUsedMem; MinUsedMem = rUsedMem
!write(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " GPU memory used: ",minUsedMem," - ",maxUsedMem," / ",rTotalMem," MBytes"
!write(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " GPU memory free: ",minFreeMem," - ",maxFreeMem," / ",rTotalMem," MBytes"
!print *," "
!call flush(6)
#endif


  !
  IF ( uspp ) then
#ifdef USE_CUDA
    spsi_d = ZERO
#else
    spsi = ZERO
#endif
  END IF
  !

#ifdef USE_CUDA
  hpsi_d = ZERO
  psi_d  = ZERO
 !$cuf kernel do(2) <<<*,*>>>
  Do j=1,nvec
      Do i=lbound(psi_d,1), ubound(psi_d,1)
        psi_d(i,j) = evc_d(i,j)
      end do
  end do

  IF ( gstart == 2 ) psi_d(1,1:nvec) = CMPLX( DBLE( psi_d(1,1:nvec) ), 0.D0 ,kind=DP) !TODO
  CALL h_psi( npwx, npw, nvec, psi_d, hpsi_d )
#else
  hpsi = ZERO
  psi  = ZERO
  psi(:,1:nvec) = evc(:,1:nvec)
  ! ... set Im[ psi(G=0) ] -  needed for numerical stability
  IF ( gstart == 2 ) psi(1,1:nvec) = CMPLX( DBLE( psi(1,1:nvec) ), 0.D0 ,kind=DP)
  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL h_psi( npwx, npw, nvec, psi, hpsi )

#endif

  !
  IF ( uspp ) THEN
#ifdef USE_CUDA
!    call compare(spsi, spsi_d, "spsi")
    CALL s_psi( npwx, npw, nvec, psi_d, spsi_d )
!    spsi = spsi_d
#else
    CALL s_psi( npwx, npw, nvec, psi, spsi )
#endif
  END IF
  !
  ! ... hr contains the projection of the hamiltonian onto the reduced
  ! ... space vr contains the eigenvectors of hr
  !
#ifdef USE_CUDA
  hr_d(:,:) = ZERO
  sr_d(:,:) = ZERO
  vr_d(:,:) = ZERO
#else
  hr(:,:) = ZERO
  sr(:,:) = ZERO
  vr(:,:) = ZERO
#endif
  !
#ifdef USE_CUDA
  ! TODO
#else
  CALL DGEMM( 'T', 'N', nbase, nbase, npw2, 2.D0 , &
              psi, npwx2, hpsi, npwx2, 0.D0, hr, nvecx )
  !
  IF ( gstart == 2 ) &
     CALL DGER( nbase, nbase, -1.D0, psi, npwx2, hpsi, npwx2, hr, nvecx )
  !
  CALL mp_sum( hr( :, 1:nbase ), intra_bgrp_comm )
#endif

  !
#ifdef USE_CUDA
  !TODO
#else
  IF ( uspp ) THEN
     !
     CALL DGEMM( 'T', 'N', nbase, nbase, npw2, 2.D0, &
                 psi, npwx2, spsi, npwx2, 0.D0, sr, nvecx )
     !
     IF ( gstart == 2 ) &
        CALL DGER( nbase, nbase, -1.D0, psi, npwx2, spsi, npwx2, sr, nvecx )
     !
  ELSE
     !
     CALL DGEMM( 'T', 'N', nbase, nbase, npw2, 2.D0, &
                 psi, npwx2, psi, npwx2, 0.D0, sr, nvecx )
     !
     IF ( gstart == 2 ) &
        CALL DGER( nbase, nbase, -1.D0, psi, npwx2, psi, npwx2, sr, nvecx )
     !
  END IF
  !
  CALL mp_sum( sr( :, 1:nbase ), intra_bgrp_comm )
#endif
  !
  IF ( lrot ) THEN
#ifdef USE_CUDA
!$cuf kernel do(1) <<<*,*>>>
     DO n = 1, nbase
        !
        e_d(n) = hr_d(n,n)
        vr_d(n,n) = 1.D0
        !
     END DO
#else
     DO n = 1, nbase
        !
        e(n) = hr(n,n)
        vr(n,n) = 1.D0
        !
     END DO
#endif
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !
#ifdef USE_CUDA
#ifdef CDIAG_CPU
  hr = hr_d
  sr = sr_d
  IF( my_bgrp_id == root_bgrp_id ) THEN
    CALL rdiaghg( nbase, nvec, hr, sr, nvecx, ew, vr )
  END IF
  IF( nbgrp > 1 ) THEN
    CALL mp_bcast( vr, root_bgrp_id, inter_bgrp_comm )
    CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
  ENDIF
  ew_d = ew
  vr_d = vr
#else
     IF( my_bgrp_id == root_bgrp_id ) THEN
        CALL rdiaghg( nbase, nvec, hr, sr, nvecx, ew, vr )
     END IF
     IF( nbgrp > 1 ) THEN
        CALL mp_bcast( vr, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
     ENDIF
#endif
!$cuf kernel do(1) <<<*,*>>>
     DO i = 1, nvec
        e_d(i) = ew_d(i)
     END DO
#else
     IF( my_bgrp_id == root_bgrp_id ) THEN
        CALL rdiaghg( nbase, nvec, hr, sr, nvecx, ew, vr )
     END IF
     IF( nbgrp > 1 ) THEN
        CALL mp_bcast( vr, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
     ENDIF
     !
     e(1:nvec) = ew(1:nvec)
#endif
     !
  END IF
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
     dav_iter = kter
     !
     CALL start_clock( 'regterg:update' )
     !
     np = 0
     !
#ifdef USE_CUDA
! Need further help to TODO:

#else
     DO n = 1, nvec
        !
        IF ( .NOT. conv(n) ) THEN
           !
           ! ... this root not yet converged ... 
           !
           np = np + 1
           !
           ! ... reorder eigenvectors so that coefficients for unconverged
           ! ... roots come first. This allows to use quick matrix-matrix 
           ! ... multiplications to set a new basis vector (see below)
           !
           IF ( np /= n ) THEN
             vr(:,np) = vr(:,n)
           ENDIF
           !
           ! ... for use in g_psi
           !
           ew(nbase+np) = e(n)
           !   
        END IF
        !
     END DO
#endif
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
#ifdef USE_CUDA
     IF ( uspp ) THEN
        !
        CALL cublasDgemm( 'N', 'N', npw2, notcnv, nbase, 1.D0, &
                    spsi_d, npwx2, vr_d, nvecx, 0.D0, psi_d(1,nb1), npwx2 )
        !
     ELSE
        !
        CALL cublasDgemm( 'N', 'N', npw2, notcnv, nbase, 1.D0, &
                    psi_d, npwx2, vr_d, nvecx, 0.D0, psi_d(1,nb1), npwx2 )
        !
     END IF

#else
     IF ( uspp ) THEN
        !
        CALL DGEMM( 'N', 'N', npw2, notcnv, nbase, 1.D0, &
                    spsi, npwx2, vr, nvecx, 0.D0, psi(1,nb1), npwx2 )
        !
     ELSE
        !
        CALL DGEMM( 'N', 'N', npw2, notcnv, nbase, 1.D0, &
                    psi, npwx2, vr, nvecx, 0.D0, psi(1,nb1), npwx2 )
        !
     END IF
#endif
     !
#ifdef USE_CUDA
!$cuf kernel do(2) <<<*,*>>>
    Do i=1,notcnv
          Do k=1,npwx
            psi_d(k,nbase+i) = - ew_d(nbase+i)*psi_d(k,nbase+i) 
          end do
    end do
#else
     DO np = 1, notcnv
        !
        psi(:,nbase+np) = - ew(nbase+np) * psi(:,nbase+np)
        !
     END DO
#endif
     !
#ifdef USE_CUDA
     CALL cublasDgemm( 'N', 'N', npw2, notcnv, nbase, 1.D0, &
                 hpsi_d, npwx2, vr_d, nvecx, 1.D0, psi_d(1,nb1), npwx2 )
#else
     CALL DGEMM( 'N', 'N', npw2, notcnv, nbase, 1.D0, &
                 hpsi, npwx2, vr, nvecx, 1.D0, psi(1,nb1), npwx2 )
#endif
     !
     CALL stop_clock( 'regterg:update' )
     !
     ! ... approximate inverse iteration
     !
#ifdef USE_CUDA
     CALL g_psi_gpu( npwx, npw, notcnv, 1, psi_d(1,nb1), ew_d(nb1) )
#else
     CALL g_psi( npwx, npw, notcnv, 1, psi(1,nb1), ew(nb1) )
#endif
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in 
     ! ... order to improve numerical stability of subspace diagonalization 
     ! ... (rdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     DO n = 1, notcnv
        !
        ew(n) = 2.D0 * ddot( npw2, psi(1,nbase+n), 1, psi(1,nbase+n), 1 )
        !
        IF ( gstart == 2 ) ew(n) = ew(n) - psi(1,nbase+n) * psi(1,nbase+n)
        !
     END DO
     !
     CALL mp_sum( ew( 1:notcnv ), intra_bgrp_comm )
     !
     DO n = 1, notcnv
        !
        psi(:,nbase+n) = psi(:,nbase+n) / SQRT( ew(n) )
        ! ... set Im[ psi(G=0) ] -  needed for numerical stability
        IF ( gstart == 2 ) psi(1,nbase+n) = CMPLX( DBLE(psi(1,nbase+n)), 0.D0 ,kind=DP)
        !
     END DO
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     CALL h_psi( npwx, npw, notcnv, psi(1,nb1), hpsi(1,nb1) )
     !
     IF ( uspp ) CALL s_psi( npwx, npw, notcnv, psi(1,nb1), spsi(1,nb1) )
     !
     ! ... update the reduced hamiltonian
     !
     CALL start_clock( 'regterg:overlap' )
     !
     CALL DGEMM( 'T', 'N', nbase+notcnv, notcnv, npw2, 2.D0, psi, &
                 npwx2, hpsi(1,nb1), npwx2, 0.D0, hr(1,nb1), nvecx )
     !
     IF ( gstart == 2 ) &
        CALL DGER( nbase+notcnv, notcnv, -1.D0, psi, &
                   npwx2, hpsi(1,nb1), npwx2, hr(1,nb1), nvecx )
     !
     CALL mp_sum( hr( :, nb1 : nb1+notcnv-1 ), intra_bgrp_comm )
     !
     IF ( uspp ) THEN
        !
        CALL DGEMM( 'T', 'N', nbase+notcnv, notcnv, npw2, 2.D0, psi, &
                    npwx2, spsi(1,nb1), npwx2, 0.D0, sr(1,nb1), nvecx )
        !
        IF ( gstart == 2 ) &
           CALL DGER( nbase+notcnv, notcnv, -1.D0, psi, &
                      npwx2, spsi(1,nb1), npwx2, sr(1,nb1), nvecx )
        !
     ELSE
        !
        CALL DGEMM( 'T', 'N', nbase+notcnv, notcnv, npw2, 2.D0, psi, &
                    npwx2, psi(1,nb1), npwx2, 0.D0, sr(1,nb1) , nvecx )
        !
        IF ( gstart == 2 ) &
           CALL DGER( nbase+notcnv, notcnv, -1.D0, psi, &
                      npwx2, psi(1,nb1), npwx2, sr(1,nb1), nvecx )
        !
     END IF
     !
     CALL mp_sum( sr( :, nb1 : nb1+notcnv-1 ), intra_bgrp_comm  )
     !
     CALL stop_clock( 'regterg:overlap' )
     !
     nbase = nbase + notcnv
     !
     DO n = 1, nbase
        !
        DO m = n + 1, nbase
           !
           hr(m,n) = hr(n,m)
           sr(m,n) = sr(n,m)
           !
        END DO
        !
     END DO
     !
     ! ... diagonalize the reduced hamiltonian
     !
     IF( my_bgrp_id == root_bgrp_id ) THEN
        CALL rdiaghg( nbase, nvec, hr, sr, nvecx, ew, vr )
     END IF
     IF( nbgrp > 1 ) THEN
        CALL mp_bcast( vr, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
     ENDIF
     !
     ! ... test for convergence
     !
     WHERE( btype(1:nvec) == 1 )
        !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < ethr ) )
        !
     ELSEWHERE
        !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < empty_ethr ) )
        !
     END WHERE
     ! ... next line useful for band parallelization of exact exchange
     IF ( nbgrp > 1 ) CALL mp_bcast(conv,root_bgrp_id,inter_bgrp_comm)
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
     e(1:nvec) = ew(1:nvec)
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     IF ( notcnv == 0 .OR. &
          nbase+notcnv > nvecx .OR. dav_iter == maxter ) THEN
        !
        CALL start_clock( 'regterg:last' )
        !
        CALL DGEMM( 'N', 'N', npw2, nvec, nbase, 1.D0, &
                    psi, npwx2, vr, nvecx, 0.D0, evc, npwx2 )
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'regterg:last' )
           !
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           WRITE( stdout, '(5X,"WARNING: ",I5, &
                &   " eigenvalues not converged in regterg")' ) notcnv
           !
           CALL stop_clock( 'regterg:last' )
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        psi(:,1:nvec) = evc(:,1:nvec)
        !
        IF ( uspp ) THEN
           !
           CALL DGEMM( 'N', 'N', npw2, nvec, nbase, 1.D0, spsi, &
                       npwx2, vr, nvecx, 0.D0, psi(1,nvec+1), npwx2 )
           !
           spsi(:,1:nvec) = psi(:,nvec+1:nvec+nvec)
           !
        END IF
        !
        CALL DGEMM( 'N', 'N', npw2, nvec, nbase, 1.D0, hpsi, &
                    npwx2, vr, nvecx, 0.D0, psi(1,nvec+1), npwx2 )
        !
        hpsi(:,1:nvec) = psi(:,nvec+1:nvec+nvec)
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        hr(:,1:nbase) = 0.D0
        sr(:,1:nbase) = 0.D0
        vr(:,1:nbase) = 0.D0
        !
        DO n = 1, nbase
           !
           hr(n,n) = e(n)
           sr(n,n) = 1.D0
           vr(n,n) = 1.D0
           !
        END DO
        !
        CALL stop_clock( 'regterg:last' )
        !
     END IF
     !
  END DO iterate
  !
  DEALLOCATE( conv )
  DEALLOCATE( ew )
  DEALLOCATE( vr )
  DEALLOCATE( hr )
  DEALLOCATE( sr )
  !
  IF ( uspp ) DEALLOCATE( spsi )
  !
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )  
  !
  CALL stop_clock( 'regterg' )
  !
  RETURN
  !
END SUBROUTINE regterg
!
!
!  Subroutine with distributed matrixes
!  (written by Carlo Cavazzoni)
!
!----------------------------------------------------------------------------
SUBROUTINE pregterg( npw, npwx, nvec, nvecx, evc, ethr, &
                    uspp, gstart, e, btype, notcnv, lrot, dav_iter )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an uspp matrix, evc is a complex vector
  ! ... (real wavefunctions with only half plane waves stored)
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE mp_bands,  ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id
  USE mp_diag,   ONLY : ortho_comm, np_ortho, me_ortho, ortho_comm_id, leg_ortho, &
                        ortho_parent_comm, ortho_cntx
  USE descriptors,      ONLY : la_descriptor, descla_init, descla_local_dims
  USE parallel_toolkit, ONLY : dsqmdst, dsqmcll, dsqmred, dsqmsym
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum
  USE cpu_gpu_interface
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx, gstart
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set
    !    (the basis set is refreshed when its dimension would exceed nvecx)
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx,nvec)
    !  evc   contains the  refined estimates of the eigenvectors
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence: root improvement is stopped,
    ! when two consecutive estimates of the root differ by less than ethr.
  LOGICAL, INTENT(IN) :: uspp
    ! if .FALSE. : S|psi> not needed
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  LOGICAL, INTENT(IN) :: lrot
    ! .TRUE. if the wfc have already be rotated
  REAL(DP), INTENT(OUT) :: e(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer  number of iterations performed
    ! number of unconverged roots
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, n, m, nb1
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
  INTEGER :: ierr
  REAL(DP), ALLOCATABLE :: ew(:)
  REAL(DP), ALLOCATABLE :: hl(:,:), sl(:,:), vl(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! eigenvectors of the Hamiltonian
    ! eigenvalues of the reduced hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi(:,:), hpsi(:,:), spsi(:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  LOGICAL, ALLOCATABLE :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr 
    ! threshold for empty bands
  INTEGER :: npw2, npwx2
  TYPE(la_descriptor) :: desc, desc_old
  INTEGER, ALLOCATABLE :: irc_ip( : )
  INTEGER, ALLOCATABLE :: nrc_ip( : )
  INTEGER, ALLOCATABLE :: rank_ip( :, : )
    ! matrix distribution descriptors
  INTEGER :: nx
    ! maximum local block dimension
  LOGICAL :: la_proc
    ! flag to distinguish procs involved in linear algebra
  INTEGER, ALLOCATABLE :: notcnv_ip( : )
  INTEGER, ALLOCATABLE :: ic_notcnv( : )
  !
  REAL(DP), EXTERNAL :: ddot
  !
  ! EXTERNAL  h_psi, s_psi, g_psi
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi> 
    ! s_psi(npwx,npw,nvec,psi,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,nvec)
    ! g_psi(npwx,npw,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
  !
  !
  CALL start_clock( 'regterg' )
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'regter', 'nvecx is too small', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  ALLOCATE( psi(  npwx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate psi ', ABS(ierr) )
  ALLOCATE( hpsi( npwx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate hpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     ALLOCATE( spsi( npwx, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( 'pregterg ',' cannot allocate spsi ', ABS(ierr) )
  END IF
  !
  ! ... Initialize the matrix descriptor
  !
  ALLOCATE( ic_notcnv( np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate ic_notcnv ', ABS(ierr) )
  ALLOCATE( notcnv_ip( np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate notcnv_ip ', ABS(ierr) )
  ALLOCATE( irc_ip( np_ortho(1) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate irc_ip ', ABS(ierr) )
  ALLOCATE( nrc_ip( np_ortho(1) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate nrc_ip ', ABS(ierr) )
  ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate rank_ip ', ABS(ierr) )
  !
  CALL desc_init( nvec, desc, irc_ip, nrc_ip  )
  !
  IF( la_proc ) THEN
     !
     ! only procs involved in the diagonalization need to allocate local 
     ! matrix block.
     !
     ALLOCATE( vl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( 'pregterg ',' cannot allocate vl ', ABS(ierr) )
     ALLOCATE( sl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( 'pregterg ',' cannot allocate sl ', ABS(ierr) )
     ALLOCATE( hl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( 'pregterg ',' cannot allocate hl ', ABS(ierr) )
     !
  ELSE
     !
     ALLOCATE( vl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( 'pregterg ',' cannot allocate vl ', ABS(ierr) )
     ALLOCATE( sl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( 'pregterg ',' cannot allocate sl ', ABS(ierr) )
     ALLOCATE( hl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( 'pregterg ',' cannot allocate hl ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate ew ', ABS(ierr) )
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( 'pregterg ',' cannot allocate conv ', ABS(ierr) )
  !
  npw2  = 2*npw
  npwx2  = 2*npwx
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
  IF ( uspp ) spsi = ZERO
  !
  hpsi = ZERO
  psi  = ZERO
  psi(:,1:nvec) = evc(:,1:nvec)
  ! ... set Im[ psi(G=0) ] -  needed for numerical stability
  IF ( gstart == 2 ) psi(1,1:nvec) = CMPLX( DBLE( psi(1,1:nvec) ), 0.D0 ,kind=DP)
  !
  ! ... hpsi contains h times the basis vectors
  !
  CALL h_psi( npwx, npw, nvec, psi, hpsi )
  !
  IF ( uspp ) CALL s_psi( npwx, npw, nvec, psi, spsi )
  !
  ! ... hl contains the projection of the hamiltonian onto the reduced
  ! ... space, vl contains the eigenvectors of hl. Remember hl, vl and sl
  ! ... are all distributed across processors, global replicated matrixes
  ! ... here are never allocated
  !
  CALL compute_distmat( hl, psi, hpsi )
  !
  IF ( uspp ) THEN
     !
     CALL compute_distmat( sl, psi, spsi )
     !
  ELSE
     !
     CALL compute_distmat( sl, psi, psi )
     !
  END IF
  !
  IF ( lrot ) THEN
     !
     CALL set_e_from_h()
     !
     CALL set_to_identity( vl, desc )
     !
  ELSE
     !
     ! ... diagonalize the reduced hamiltonian
     !     Calling block parallel algorithm
     !
     IF( my_bgrp_id == root_bgrp_id ) THEN
        CALL prdiaghg( nbase, hl, sl, nx, ew, vl, desc )
     END IF
     IF( nbgrp > 1 ) THEN
        CALL mp_bcast( vl, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
     ENDIF
     !
     e(1:nvec) = ew(1:nvec)
     !
  END IF

  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter
     !
     dav_iter = kter
     !
     CALL start_clock( 'regterg:update' )
     !
     CALL reorder_v()
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL hpsi_dot_v()
     !
     CALL stop_clock( 'regterg:update' )
     !
     ! ... approximate inverse iteration
     !
     CALL g_psi( npwx, npw, notcnv, 1, psi(1,nb1), ew(nb1) )
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in 
     ! ... order to improve numerical stability of subspace diagonalization 
     ! ... (cdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     DO n = 1, notcnv
        !
        ew(n) = 2.D0 * ddot( npw2, psi(1,nbase+n), 1, psi(1,nbase+n), 1 )
        !
        IF ( gstart == 2 ) ew(n) = ew(n) - psi(1,nbase+n) * psi(1,nbase+n)
        !
     END DO
     !
     CALL mp_sum( ew( 1:notcnv ), intra_bgrp_comm )
     !
     DO n = 1, notcnv
        !
        psi(:,nbase+n) = psi(:,nbase+n) / SQRT( ew(n) )
        ! ... set Im[ psi(G=0) ] -  needed for numerical stability
        IF ( gstart == 2 ) psi(1,nbase+n) = CMPLX( DBLE(psi(1,nbase+n)), 0.D0 ,kind=DP)
        !
     END DO
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     CALL h_psi( npwx, npw, notcnv, psi(1,nb1), hpsi(1,nb1) )
     !
     IF ( uspp ) CALL s_psi( npwx, npw, notcnv, psi(1,nb1), spsi(1,nb1) )
     !
     ! ... update the reduced hamiltonian
     !
     ! we need to save the old descriptor in order to redistribute matrices 
     !
     desc_old = desc
     !
     ! ... RE-Initialize the matrix descriptor
     !
     CALL desc_init( nbase+notcnv, desc, irc_ip, nrc_ip  )
     !
     IF( la_proc ) THEN

        !  redistribute hl and sl (see dsqmred), since the dimension of the subspace has changed
        !
        vl = hl
        DEALLOCATE( hl )
        ALLOCATE( hl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( 'pregterg ',' cannot allocate hl ', ABS(ierr) )
        CALL dsqmred( nbase, vl, desc_old%nrcx, desc_old, nbase+notcnv, hl, nx, desc )

        vl = sl
        DEALLOCATE( sl )
        ALLOCATE( sl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( 'pregterg ',' cannot allocate sl ', ABS(ierr) )
        CALL dsqmred( nbase, vl, desc_old%nrcx, desc_old, nbase+notcnv, sl, nx, desc )

        DEALLOCATE( vl )
        ALLOCATE( vl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( 'pregterg ',' cannot allocate vl ', ABS(ierr) )

     END IF
     !
     CALL start_clock( 'regterg:overlap' )
     !
     CALL update_distmat( hl, psi, hpsi )
     !
     IF ( uspp ) THEN
        !
        CALL update_distmat( sl, psi, spsi )
        !
     ELSE
        !
        CALL update_distmat( sl, psi, psi )
        !
     END IF
     !
     CALL stop_clock( 'regterg:overlap' )
     !
     nbase = nbase + notcnv
     !
     ! ... diagonalize the reduced hamiltonian
     !     Call block parallel algorithm
     !
     IF( my_bgrp_id == root_bgrp_id ) THEN
        CALL prdiaghg( nbase, hl, sl, nx, ew, vl, desc )
     END IF
     IF( nbgrp > 1 ) THEN
        CALL mp_bcast( vl, root_bgrp_id, inter_bgrp_comm )
        CALL mp_bcast( ew, root_bgrp_id, inter_bgrp_comm )
     ENDIF
     !
     ! ... test for convergence
     !
     WHERE( btype(1:nvec) == 1 )
        !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < ethr ) )
        !
     ELSEWHERE
        !
        conv(1:nvec) = ( ( ABS( ew(1:nvec) - e(1:nvec) ) < empty_ethr ) )
        !
     END WHERE
     ! ... next line useful for band parallelization of exact exchange
     IF ( nbgrp > 1 ) CALL mp_bcast(conv,root_bgrp_id,inter_bgrp_comm)
     !
     notcnv = COUNT( .NOT. conv(:) )
     !
     e(1:nvec) = ew(1:nvec)
     !
     ! ... if overall convergence has been achieved, or the dimension of
     ! ... the reduced basis set is becoming too large, or in any case if
     ! ... we are at the last iteration refresh the basis set. i.e. replace
     ! ... the first nvec elements with the current estimate of the
     ! ... eigenvectors;  set the basis dimension to nvec.
     !
     IF ( notcnv == 0 .OR. nbase+notcnv > nvecx .OR. dav_iter == maxter ) THEN
        !
        CALL start_clock( 'regterg:last' )
        !
        CALL refresh_evc()       
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'regterg:last' )
           !
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           WRITE( stdout, '(5X,"WARNING: ",I5, &
                &   " eigenvalues not converged")' ) notcnv
           !
           CALL stop_clock( 'regterg:last' )
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        psi(:,1:nvec) = evc(:,1:nvec)
        !
        IF ( uspp ) THEN
           !
           CALL refresh_spsi()
           ! 
        END IF
        !
        CALL refresh_hpsi()
        !
        ! ... refresh the reduced hamiltonian
        !
        nbase = nvec
        !
        CALL desc_init( nvec, desc, irc_ip, nrc_ip  )
        !
        IF( la_proc ) THEN
           !
           ! note that nx has been changed by desc_init
           ! we need to re-alloc with the new size.
           !
           DEALLOCATE( vl, hl, sl )
           ALLOCATE( vl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( 'pregterg ',' cannot allocate vl ', ABS(ierr) )
           ALLOCATE( hl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( 'pregterg ',' cannot allocate hl ', ABS(ierr) )
           ALLOCATE( sl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( 'pregterg ',' cannot allocate sl ', ABS(ierr) )
           !
        END IF
        !
        CALL set_h_from_e( )
        !
        CALL set_to_identity( vl, desc )
        CALL set_to_identity( sl, desc )
        !
        CALL stop_clock( 'regterg:last' )
        !
     END IF
     !
  END DO iterate
  !
  DEALLOCATE( vl, hl, sl )
  !
  DEALLOCATE( rank_ip )
  DEALLOCATE( ic_notcnv )
  DEALLOCATE( irc_ip )
  DEALLOCATE( nrc_ip )
  DEALLOCATE( notcnv_ip )
  DEALLOCATE( conv )
  DEALLOCATE( ew )
  !
  IF ( uspp ) DEALLOCATE( spsi )
  !
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )  
  !
  CALL stop_clock( 'regterg' )
  !
  RETURN
  !
  !
CONTAINS
  !
  !
  SUBROUTINE desc_init( nsiz, desc, irc_ip, nrc_ip )
     !
     INTEGER, INTENT(IN)  :: nsiz
     TYPE(la_descriptor), INTENT(OUT) :: desc
     INTEGER, INTENT(OUT) :: irc_ip(:)
     INTEGER, INTENT(OUT) :: nrc_ip(:)

     INTEGER :: i, j, rank
     !
     CALL descla_init( desc, nsiz, nsiz, np_ortho, me_ortho, ortho_comm, ortho_cntx, ortho_comm_id )
     !
     nx = desc%nrcx
     !
     DO j = 0, desc%npc - 1
        CALL descla_local_dims( irc_ip( j + 1 ), nrc_ip( j + 1 ), desc%n, desc%nx, np_ortho(1), j )
        DO i = 0, desc%npr - 1
           CALL GRID2D_RANK( 'R', desc%npr, desc%npc, i, j, rank )
           rank_ip( i+1, j+1 ) = rank * leg_ortho
        END DO
     END DO
     !
     la_proc = .FALSE.
     IF( desc%active_node > 0 ) la_proc = .TRUE.
     !
     RETURN
  END SUBROUTINE desc_init
  !
  !
  SUBROUTINE set_to_identity( distmat, desc )
     TYPE(la_descriptor), INTENT(IN)  :: desc
     REAL(DP), INTENT(OUT) :: distmat(:,:)
     INTEGER :: i
     distmat = 0_DP
     IF( desc%myc == desc%myr .AND. desc%active_node > 0 ) THEN
        DO i = 1, desc%nc
           distmat( i, i ) = 1_DP
        END DO
     END IF 
     RETURN
  END SUBROUTINE set_to_identity
  !
  !
  SUBROUTINE reorder_v()
     !
     INTEGER :: ipc, ipr
     INTEGER :: nc, ic
     INTEGER :: nl, npl
     !
     np = 0
     !
     notcnv_ip = 0
     !
     n = 0
     !
     DO ipc = 1, desc%npc
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        npl = 0
        !
        IF( ic <= nvec ) THEN
           !
           DO nl = 1, min( nvec - ic + 1, nc )
              !
              n  = n  + 1
              !
              IF ( .NOT. conv(n) ) THEN
                 !
                 ! ... this root not yet converged ... 
                 !
                 np  = np  + 1
                 npl = npl + 1
                 IF( npl == 1 ) ic_notcnv( ipc ) = np
                 !
                 ! ... reorder eigenvectors so that coefficients for unconverged
                 ! ... roots come first. This allows to use quick matrix-matrix 
                 ! ... multiplications to set a new basis vector (see below)
                 !
                 notcnv_ip( ipc ) = notcnv_ip( ipc ) + 1
                 !
                 IF ( npl /= nl ) THEN
                    IF( la_proc .AND. desc%myc == ipc-1 ) THEN
                       vl( :, npl) = vl( :, nl )
                    END IF
                 END IF
                 !
                 ! ... for use in g_psi
                 !
                 ew(nbase+np) = e(n)
                 !   
              END IF
              !
           END DO
           !
        END IF
        !
     END DO
     !
  END SUBROUTINE reorder_v
  !
  !
  SUBROUTINE hpsi_dot_v()
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, notcl, root, np
     REAL(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP), ALLOCATABLE :: ptmp( :, : )
     REAL(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     ALLOCATE( ptmp( npwx, nx ) )

     DO ipc = 1, desc%npc
        !
        IF( notcnv_ip( ipc ) > 0 ) THEN

           notcl = notcnv_ip( ipc )
           ic    = ic_notcnv( ipc ) 

           ptmp = 0.0d0
           beta = 0.0d0

           DO ipr = 1, desc%npr
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == desc%myr .AND. ipc-1 == desc%myc .AND. la_proc ) THEN
                 vtmp(:,1:notcl) = vl(:,1:notcl)
              END IF

              CALL mp_bcast( vtmp(:,1:notcl), root, ortho_parent_comm )
              ! 
              IF ( uspp ) THEN
                 !
                 CALL DGEMM( 'N', 'N', npw2, notcl, nr, 1.D0, &
                    spsi( 1, ir ), npwx2, vtmp, nx, beta, psi(1,nb1+ic-1), npwx2 )
                 !
              ELSE
                 !
                 CALL DGEMM( 'N', 'N', npw2, notcl, nr, 1.D0, &
                    psi( 1, ir ), npwx2, vtmp, nx, beta, psi(1,nb1+ic-1), npwx2 )
                 !
              END IF
              !
              CALL DGEMM( 'N', 'N', npw2, notcl, nr, 1.D0, &
                      hpsi( 1, ir ), npwx2, vtmp, nx, 1.D0, ptmp, npwx2 )

              beta = 1.0d0

           END DO

           DO np = 1, notcl
              !
              psi(:,nbase+np+ic-1) = ptmp(:,np) - ew(nbase+np+ic-1) * psi(:,nbase+np+ic-1)
              !
           END DO
           !
        END IF
        !
     END DO

     
     DEALLOCATE( vtmp )
     DEALLOCATE( ptmp )

     RETURN
  END SUBROUTINE hpsi_dot_v
  !
  !
  SUBROUTINE refresh_evc( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), ALLOCATABLE :: vtmp( :, : )
     REAL(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, desc%npc
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = 0.0d0

           DO ipr = 1, desc%npr
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == desc%myr .AND. ipc-1 == desc%myc .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 ! 
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          psi(1,ir), npwx2, vl, nx, beta, evc(1,ic), npwx2 )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          psi(1,ir), npwx2, vtmp, nx, beta, evc(1,ic), npwx2 )
              END IF
              ! 

              beta = 1.0d0

           END DO
           !
        END IF
        !
     END DO
     !
     DEALLOCATE( vtmp )

     RETURN
  END SUBROUTINE refresh_evc
  !
  !
  SUBROUTINE refresh_spsi( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), ALLOCATABLE :: vtmp( :, : )
     REAL(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, desc%npc
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = 0_DP
           !
           DO ipr = 1, desc%npr
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == desc%myr .AND. ipc-1 == desc%myc .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 ! 
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          spsi(1,ir), npwx2, vl, nx, beta, psi(1,nvec+ic), npwx2 )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          spsi(1,ir), npwx2, vtmp, nx, beta, psi(1,nvec+ic), npwx2 )
              END IF
              ! 
              beta = 1_DP

           END DO
           !
        END IF
        !
     END DO
     !
     spsi(:,1:nvec) = psi(:,nvec+1:nvec+nvec)
     !
     DEALLOCATE( vtmp )

     RETURN
  END SUBROUTINE refresh_spsi
  !
  !
  !
  SUBROUTINE refresh_hpsi( )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), ALLOCATABLE :: vtmp( :, : )
     REAL(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     !
     DO ipc = 1, desc%npc
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic <= nvec ) THEN
           !
           nc = min( nc, nvec - ic + 1 )
           !
           beta = 0.0d0
           !
           DO ipr = 1, desc%npr
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              IF( ipr-1 == desc%myr .AND. ipc-1 == desc%myc .AND. la_proc ) THEN
                 !
                 !  this proc sends his block
                 ! 
                 CALL mp_bcast( vl(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          hpsi(1,ir), npwx2, vl, nx, beta, psi(1,nvec+ic), npwx2 )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL DGEMM( 'N', 'N', npw2, nc, nr, 1.D0, &
                          hpsi(1,ir), npwx2, vtmp, nx, beta, psi(1,nvec+ic), npwx2 )
              END IF
              ! 
              beta = 1.0d0

           END DO
           !
        END IF
        !
     END DO
     !
     DEALLOCATE( vtmp )

     hpsi(:,1:nvec) = psi(:,nvec+1:nvec+nvec)

     RETURN
  END SUBROUTINE refresh_hpsi
  !
  !

  SUBROUTINE compute_distmat( dm, v, w )
     !
     !  This subroutine compute <vi|wj> and store the
     !  result in distributed matrix dm 
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root
     REAL(DP), INTENT(OUT) :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     REAL(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = 0.0d0
     !
     DO ipc = 1, desc%npc !  loop on column procs 
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        DO ipr = 1, ipc ! use symmetry for the loop on row procs
           !
           nr = nrc_ip( ipr )
           ir = irc_ip( ipr )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL DGEMM( 'T', 'N', nr, nc, npw2, 2.D0 , &
                       v(1,ir), npwx2, w(1,ic), npwx2, 0.D0, work, nx )

           IF ( gstart == 2 ) &
              CALL DGER( nr, nc, -1.D0, v(1,ir), npwx2, w(1,ic), npwx2, work, nx )

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, ortho_parent_comm )

        END DO
        !
     END DO
     IF (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) dm = dm/nbgrp
     !
     CALL dsqmsym( nbase, dm, nx, desc )
     !
     DEALLOCATE( work )
     !
     RETURN
  END SUBROUTINE compute_distmat
  !
  !
  SUBROUTINE update_distmat( dm, v, w )
     !
     INTEGER :: ipc, ipr
     INTEGER :: nr, nc, ir, ic, root, icc, ii
     REAL(DP)    :: dm( :, : )
     COMPLEX(DP) :: v(:,:), w(:,:)
     REAL(DP), ALLOCATABLE :: vtmp( :, : )

     ALLOCATE( vtmp( nx, nx ) )
     !
     vtmp = 0.0d0
     !
     DO ipc = 1, desc%npc
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        IF( ic+nc-1 >= nb1 ) THEN

           nc = MIN( nc, ic+nc-1 - nb1 + 1 )
           IF( ic >= nb1 ) THEN
              ii = ic
              icc = 1
           ELSE
              ii = nb1
              icc = nb1-ic+1
           END IF

           DO ipr = 1, ipc ! desc%npr use symmetry
              !
              nr = nrc_ip( ipr )
              ir = irc_ip( ipr )
              !
              root = rank_ip( ipr, ipc )

              CALL DGEMM( 'T', 'N', nr, nc, npw2, 2.D0, v( 1, ir ), &
                          npwx2, w(1,ii), npwx2, 0.D0, vtmp, nx )
              !
              IF ( gstart == 2 ) &
                 CALL DGER( nr, nc, -1.D0, v( 1, ir ), npwx2, w(1,ii), npwx2, vtmp, nx )
              IF (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) vtmp = vtmp/nbgrp

              IF(  (desc%active_node > 0) .AND. (ipr-1 == desc%myr) .AND. (ipc-1 == desc%myc) ) THEN
                 CALL mp_root_sum( vtmp(:,1:nc), dm(:,icc:icc+nc-1), root, ortho_parent_comm )
              ELSE
                 CALL mp_root_sum( vtmp(:,1:nc), dm, root, ortho_parent_comm )
              END IF


           END DO
           !
        END IF
        !
     END DO
     !
     CALL dsqmsym( nbase+notcnv, dm, nx, desc )
     !
     DEALLOCATE( vtmp )
     RETURN
  END SUBROUTINE update_distmat
  !
  !
  !
  SUBROUTINE set_e_from_h()
     INTEGER :: nc, ic, i
     e(1:nbase) = 0.0d0
     IF( desc%myc == desc%myr .AND. la_proc ) THEN
        nc = desc%nc
        ic = desc%ic
        DO i = 1, nc
           e( i + ic - 1 ) = hl( i, i )
        END DO
     END IF
     CALL mp_sum( e(1:nbase), ortho_parent_comm )
     RETURN
  END SUBROUTINE set_e_from_h
  !
  SUBROUTINE set_h_from_e()
     INTEGER :: nc, ic, i
     IF( la_proc ) THEN
        hl = 0.0d0
        IF( desc%myc == desc%myr ) THEN
           nc = desc%nc
           ic = desc%ic
           DO i = 1, nc
              hl(i,i) = e( i + ic - 1 )
           END DO
        END IF
     END IF
     RETURN
  END SUBROUTINE set_h_from_e
  !
END SUBROUTINE pregterg
