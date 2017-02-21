!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

!#define CDIAG_CPU

#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
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
!
!----------------------------------------------------------------------------
SUBROUTINE cegterg( npw, npwx, nvec, nvecx, npol, evc, evc_d, ethr, &
                    uspp, e, e_d, btype, notcnv, lrot, dav_iter )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an overlap matrix, evc is a complex vector
  !
  USE kinds,         ONLY : DP
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
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set :
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! umber of spin polarizations
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx,npol,nvec)
    !  evc contains the  refined estimates of the eigenvectors  
#ifdef USE_CUDA
  COMPLEX(DP), DEVICE, INTENT(INOUT) :: evc_d(npwx,npol,nvec)
  REAL(DP), DEVICE, INTENT(OUT) :: e_d(nvec)
#else
  COMPLEX(DP), OPTIONAL :: evc_d(:)
#endif
    !  evc_d contains the refined estimates of the eigenvectors  
  REAL(DP), INTENT(IN) :: ethr
    ! energy threshold for convergence :
    !   root improvement is stopped, when two consecutive estimates of the root
    !   differ by less than ethr.
  LOGICAL, INTENT(IN) :: uspp
    ! if .FALSE. : do not calculate S|psi>
  INTEGER, INTENT(IN) :: btype(nvec)
    ! band type ( 1 = occupied, 0 = empty )
  LOGICAL, INTENT(IN) :: lrot
    ! .TRUE. if the wfc have already been rotated
  REAL(DP), INTENT(OUT) :: e(nvec)
    ! contains the estimated roots.
  INTEGER, INTENT(OUT) :: dav_iter, notcnv
    ! integer number of iterations performed
    ! number of unconverged roots
  !
  ! ... LOCAL variables
  !
  INTEGER, PARAMETER :: maxter = 20
    ! maximum number of iterations
  !
  INTEGER :: kter, nbase, np, kdim, kdmx, n, m, nb1, nbn
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! adapted npw and npwx
    ! do-loop counters
  INTEGER :: ierr
  COMPLEX(DP), ALLOCATABLE :: hc(:,:), sc(:,:), vc(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! the eigenvectors of the Hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi(:,:,:), hpsi(:,:,:), spsi(:,:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  REAL(DP), ALLOCATABLE :: ew(:)
    ! eigenvalues of the reduced hamiltonian
  LOGICAL, ALLOCATABLE  :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr
    ! threshold for empty bands
  !
#ifndef USE_CUDA
  REAL(DP), EXTERNAL :: ddot
#endif

#ifdef USE_CUDA
  attributes(pinned) :: vc
  COMPLEX(DP), DEVICE, ALLOCATABLE :: hc_d(:,:), sc_d(:,:), vc_d(:,:), vc_temp_d(:,:)
  COMPLEX(DP), PINNED, ALLOCATABLE :: comm_h_c(:,:)
  COMPLEX(DP), DEVICE, ALLOCATABLE :: psi_d(:,:,:), hpsi_d(:,:,:), spsi_d(:,:,:)
  REAL(DP), DEVICE, ALLOCATABLE :: ew_d(:)
  REAL(DP), PINNED, ALLOCATABLE :: comm_h_r(:)
  LOGICAL, DEVICE, ALLOCATABLE :: conv_d(:)
  INTEGER :: conv_idx(nvec)
  INTEGER, DEVICE :: conv_idx_d(nvec)
  INTEGER :: istat, i, j, k
  INTEGER(KIND=8) :: freeMem,totalMem
  REAL(DP) :: rFreeMem,rUsedMem,rTotalMem,minUsedMem,maxUsedMem,minFreeMem,maxFreeMem
  REAL(DP) :: ew_temp1, ew_temp2
!  REAL(DP) :: timer !, it_timer
#endif

  !
  ! EXTERNAL  h_psi,    s_psi,    g_psi
    ! h_psi(npwx,npw,nvec,psi,hpsi)
    !     calculates H|psi>
    ! s_psi(npwx,npw,nvec,spsi)
    !     calculates S|psi> (if needed)
    !     Vectors psi,hpsi,spsi are dimensioned (npwx,npol,nvec)
    ! g_psi(npwx,npw,notcnv,psi,e)
    !    calculates (diag(h)-e)^-1 * psi, diagonal approx. to (h-e)^-1*psi
    !    the first nvec columns contain the trial eigenvectors
!  print *,"cegterg: ",npw, npwx, nvec, nvecx, npol
!  timer = MPI_Wtime()
  !
!  CALL read_compare( evc_d, "evc" )
  !
  CALL start_clock( 'cegterg' )
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'cegterg', 'nvecx is too small', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx*npol
     kdmx = npwx*npol
     !
  END IF
  !
  ALLOCATE(  psi( npwx, npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate psi ', ABS(ierr) )
  ALLOCATE( hpsi( npwx, npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate hpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     ALLOCATE( spsi( npwx, npol, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' cegterg ',' cannot allocate spsi ', ABS(ierr) )
#ifdef USE_CUDA
     ALLOCATE( spsi_d( npwx, npol, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' cegterg ',' cannot allocate spsi_d ', ABS(ierr) )
#endif
  END IF
  !
  ALLOCATE( sc( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate sc ', ABS(ierr) )
  ALLOCATE( hc( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate hc ', ABS(ierr) )
  ALLOCATE( vc( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate vc ', ABS(ierr) )
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate ew ', ABS(ierr) )
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate conv ', ABS(ierr) )
  !
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.

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

  ALLOCATE( sc_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate sc_d ', ABS(ierr) )
  ALLOCATE( hc_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate hc_d ', ABS(ierr) )
  ALLOCATE( vc_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate vc_d ', ABS(ierr) )

  ALLOCATE( vc_temp_d( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate vc_temp_d ', ABS(ierr) )

  ALLOCATE( comm_h_c( nvecx, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate comm_h_c ', ABS(ierr) )

#if 0
      istat=CudaMemGetInfo(freeMem,totalMem)
      rTotalMem = totalMem/(10.**6)
      rFreeMem = freeMem/(10.**6);            MinFreeMem = rFreeMem; MaxFreeMem = rFreeMem
      rUsedMem = (totalMem-freeMem)/(10.**6); MaxUsedMem = rUsedMem; MinUsedMem = rUsedMem
write(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " GPU memory used: ",minUsedMem," - ",maxUsedMem," / ",rTotalMem," MBytes"
write(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " GPU memory free: ",minFreeMem," - ",maxFreeMem," / ",rTotalMem," MBytes"
print *," "
call flush(6)
print *,"psi_d: ",npwx,npol,nvecx
call flush(6)
#endif

  ALLOCATE(  psi_d( npwx, npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' cegterg ',' cannot allocate psi_d ', ABS(ierr) )

#if 0
      istat=CudaMemGetInfo(freeMem,totalMem)
      rTotalMem = totalMem/(10.**6)
      rFreeMem = freeMem/(10.**6);            MinFreeMem = rFreeMem; MaxFreeMem = rFreeMem
      rUsedMem = (totalMem-freeMem)/(10.**6); MaxUsedMem = rUsedMem; MinUsedMem = rUsedMem
write(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " GPU memory used: ",minUsedMem," - ",maxUsedMem," / ",rTotalMem," MBytes"
write(*,"(A20,F7.1,A3,F7.1,A3,F7.1,A8)") " GPU memory free: ",minFreeMem," - ",maxFreeMem," / ",rTotalMem," MBytes"
print *," "
call flush(6)
#endif

  ALLOCATE( hpsi_d( npwx, npol, nvecx ), STAT=ierr )
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
 !$cuf kernel do(3) <<<*,*>>>
  Do k=1,nvec
    Do j=lbound(psi_d,2), ubound(psi_d,2)
      Do i=lbound(psi_d,1), ubound(psi_d,1)
        psi_d(i,j,k) = evc_d(i,j,k)
      end do
    end do
  end do

  CALL h_psi( npwx, npw, nvec, psi_d, hpsi_d )
#else

  hpsi = ZERO
  psi  = ZERO
  psi(:,:,1:nvec) = evc(:,:,1:nvec)

  CALL h_psi( npwx, npw, nvec, psi, hpsi )

#endif

!  CALL read_compare( hpsi_d, "hpsi" )

  !
  ! ... hpsi contains h times the basis vectors
  !
  !CALL h_psi( npwx, npw, nvec, psi, hpsi )
  !hpsi_d = hpsi
  !    Do k=1,nvec
  !      Call compare( hpsi(:,:,k), hpsi_d(:,:,k),"hpsi")
  !    end do
  !hpsi = hpsi_d
  
  !
  ! ... spsi contains s times the basis vectors
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
  ! ... hc contains the projection of the hamiltonian onto the reduced 
  ! ... space vc contains the eigenvectors of hc
  !
#ifdef USE_CUDA
  hc_d(:,:) = ZERO
  sc_d(:,:) = ZERO
  vc_d(:,:) = ZERO
  ew_d(:) = 0.d0
#else
  hc(:,:) = ZERO
  sc(:,:) = ZERO
  vc(:,:) = ZERO
  ew(:) = 0.d0
#endif
  !
#ifdef USE_CUDA
  CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
              psi_d, kdmx, hpsi_d, kdmx, ZERO, hc_d, nvecx )

  comm_h_c( :, 1:nbase ) = hc_d( :, 1:nbase )
  CALL mp_sum( comm_h_c( :, 1:nbase ), intra_bgrp_comm )
  hc_d( :, 1:nbase ) = comm_h_c( :, 1:nbase )
#else
  CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
              psi, kdmx, hpsi, kdmx, ZERO, hc, nvecx )

  CALL mp_sum( hc( :, 1:nbase ), intra_bgrp_comm )

#endif

#ifdef USE_CUDA
  IF ( uspp ) THEN
     !
     CALL cublasZgemm( 'C', 'N', nbase, nbase, kdim, ONE, &
                 psi_d, kdmx, spsi_d, kdmx, ZERO, sc_d, nvecx )
     !     
  ELSE
     !
     CALL cublasZgemm( 'C', 'N', nbase, nbase, kdim, ONE, &
                 psi_d, kdmx, psi_d, kdmx, ZERO, sc_d, nvecx )
     !
  END IF

 comm_h_c(:,1:nbase) = sc_d(:,1:nbase)
 CALL mp_sum( comm_h_c( :, 1:nbase ), intra_bgrp_comm )
 sc_d(:,1:nbase) = comm_h_c(:,1:nbase)
#else
  !
  IF ( uspp ) THEN
     !
     CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
                 psi, kdmx, spsi, kdmx, ZERO, sc, nvecx )
     !     
  ELSE
     !
     CALL ZGEMM( 'C', 'N', nbase, nbase, kdim, ONE, &
                 psi, kdmx, psi, kdmx, ZERO, sc, nvecx )
     !
  END IF

 CALL mp_sum( sc( :, 1:nbase ), intra_bgrp_comm )

#endif
  !
  !
  IF ( lrot ) THEN
!   print *,"LROT !!!!!!!!!!!!!!!!!!!"
     !
#ifdef USE_CUDA
!$cuf kernel do(1) <<<*,*>>>
     DO i = 1, nbase
        !
        e_d(i) = REAL( hc_d(i,i) )
        !
        vc_d(i,i) = ONE
        !
     END DO
#else

     DO n = 1, nbase
        !
        e(n) = REAL( hc(n,n) )
        !
        vc(n,n) = ONE
        !
     END DO

#endif
     !
  ELSE
!  print *,"NO LROT!!"
     !
     ! ... diagonalize the reduced hamiltonian
     !
#ifdef USE_CUDA
#ifdef CDIAG_CPU

     hc = hc_d
     sc = sc_d
     CALL cdiaghg( nbase, nvec, hc, sc, nvecx, ew, vc )
     ew_d = ew
     vc_d = vc
#else
     CALL cdiaghg( nbase, nvec, hc_d, sc_d, nvecx, ew_d, vc_d )
#endif

!$cuf kernel do(1) <<<*,*>>>
     DO i = 1, nvec
        e_d(i) = ew_d(i)
     END DO

#else

     CALL cdiaghg( nbase, nvec, hc, sc, nvecx, ew, vc )
     !
     e(1:nvec) = ew(1:nvec)

#endif
     !
  END IF
  !
  ! ... iterate
  !
  iterate: DO kter = 1, maxter

!   print *,"starting kter = ",kter 
     !
     dav_iter = kter
     !
     !it_timer = MPI_Wtime()
     CALL start_clock( 'cegterg:update' )
     !
     np = 0
     !
#ifdef USE_CUDA
     DO n = 1, nvec
        conv_idx(n) = -1
        IF ( .NOT. conv(n) ) THEN
           np = np + 1
           conv_idx(n) = np
        END IF
     END DO

    conv_idx_d = conv_idx
!vc_temp_d = vc_d

!$cuf kernel do(2) <<<*,*>>>
    Do j=1,nvec
      Do k=1,nvecx
        vc_temp_d(k,j) = vc_d(k,j)
      end do
    end do


!$cuf kernel do(2) <<<*,*>>>
    Do j=1,nvec
       Do k=1,nvecx
          if(conv_idx_d(j) /= -1) then
            vc_d(k,conv_idx_d(j)) = vc_temp_d(k,j)
            if(k==1) ew_d(nbase+conv_idx_d(j)) = e_d(j)
          end if
       end do
    end do
#else

     DO n = 1, nvec
        !
        conv_idx(n) = -1
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
           conv_idx(n) = np

           IF ( np /= n ) then
              vc(:,np) = vc(:,n)
           END IF
           !
           ! ... for use in g_psi
           !
           ew(nbase+np) = e(n)
           !
        END IF
        !
     END DO

#endif
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
#ifdef USE_CUDA
     IF ( uspp ) THEN
        !
        CALL cublasZgemm( 'N', 'N', kdim, notcnv, nbase, ONE, spsi_d, &
                    kdmx, vc_d, nvecx, ZERO, psi_d(1,1,nb1), kdmx )
        !     
     ELSE
        !
        CALL cublasZgemm( 'N', 'N', kdim, notcnv, nbase, ONE, psi_d, &
                    kdmx, vc_d, nvecx, ZERO, psi_d(1,1,nb1), kdmx )
        !
     END IF

#else
     IF ( uspp ) THEN
        !
        CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, spsi, &
                    kdmx, vc, nvecx, ZERO, psi(1,1,nb1), kdmx )
        !     
     ELSE
        !
        CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, psi, &
                    kdmx, vc, nvecx, ZERO, psi(1,1,nb1), kdmx )
        !
     END IF
#endif
     !
#ifdef USE_CUDA
!$cuf kernel do(3) <<<*,*>>>
    Do i=1,notcnv
       Do j=1,npol
          Do k=1,npwx
            psi_d(k,j,nbase+i) = - ew_d(nbase+i)*psi_d(k,j,nbase+i) 
          end do
       end do
    end do

#else
     DO np = 1, notcnv
        !
        psi(:,:,nbase+np) = - ew(nbase+np)*psi(:,:,nbase+np)
        !
     END DO
#endif

#ifdef USE_CUDA
     CALL cublasZgemm( 'N', 'N', kdim, notcnv, nbase, ONE, hpsi_d, &
                 kdmx, vc_d, nvecx, ONE, psi_d(1,1,nb1), kdmx )
#else
     CALL ZGEMM( 'N', 'N', kdim, notcnv, nbase, ONE, hpsi, &
                 kdmx, vc, nvecx, ONE, psi(1,1,nb1), kdmx )
#endif
     !
     CALL stop_clock( 'cegterg:update' )
     !
     ! ... approximate inverse iteration
     !
#ifdef USE_CUDA
     CALL g_psi_gpu( npwx, npw, notcnv, npol, psi_d(1,1,nb1), ew_d(nb1) )
#else

     CALL g_psi( npwx, npw, notcnv, npol, psi(1,1,nb1), ew(nb1) )
#endif

#ifdef USE_CUDA
     DO n = 1, notcnv
        !
        nbn = nbase + n
        !
        IF ( npol == 1 ) THEN
           !
#if 0
!$cuf kernel do(1) <<<*,*>>>
           do i = 1, npw
             ew_d(n) = ew_d(n) +  REAL(psi_d(i,1,nbn)) *  REAL(psi_d(i,1,nbn)) &
                               + AIMAG(psi_d(i,1,nbn)) * AIMAG(psi_d(i,1,nbn))
           end do
#else
           call myDdot( 2*npw, psi_d(1,1,nbn), ew_temp1 )
           !call myDdot( 2*npw, psi_d(1,2,nbn), ew_temp2 )
           ew(n) = ew_temp1 !+ ew_temp2
           !print *,"ew ",ew(n),ew_temp1
#endif
           !
        ELSE
           !
#if 0
!$cuf kernel do(1) <<<*,*>>>
           do i = 1, npw
             ew_d(n) = ew_d(n) +  REAL(psi_d(i,1,nbn)) *  REAL(psi_d(i,1,nbn)) &
                               + AIMAG(psi_d(i,1,nbn)) * AIMAG(psi_d(i,1,nbn)) &
                               +  REAL(psi_d(i,2,nbn)) *  REAL(psi_d(i,2,nbn)) &
                               + AIMAG(psi_d(i,2,nbn)) * AIMAG(psi_d(i,2,nbn))

           end do
#else
           call myDdot( 2*npw, psi_d(1,1,nbn), ew_temp1 )
           call myDdot( 2*npw, psi_d(1,2,nbn), ew_temp2 )
           ew(n) = ew_temp1 + ew_temp2
#endif
           !
        END IF
        !
     END DO
     !ew_d(1:notcnv) = ew(1:notcnv)
     !comm_h_r(1:notcnv) = ew_d(1:notcnv)
     !comm_h_r(1:notcnv) = ew(1:notcnv)
     !
     !CALL mp_sum( comm_h_r( 1:notcnv ), intra_bgrp_comm )
     !
     !ew(1:notcnv) = comm_h_r(1:notcnv)
     CALL mp_sum( ew( 1:notcnv ), intra_bgrp_comm )
     ew_d(1:notcnv) = ew(1:notcnv)
#else
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in
     ! ... order to improve numerical stability of subspace diagonalization
     ! ... (cdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     DO n = 1, notcnv
        !
        nbn = nbase + n
        !
        IF ( npol == 1 ) THEN
           !
           ew(n) = ddot( 2*npw, psi(1,1,nbn), 1, psi(1,1,nbn), 1 )
           !
        ELSE
           !
           ew(n) = ddot( 2*npw, psi(1,1,nbn), 1, psi(1,1,nbn), 1 ) + &
                   ddot( 2*npw, psi(1,2,nbn), 1, psi(1,2,nbn), 1 )
           !
        END IF
        !
     END DO
     !
     CALL mp_sum( ew( 1:notcnv ), intra_bgrp_comm )
     !
#endif

#ifdef USE_CUDA

!ew_d(1:notcnv) = ew(1:notcnv)

!$cuf kernel do(3) <<<*,*>>>
    Do i=1,notcnv
       Do j=1,npol
          Do k=1,npwx
            psi_d(k,j,nbase+i) = psi_d(k,j,nbase+i)/SQRT( ew_d(i) )
          end do
       end do
    end do
#else
     DO n = 1, notcnv
        !
        psi(:,:,nbase+n) = psi(:,:,nbase+n) / SQRT( ew(n) )
        !
     END DO
#endif

     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     !
#ifdef USE_CUDA

  CALL h_psi( npwx, npw, notcnv, psi_d(:,:,nb1), hpsi_d(:,:,nb1) )
#else
     CALL h_psi( npwx, npw, notcnv, psi(1,1,nb1), hpsi(1,1,nb1) )
#endif
     !
  IF ( uspp ) THEN
#ifdef USE_CUDA
    CALL s_psi( npwx, npw, notcnv, psi_d(1,1,nb1), spsi_d(1,1,nb1) )
#else
    CALL s_psi( npwx, npw, notcnv, psi(1,1,nb1), spsi(1,1,nb1) )
#endif
 END IF
     !
     ! ... update the reduced hamiltonian
     !
     CALL start_clock( 'cegterg:overlap' )
     !
#ifdef USE_CUDA
     CALL cublasZgemm( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi_d, &
                 kdmx, hpsi_d(1,1,nb1), kdmx, ZERO, hc_d(1,nb1), nvecx )

      comm_h_c( :, nb1:nb1+notcnv-1 ) =  hc_d( :, nb1:nb1+notcnv-1 )
     !
     CALL mp_sum( comm_h_c( :, nb1:nb1+notcnv-1 ), intra_bgrp_comm )
     !
      hc_d( :, nb1:nb1+notcnv-1 ) =  comm_h_c( :, nb1:nb1+notcnv-1 )
#else
     CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi, &
                 kdmx, hpsi(1,1,nb1), kdmx, ZERO, hc(1,nb1), nvecx )

     CALL mp_sum( hc( :, nb1:nb1+notcnv-1 ), intra_bgrp_comm )
#endif

#ifdef USE_CUDA
     IF ( uspp ) THEN
        !
        CALL cublasZgemm( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi_d, &
                    kdmx, spsi_d(1,1,nb1), kdmx, ZERO, sc_d(1,nb1), nvecx )
        !     
     ELSE
        !
        CALL cublasZgemm( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi_d, &
                    kdmx, psi_d(1,1,nb1), kdmx, ZERO, sc_d(1,nb1), nvecx )
        !
     END IF

     comm_h_c( :, nb1:nb1+notcnv-1 ) = sc_d( :, nb1:nb1+notcnv-1 )
     !
     CALL mp_sum( comm_h_c( :, nb1:nb1+notcnv-1 ), intra_bgrp_comm )
     !
     sc_d( :, nb1:nb1+notcnv-1 ) = comm_h_c( :, nb1:nb1+notcnv-1 )
#else
     IF ( uspp ) THEN
        !
        CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi, &
                    kdmx, spsi(1,1,nb1), kdmx, ZERO, sc(1,nb1), nvecx )
        !     
     ELSE
        !
        CALL ZGEMM( 'C', 'N', nbase+notcnv, notcnv, kdim, ONE, psi, &
                    kdmx, psi(1,1,nb1), kdmx, ZERO, sc(1,nb1), nvecx )
        !
     END IF
     CALL mp_sum( sc( :, nb1:nb1+notcnv-1 ), intra_bgrp_comm )
#endif
     CALL stop_clock( 'cegterg:overlap' )
     !

     nbase = nbase + notcnv
     !
#ifdef USE_CUDA
!$cuf kernel do(1) <<<*,*>>>
    Do i=1,nbase
        hc_d(i,i) = CMPLX( REAL( hc_d(i,i) ), 0.D0 ,kind=DP)
        sc_d(i,i) = CMPLX( REAL( sc_d(i,i) ), 0.D0 ,kind=DP)
        Do j=i+1,nbase
           hc_d(j,i) = CONJG( hc_d(i,j) )
           sc_d(j,i) = CONJG( sc_d(i,j) )
        end do
    end do
#else
     DO n = 1, nbase
        !
        ! ... the diagonal of hc and sc must be strictly real 
        !
        hc(n,n) = CMPLX( REAL( hc(n,n) ), 0.D0 ,kind=DP)
        sc(n,n) = CMPLX( REAL( sc(n,n) ), 0.D0 ,kind=DP)
        !
        DO m = n + 1, nbase
           !
           hc(m,n) = CONJG( hc(n,m) )
           sc(m,n) = CONJG( sc(n,m) )
           !
        END DO
        !
     END DO
#endif

     !
     ! ... diagonalize the reduced hamiltonian
     !
#ifdef USE_CUDA

#ifdef CDIAG_CPU

     hc = hc_d
     sc = sc_d
     CALL cdiaghg( nbase, nvec, hc, sc, nvecx, ew, vc )
     ew_d = ew
     vc_d = vc
#else
     CALL cdiaghg( nbase, nvec, hc_d, sc_d, nvecx, ew_d, vc_d )
#endif
#else
     CALL cdiaghg( nbase, nvec, hc, sc, nvecx, ew, vc )
#endif
#ifdef USE_CUDA
     ew = ew_d
     e = e_d
#endif
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

! TODO: Do we need convergence test on GPU? LEave off for now
!#ifdef USE_CUDA
#if 0
    !$cuf kernel do(1) <<<*,*>>>
    Do i=1,nvec
      if( btype_d(i) == 1 ) then
        conv_d(i) = ( ABS( ew_d(i) - e_d(i) ) < ethr )
      else
        conv_d(i) = ( ABS( ew_d(i) - e_d(i) ) < empty_ethr )
      end if
    end do

   ! call compare(conv, conv_d, "conv test")
#endif
     ! ... next line useful for band parallelization of exact exchange
     IF ( nbgrp > 1 ) then
       print *,"NBGRP > 1 !!!!!"
       CALL mp_bcast(conv,root_bgrp_id,inter_bgrp_comm)
     END IF
     !
     notcnv = COUNT( .NOT. conv(:) )
     !it_timer = MPI_Wtime() - it_timer

     print *,"kter: ",kter," notcnv: ",notcnv," nbase: ",nbase !,it_timer,"sec"
 
     !
#ifdef USE_CUDA
    !$cuf kernel do(1) <<<*,*>>>
    Do i=1,nvec
      e_d(i) = ew_d(i)
    end do
#else
     e(1:nvec) = ew(1:nvec)
#endif

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
        CALL start_clock( 'cegterg:last' )
        !
#ifdef USE_CUDA
        CALL cublasZgemm( 'N', 'N', kdim, nvec, nbase, ONE, &
                    psi_d, kdmx, vc_d, nvecx, ZERO, evc_d, kdmx )
#else
        CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, &
                    psi, kdmx, vc, nvecx, ZERO, evc, kdmx )
#endif
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'cegterg:last' )
           !
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           !!!WRITE( stdout, '(5X,"WARNING: ",I5, &
           !!!     &   " eigenvalues not converged")' ) notcnv
           !
           CALL stop_clock( 'cegterg:last' )
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
    Do i=1,nvec
      Do j=lbound(psi_d,2),ubound(psi_d,2)
        Do k=lbound(psi_d,1),ubound(psi_d,1)
          psi_d(k,j,i) = evc_d(k,j,i)
        end do
      end do
    end do
#else
        psi(:,:,1:nvec) = evc(:,:,1:nvec)
#endif
        !
        IF ( uspp ) THEN
           !
#ifdef USE_CUDA
           CALL cublasZgemm( 'N', 'N', kdim, nvec, nbase, ONE, spsi_d, &
                       kdmx, vc_d, nvecx, ZERO, psi_d(1,1,nvec+1), kdmx )
#else
           CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, spsi, &
                       kdmx, vc, nvecx, ZERO, psi(1,1,nvec+1), kdmx )
#endif

#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
    Do i=1,nvec
      Do j=lbound(psi_d,2),ubound(psi_d,2)
        Do k=lbound(psi_d,1),ubound(psi_d,1)
          spsi_d(k,j,i) = psi_d(k,j,i+nvec)
        end do
      end do
    end do
#else
           spsi(:,:,1:nvec) = psi(:,:,nvec+1:nvec+nvec)
#endif
           !
        END IF
        !
#ifdef USE_CUDA
        CALL cublasZgemm( 'N', 'N', kdim, nvec, nbase, ONE, hpsi_d, &
                    kdmx, vc_d, nvecx, ZERO, psi_d(1,1,nvec+1), kdmx )
#else
        CALL ZGEMM( 'N', 'N', kdim, nvec, nbase, ONE, hpsi, &
                    kdmx, vc, nvecx, ZERO, psi(1,1,nvec+1), kdmx )
#endif
        !
#ifdef USE_CUDA
    !$cuf kernel do(3) <<<*,*>>>
    Do i=1,nvec
      Do j=lbound(psi_d,2),ubound(psi_d,2)
        Do k=lbound(psi_d,1),ubound(psi_d,1)
          hpsi_d(k,j,i) = psi_d(k,j,i+nvec)
        end do
      end do
    end do
#else
        hpsi(:,:,1:nvec) = psi(:,:,nvec+1:nvec+nvec)
#endif

        !
        ! ... refresh the reduced hamiltonian 
        !
        nbase = nvec
        !
#ifdef USE_CUDA
        hc_d(:,1:nbase) = ZERO
        sc_d(:,1:nbase) = ZERO
        vc_d(:,1:nbase) = ZERO
#else
        hc(:,1:nbase) = ZERO
        sc(:,1:nbase) = ZERO
        vc(:,1:nbase) = ZERO
#endif
        !
#ifdef USE_CUDA
    !$cuf kernel do(1) <<<*,*>>>
        DO i = 1, nbase
           !
!           hc(n,n) = REAL( e(n) )
           hc_d(i,i) = CMPLX( e_d(i), 0.0_DP ,kind=DP)
           !
           sc_d(i,i) = ONE
           vc_d(i,i) = ONE
           !
        END DO
#else
        DO n = 1, nbase
           !
!           hc(n,n) = REAL( e(n) )
           hc(n,n) = CMPLX( e(n), 0.0_DP ,kind=DP)
           !
           sc(n,n) = ONE
           vc(n,n) = ONE
           !
        END DO
#endif
        !
        CALL stop_clock( 'cegterg:last' )
        !
     END IF
     !
  END DO iterate
  !

  !copy outputs to CPU (not sure if this is needed)
#ifdef USE_CUDA
  evc = evc_d
  e = e_d
#endif

  DEALLOCATE( conv )
  DEALLOCATE( ew )
  DEALLOCATE( vc )
  DEALLOCATE( hc )
  DEALLOCATE( sc )

#ifdef USE_CUDA
  DEALLOCATE( conv_d )
  DEALLOCATE( ew_d )
  DEALLOCATE( vc_d )
  DEALLOCATE( hc_d )
  DEALLOCATE( sc_d )
#endif
  !
  IF ( uspp ) then
   DEALLOCATE( spsi )
#ifdef USE_CUDA
   DEALLOCATE( spsi_d )
#endif
  end if
  !
  DEALLOCATE( hpsi )
  DEALLOCATE( psi )

#ifdef USE_CUDA
  DEALLOCATE( hpsi_d )
  DEALLOCATE( psi_d )
#endif
  !
!  timer = MPI_Wtime() - timer
!  print *,"CEGTERG in ",timer," sec"
  CALL stop_clock( 'cegterg' )
  !
!  CALL read_compare( e_d, "e_o" )
!  CALL read_compare( evc_d, "evc_o" )
  !
  RETURN
  !
END SUBROUTINE cegterg

!
!  Subroutine with distributed matrixes
!  (written by Carlo Cavazzoni)
!
!----------------------------------------------------------------------------
SUBROUTINE pcegterg( npw, npwx, nvec, nvecx, npol, evc, ethr, &
                    uspp, e, btype, notcnv, lrot, dav_iter )
  !----------------------------------------------------------------------------
  !
  ! ... iterative solution of the eigenvalue problem:
  !
  ! ... ( H - e S ) * evc = 0
  !
  ! ... where H is an hermitean operator, e is a real scalar,
  ! ... S is an uspp matrix, evc is a complex vector
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE mp_bands,  ONLY : intra_bgrp_comm, inter_bgrp_comm, root_bgrp_id, nbgrp, my_bgrp_id
  USE mp_diag,   ONLY : ortho_comm, np_ortho, me_ortho, ortho_comm_id, leg_ortho, &
                        ortho_parent_comm, ortho_cntx
  USE descriptors,      ONLY : la_descriptor, descla_init , descla_local_dims
  USE parallel_toolkit, ONLY : zsqmred, zsqmher, zsqmdst
  USE mp,               ONLY : mp_bcast, mp_root_sum, mp_sum, mp_barrier
  USE cpu_gpu_interface
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: npw, npwx, nvec, nvecx, npol
    ! dimension of the matrix to be diagonalized
    ! leading dimension of matrix evc, as declared in the calling pgm unit
    ! integer number of searched low-lying roots
    ! maximum dimension of the reduced basis set
    !    (the basis set is refreshed when its dimension would exceed nvecx)
    ! number of spin polarizations
  COMPLEX(DP), INTENT(INOUT) :: evc(npwx,npol,nvec)
    !  evc   contains the  refined estimates of the eigenvectors
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
  INTEGER :: kter, nbase, np, kdim, kdmx, n, nb1, nbn
    ! counter on iterations
    ! dimension of the reduced basis
    ! counter on the reduced basis vectors
    ! do-loop counters
  INTEGER :: ierr
  REAL(DP), ALLOCATABLE :: ew(:)
  COMPLEX(DP), ALLOCATABLE :: hl(:,:), sl(:,:), vl(:,:)
    ! Hamiltonian on the reduced basis
    ! S matrix on the reduced basis
    ! eigenvectors of the Hamiltonian
    ! eigenvalues of the reduced hamiltonian
  COMPLEX(DP), ALLOCATABLE :: psi(:,:,:), hpsi(:,:,:), spsi(:,:,:)
    ! work space, contains psi
    ! the product of H and psi
    ! the product of S and psi
  LOGICAL, ALLOCATABLE :: conv(:)
    ! true if the root is converged
  REAL(DP) :: empty_ethr 
    ! threshold for empty bands
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
  CALL start_clock( 'cegterg' )
  !
  IF ( nvec > nvecx / 2 ) CALL errore( 'pcegterg', 'nvecx is too small', 1 )
  !
  ! ... threshold for empty bands
  !
  empty_ethr = MAX( ( ethr * 5.D0 ), 1.D-5 )
  !
  IF ( npol == 1 ) THEN
     !
     kdim = npw
     kdmx = npwx
     !
  ELSE
     !
     kdim = npwx*npol
     kdmx = npwx*npol
     !
  END IF

  ALLOCATE(  psi( npwx, npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate psi ', ABS(ierr) )
  !
  ALLOCATE( hpsi( npwx, npol, nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate hpsi ', ABS(ierr) )
  !
  IF ( uspp ) THEN
     ALLOCATE( spsi( npwx, npol, nvecx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate spsi ', ABS(ierr) )
  END IF
  !
  ! ... Initialize the matrix descriptor
  !
  ALLOCATE( ic_notcnv( np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate ic_notcnv ', ABS(ierr) )
  !
  ALLOCATE( notcnv_ip( np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate notcnv_ip ', ABS(ierr) )
  !
  ALLOCATE( irc_ip( np_ortho(1) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate irc_ip ', ABS(ierr) )
  !
  ALLOCATE( nrc_ip( np_ortho(1) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate nrc_ip ', ABS(ierr) )
  !
  ALLOCATE( rank_ip( np_ortho(1), np_ortho(2) ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate rank_ip ', ABS(ierr) )
  !
  CALL desc_init( nvec, desc, irc_ip, nrc_ip )
  !
  IF( la_proc ) THEN
     !
     ! only procs involved in the diagonalization need to allocate local 
     ! matrix block.
     !
     ALLOCATE( vl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )
     !
     ALLOCATE( sl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )
     !
     ALLOCATE( hl( nx , nx ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )
     !
  ELSE
     !
     ALLOCATE( vl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )
     !
     ALLOCATE( sl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )
     !
     ALLOCATE( hl( 1 , 1 ), STAT=ierr )
     IF( ierr /= 0 ) &
        CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )
     !
  END IF
  !
  ALLOCATE( ew( nvecx ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate ew ', ABS(ierr) )
  !
  ALLOCATE( conv( nvec ), STAT=ierr )
  IF( ierr /= 0 ) &
     CALL errore( ' pcegterg ',' cannot allocate conv ', ABS(ierr) )
  !
  notcnv = nvec
  nbase  = nvec
  conv   = .FALSE.
  !
  IF ( uspp ) spsi = ZERO
  !
  hpsi = ZERO
  psi  = ZERO
  psi(:,:,1:nvec) = evc(:,:,1:nvec)
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
     CALL pcdiaghg( nbase, hl, sl, nx, ew, vl, desc )
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
     CALL start_clock( 'cegterg:update' )
     !
     CALL reorder_v()
     !
     nb1 = nbase + 1
     !
     ! ... expand the basis set with new basis vectors ( H - e*S )|psi> ...
     !
     CALL hpsi_dot_v()
     !
     CALL stop_clock( 'cegterg:update' )
     !
     ! ... approximate inverse iteration
     !
     CALL g_psi( npwx, npw, notcnv, npol, psi(1,1,nb1), ew(nb1) )
     !
     ! ... "normalize" correction vectors psi(:,nb1:nbase+notcnv) in 
     ! ... order to improve numerical stability of subspace diagonalization 
     ! ... (cdiaghg) ew is used as work array :
     !
     ! ...         ew = <psi_i|psi_i>,  i = nbase + 1, nbase + notcnv
     !
     DO n = 1, notcnv
        !
        nbn = nbase + n
        !
        IF ( npol == 1 ) THEN
           !
           ew(n) = ddot( 2*npw, psi(1,1,nbn), 1, psi(1,1,nbn), 1 )
           !
        ELSE
           !
           ew(n) = ddot( 2*npw, psi(1,1,nbn), 1, psi(1,1,nbn), 1 ) + &
                   ddot( 2*npw, psi(1,2,nbn), 1, psi(1,2,nbn), 1 )
           !
        END IF
        !
     END DO
     !
     CALL mp_sum( ew( 1:notcnv ), intra_bgrp_comm )
     !
     DO n = 1, notcnv
        !
        psi(:,:,nbase+n) = psi(:,:,nbase+n) / SQRT( ew(n) )
        !
     END DO
     !
     ! ... here compute the hpsi and spsi of the new functions
     !
     CALL h_psi( npwx, npw, notcnv, psi(1,1,nb1), hpsi(1,1,nb1) )
     !
     IF ( uspp ) &
        CALL s_psi( npwx, npw, notcnv, psi(1,1,nb1), spsi(1,1,nb1) )
     !
     ! ... update the reduced hamiltonian
     !
     ! we need to save the old descriptor in order to redistribute matrices 
     !
     desc_old = desc
     !
     ! ... RE-Initialize the matrix descriptor
     !
     CALL desc_init( nbase+notcnv, desc, irc_ip, nrc_ip )
     !
     IF( la_proc ) THEN

        !  redistribute hl and sl (see dsqmred), since the dimension of the subspace has changed
        !
        vl = hl
        DEALLOCATE( hl )
        ALLOCATE( hl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )

        CALL zsqmred( nbase, vl, desc_old%nrcx, desc_old, nbase+notcnv, hl, nx, desc )

        vl = sl
        DEALLOCATE( sl )
        ALLOCATE( sl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )

        CALL zsqmred( nbase, vl, desc_old%nrcx, desc_old, nbase+notcnv, sl, nx, desc )

        DEALLOCATE( vl )
        ALLOCATE( vl( nx , nx ), STAT=ierr )
        IF( ierr /= 0 ) &
           CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )

     END IF
     !
     CALL start_clock( 'cegterg:overlap' )
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
     CALL stop_clock( 'cegterg:overlap' )
     !
     nbase = nbase + notcnv
     !
     ! ... diagonalize the reduced hamiltonian
     !     Call block parallel algorithm
     !
     CALL pcdiaghg( nbase, hl, sl, nx, ew, vl, desc )
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
        CALL start_clock( 'cegterg:last' )
        !
        CALL refresh_evc()       
        !
        IF ( notcnv == 0 ) THEN
           !
           ! ... all roots converged: return
           !
           CALL stop_clock( 'cegterg:last' )
           !
           EXIT iterate
           !
        ELSE IF ( dav_iter == maxter ) THEN
           !
           ! ... last iteration, some roots not converged: return
           !
           !!!WRITE( stdout, '(5X,"WARNING: ",I5, &
           !!!     &   " eigenvalues not converged")' ) notcnv
           !
           CALL stop_clock( 'cegterg:last' )
           !
           EXIT iterate
           !
        END IF
        !
        ! ... refresh psi, H*psi and S*psi
        !
        psi(:,:,1:nvec) = evc(:,:,1:nvec)
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
        CALL desc_init( nvec, desc, irc_ip, nrc_ip )
        !
        IF( la_proc ) THEN
           !
           ! note that nx has been changed by desc_init
           ! we need to re-alloc with the new size.
           !
           DEALLOCATE( vl, hl, sl )
           ALLOCATE( vl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( ' pcegterg ',' cannot allocate vl ', ABS(ierr) )
           ALLOCATE( hl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( ' pcegterg ',' cannot allocate hl ', ABS(ierr) )
           ALLOCATE( sl( nx, nx ), STAT=ierr )
           IF( ierr /= 0 ) &
              CALL errore( ' pcegterg ',' cannot allocate sl ', ABS(ierr) )
           !
        END IF
        !
        CALL set_h_from_e( )
        !
        CALL set_to_identity( vl, desc )
        CALL set_to_identity( sl, desc )
        !
        CALL stop_clock( 'cegterg:last' )
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
  CALL stop_clock( 'cegterg' )
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
     CALL descla_init( desc, nsiz, nsiz, np_ortho, me_ortho, ortho_comm, ortho_cntx,  ortho_comm_id )
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
     COMPLEX(DP), INTENT(OUT) :: distmat(:,:)
     INTEGER :: i
     distmat = ( 0_DP , 0_DP )
     IF( desc%myc == desc%myr .AND. desc%active_node > 0 ) THEN
        DO i = 1, desc%nc
           distmat( i, i ) = ( 1_DP , 0_DP )
        END DO
     END IF 
     RETURN
  END SUBROUTINE set_to_identity
  !
  !
  SUBROUTINE reorder_v()
     !
     INTEGER :: ipc
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
     INTEGER :: nr, ir, ic, notcl, root, np
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP), ALLOCATABLE :: ptmp( :, :, : )
     COMPLEX(DP) :: beta

     ALLOCATE( vtmp( nx, nx ) )
     ALLOCATE( ptmp( npwx, npol, nx ) )

     DO ipc = 1, desc%npc
        !
        IF( notcnv_ip( ipc ) > 0 ) THEN

           notcl = notcnv_ip( ipc )
           ic    = ic_notcnv( ipc ) 

           ptmp = ZERO
           beta = ZERO

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
                 CALL ZGEMM( 'N', 'N', kdim, notcl, nr, ONE, &
                    spsi( 1, 1, ir ), kdmx, vtmp, nx, beta, psi(1,1,nb1+ic-1), kdmx )
                 !
              ELSE
                 !
                 CALL ZGEMM( 'N', 'N', kdim, notcl, nr, ONE, &
                    psi( 1, 1, ir ), kdmx, vtmp, nx, beta, psi(1,1,nb1+ic-1), kdmx )
                 !
              END IF
              !
              CALL ZGEMM( 'N', 'N', kdim, notcl, nr, ONE, &
                      hpsi( 1, 1, ir ), kdmx, vtmp, nx, ONE, ptmp, kdmx )

              beta = ONE

           END DO

           DO np = 1, notcl
              !
              psi(:,:,nbase+np+ic-1) = ptmp(:,:,np) - ew(nbase+np+ic-1) * psi(:,:,nbase+np+ic-1)
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
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

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
           beta = ZERO

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
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          psi(1,1,ir), kdmx, vl, nx, beta, evc(1,1,ic), kdmx )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          psi(1,1,ir), kdmx, vtmp, nx, beta, evc(1,1,ic), kdmx )
              END IF
              ! 

              beta = ONE

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
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

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
           beta = ZERO
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
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          spsi(1,1,ir), kdmx, vl, nx, beta, psi(1,1,nvec+ic), kdmx )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          spsi(1,1,ir), kdmx, vtmp, nx, beta, psi(1,1,nvec+ic), kdmx )
              END IF
              ! 
              beta = ONE

           END DO
           !
        END IF
        !
     END DO
     !
     spsi(:,:,1:nvec) = psi(:,:,nvec+1:nvec+nvec)
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
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )
     COMPLEX(DP) :: beta

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
           beta = ZERO
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
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          hpsi(1,1,ir), kdmx, vl, nx, beta, psi(1,1,nvec+ic), kdmx )
              ELSE
                 !
                 !  all other procs receive
                 ! 
                 CALL mp_bcast( vtmp(:,1:nc), root, ortho_parent_comm )
                 CALL ZGEMM( 'N', 'N', kdim, nc, nr, ONE, &
                          hpsi(1,1,ir), kdmx, vtmp, nx, beta, psi(1,1,nvec+ic), kdmx )
              END IF
              ! 
              beta = ONE

           END DO
           !
        END IF
        !
     END DO
     !
     DEALLOCATE( vtmp )

     hpsi(:,:,1:nvec) = psi(:,:,nvec+1:nvec+nvec)

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
     COMPLEX(DP), INTENT(OUT) :: dm( :, : )
     COMPLEX(DP) :: v(:,:,:), w(:,:,:)
     COMPLEX(DP), ALLOCATABLE :: work( :, : )
     !
     ALLOCATE( work( nx, nx ) )
     !
     work = ZERO
     !
     !  Only upper triangle is computed, then the matrix is hermitianized
     !
     DO ipc = 1, desc%npc !  loop on column procs 
        !
        nc = nrc_ip( ipc )
        ic = irc_ip( ipc )
        !
        DO ipr = 1, ipc ! desc%npr ! ipc ! use symmetry for the loop on row procs
           !
           nr = nrc_ip( ipr )
           ir = irc_ip( ipr )
           !
           !  rank of the processor for which this block (ipr,ipc) is destinated
           !
           root = rank_ip( ipr, ipc )

           ! use blas subs. on the matrix block

           CALL ZGEMM( 'C', 'N', nr, nc, kdim, ONE , &
                       v(1,1,ir), kdmx, w(1,1,ic), kdmx, ZERO, work, nx )

           ! accumulate result on dm of root proc.

           CALL mp_root_sum( work, dm, root, ortho_parent_comm )

        END DO
        !
     END DO
     if (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) dm = dm/nbgrp
     !
     !  The matrix is hermitianized using upper triangle
     !
     CALL zsqmher( nbase, dm, nx, desc )
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
     COMPLEX(DP) :: dm( :, : )
     COMPLEX(DP) :: v(:,:,:), w(:,:,:)
     COMPLEX(DP), ALLOCATABLE :: vtmp( :, : )

     ALLOCATE( vtmp( nx, nx ) )
     !
     vtmp = ZERO
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

              CALL ZGEMM( 'C', 'N', nr, nc, kdim, ONE, v( 1, 1, ir ), &
                          kdmx, w(1,1,ii), kdmx, ZERO, vtmp, nx )
              IF (ortho_parent_comm.ne.intra_bgrp_comm .and. nbgrp > 1) vtmp = vtmp/nbgrp
              !
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
     CALL zsqmher( nbase+notcnv, dm, nx, desc )
     !
     DEALLOCATE( vtmp )
     RETURN
  END SUBROUTINE update_distmat
  !
  !
  !
  SUBROUTINE set_e_from_h()
     INTEGER :: nc, ic, i
     e(1:nbase) = 0_DP
     IF( desc%myc == desc%myr .AND. la_proc ) THEN
        nc = desc%nc
        ic = desc%ic
        DO i = 1, nc
           e( i + ic - 1 ) = REAL( hl( i, i ) )
        END DO
     END IF
     CALL mp_sum( e(1:nbase), ortho_parent_comm )
     RETURN
  END SUBROUTINE set_e_from_h
  !
  SUBROUTINE set_h_from_e()
     INTEGER :: nc, ic, i
     IF( la_proc ) THEN
        hl = ZERO
        IF( desc%myc == desc%myr ) THEN
           nc = desc%nc
           ic = desc%ic
           DO i = 1, nc
              hl(i,i) = CMPLX( e( i + ic - 1 ), 0_DP ,kind=DP)
           END DO
        END IF
     END IF
     RETURN
  END SUBROUTINE set_h_from_e
  !
END SUBROUTINE pcegterg
