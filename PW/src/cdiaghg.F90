!
! Copyright (C) 2001-2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
#define ZERO ( 0.D0, 0.D0 )
#define ONE  ( 1.D0, 0.D0 )
!
#ifdef USE_GPU
#define MY_ROUTINE(x)  x##_gpu
#else
#define MY_ROUTINE(x)  x##_cpu
!MODULE ONLY COMPILED ONCE ( IN CPU PATH )
module zhegvx_module
#ifdef USE_CUDA
  use cudafor
#endif
  use kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  INTEGER :: first_time_zhegvx=1
!  INTEGER :: nb, lwork, lrwork, liwork
  INTEGER,     ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP),    ALLOCATABLE :: rwork(:), sdiag(:), hdiag(:)
  COMPLEX(DP), ALLOCATABLE :: work(:)
  COMPLEX(DP), ALLOCATABLE :: h_temp(:,:),s_temp(:,:)
#ifdef USE_CUDA
  ATTRIBUTES( PINNED ) :: iwork, ifail, rwork, sdiag, hdiag, work
  REAL(DP), ALLOCATABLE, PINNED :: e_h(:)
  COMPLEX(DP), ALLOCATABLE, PINNED :: v_h(:,:)
  REAL(DP),    DEVICE, ALLOCATABLE :: rwork_d(:), sdiag_d(:), hdiag_d(:)
  COMPLEX(DP), DEVICE, ALLOCATABLE :: work_d(:)
#endif
END MODULE zhegvx_module
#endif 

!----------------------------------------------------------------------------
SUBROUTINE MY_ROUTINE( cdiaghg )( n, m, h, s, ldh, e, v )
  !----------------------------------------------------------------------------
  !
  ! ... calculates eigenvalues and eigenvectors of the generalized problem
  ! ... Hv=eSv, with H hermitean matrix, S overlap matrix.
  ! ... On output both matrix are unchanged
  !
  ! ... LAPACK version - uses both ZHEGV and ZHEGVX
  !
  USE kinds,            ONLY : DP
  USE mp,               ONLY : mp_bcast, mp_sum, mp_barrier, mp_max
  USE mp_bands,         ONLY : me_bgrp, root_bgrp, intra_bgrp_comm
  USE zhegvx_module,    ONLY : iwork, ifail, rwork, work, first_time_zhegvx
#ifdef USE_GPU
  USE cudafor
  USE gpu_routines, ONLY: cufMemcpy2D, copy_diag, restore_upper_tri
  USE zhegvdx_gpu
  USE zhegvx_module,    ONLY : h_temp, s_temp, hdiag,hdiag_d, sdiag,sdiag_d
  USE zhegvx_module,    ONLY : e_h, v_h, rwork_d, work_d
#else
  USE zhegvx_module,    ONLY : h_temp, s_temp, hdiag, sdiag
#endif
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n, m, ldh
    ! dimension of the matrix to be diagonalized
    ! number of eigenstates to be calculate
    ! leading dimension of h, as declared in the calling pgm unit
  COMPLEX(DP), INTENT(INOUT) :: h(ldh,n), s(ldh,n)
    ! actually intent(in) but compilers don't know and complain
    ! matrix to be diagonalized
    ! overlap matrix
  REAL(DP), INTENT(OUT) :: e(n)
    ! eigenvalues
  COMPLEX(DP), INTENT(OUT) :: v(ldh,m)
    ! eigenvectors (column-wise)
#ifdef USE_GPU
  ATTRIBUTES( DEVICE ) :: h, s, e, v
#endif
  !
  INTEGER                  :: lwork, lrwork, liwork, nb, mm, info, i, j
  INTEGER                  :: max_lwork, max_lrwork, max_liwork
#ifdef USE_GPU
  INTEGER                  :: lwork_d, lrwork_d
  INTEGER                  :: max_lwork_d, max_lrwork_d
  INTEGER                  :: cuf_i, ii, istat, cpu_path
  REAL(DP)                 :: rFreeMem,rNeedMem
  INTEGER(KIND=8)          :: freeMem,totalMem,needMem
#endif
    ! mm = number of calculated eigenvectors
  REAL(DP)                 :: abstol
  LOGICAL                  :: all_eigenvalues
  INTEGER,  EXTERNAL       :: ILAENV
    ! ILAENV returns optimal block size "nb"
  !
  CALL start_clock( 'cdiaghg' )

#ifdef USE_GPU
  IF(first_time_zhegvx) THEN
     ALLOCATE( e_h(ldh), v_h(ldh,ldh) )
  ELSE IF( size( e_h ) < ldh .or. size( v_h ) < ldh*ldh ) THEN
     DEALLOCATE( e_h, v_h )
     ALLOCATE( e_h(ldh), v_h(ldh,ldh) )
  ENDIF
#endif  
  !
  ! ... only the first processor diagonalizes the matrix
  !
  IF ( me_bgrp == root_bgrp ) THEN
     !
     ! ... save the diagonal of input S (it will be overwritten)
     !
#ifdef USE_GPU
  ! Host workspace
  lwork = n
  lrwork = 1+5*n+2*n*n
  liwork = 3+5*n

  max_lwork = ldh
  max_lrwork = 1 + 5*ldh + 2*ldh*ldh
  max_liwork = 3 + 5*ldh

  ! Device workspace
  lwork_d  = 2*64*64 + 65*n
  lrwork_d = n
  max_lwork_d  = 2*64*64 + 65*ldh
  max_lrwork_d = ldh

  cpu_path=0

#if 0
  istat=CudaMemGetInfo(freeMem,totalMem)
  rFreeMem = freeMem/(10.**6)
  needMem = 8*2*ldh + 16*2*ldh*ldh + 16*max_lwork_d + 8*max_lrwork_d
  if(allocated(sdiag_d)) needMem = needMem - 8*size(sdiag_d)
  if(allocated(hdiag_d)) needMem = needMem - 8*size(hdiag_d)
  if(allocated(work_d)) needMem = needMem - 16*size(work_d)
  if(allocated(rwork_d)) needMem = needMem - 8*size(rwork_d)
  rNeedMem = needMem/(10.**6)

  if(rNeedMem > rFreeMem) then
    write(*,"(A16,F8.1,A16,F8.1,A26)") "ZHEGVX: GPU has",rFreeMem,"MB available / ",rNeedMem,"MB required --> using CPU"
    cpu_path=1
  endif
#endif
  
  if(cpu_path==1) then
    if(allocated(hdiag_d)) deallocate(hdiag_d)
    if(allocated(sdiag_d)) deallocate(sdiag_d)
    if(allocated(work_d)) deallocate(work_d)
    if(allocated(rwork_d)) deallocate(rwork_d)
#endif
     nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
     lwork = max(2*n, (nb+1)*n)
     lrwork = 7*n
     liwork = 5*n
     nb = max(nb, ILAENV( 1, 'ZHETRD', 'U', ldh, -1, -1, -1 ))
     max_lwork = 2*max(2*ldh, (nb+1)*ldh)
     max_lrwork = 2*7*ldh
     max_liwork = 2*5*ldh
#ifdef USE_GPU
  endif
#endif

  IF( ALLOCATED( work   ) .and. SIZE( work   ) < lwork  ) DEALLOCATE( work   )
  IF( ALLOCATED( rwork  ) .and. SIZE( rwork  ) < lrwork ) DEALLOCATE( rwork  )
  IF( ALLOCATED( iwork  ) .and. SIZE( iwork  ) < liwork ) DEALLOCATE( iwork  )
  IF( ALLOCATED( ifail  ) .and. SIZE( ifail  ) < n      ) DEALLOCATE( ifail  )
  IF( ALLOCATED( hdiag  ) .and. SIZE( hdiag  ) < n      ) DEALLOCATE( hdiag  )
  IF( ALLOCATED( sdiag  ) .and. SIZE( sdiag  ) < n      ) DEALLOCATE( sdiag  )
  IF( ALLOCATED( h_temp ) .and. SIZE( h_temp ) < n*n    ) DEALLOCATE( h_temp )
  IF( ALLOCATED( s_temp ) .and. SIZE( s_temp ) < n*n    ) DEALLOCATE( s_temp )

  IF( .not. ALLOCATED( work   ) ) ALLOCATE( work( max_lwork )   )
  IF( .not. ALLOCATED( rwork  ) ) ALLOCATE( rwork( max_lrwork )  )
  IF( .not. ALLOCATED( iwork  ) ) ALLOCATE( iwork( max_liwork )  )
  IF( .not. ALLOCATED( ifail  ) ) ALLOCATE( ifail( 2*ldh )  )
  IF( .not. ALLOCATED( hdiag  ) ) ALLOCATE( hdiag( 2*ldh )  )
  IF( .not. ALLOCATED( sdiag  ) ) ALLOCATE( sdiag( 2*ldh )  )
  IF( .not. ALLOCATED( h_temp ) ) ALLOCATE( h_temp(ldh,2*ldh) )
  IF( .not. ALLOCATED( s_temp ) ) ALLOCATE( s_temp(ldh,2*ldh) )

#ifdef USE_GPU
  IF( cpu_path==0 ) THEN
     IF( ALLOCATED( work_d   ) .and. SIZE( work_d   ) < lwork_d  ) DEALLOCATE( work_d   )
     IF( ALLOCATED( rwork_d  ) .and. SIZE( rwork_d  ) < lrwork_d ) DEALLOCATE( rwork_d  )
     IF( ALLOCATED( hdiag_d  ) .and. SIZE( hdiag_d  ) < n      ) DEALLOCATE( hdiag_d  )
     IF( ALLOCATED( sdiag_d  ) .and. SIZE( sdiag_d  ) < n      ) DEALLOCATE( sdiag_d  )

     IF( .not. ALLOCATED( work_d   ) ) ALLOCATE( work_d( max_lwork_d )   )
     IF( .not. ALLOCATED( rwork_d  ) ) ALLOCATE( rwork_d( max_lrwork_d ) )
     IF( .not. ALLOCATED( hdiag_d  ) ) ALLOCATE( hdiag_d( 2*ldh )        )
     IF( .not. ALLOCATED( sdiag_d  ) ) ALLOCATE( sdiag_d( 2*ldh )        )
  ENDIF
#endif

#ifdef USE_GPU
  if( cpu_path==0 ) then
   !call cufMemcpy2D( s_temp_d, size( s_temp_d, 1 ), s, ldh, ldh, n )
    call copy_diag(sdiag_d, s, ldh, n)
  else
#endif
    s_temp(1:ldh,1:n) = s(1:ldh,1:n)  
    !DO i = 1, n
    !  sdiag(i) = DBLE( s(i,i) )
    !END DO
#ifdef USE_GPU
  endif   
#endif

#ifndef USE_GPU
    !
    all_eigenvalues = ( m == n )
    !
    !
    IF ( all_eigenvalues ) THEN
       !
       IF( SIZE(rwork) < 3*n-2) THEN
         DEALLOCATE( rwork )
         ALLOCATE( rwork( 3*n - 2 ) )
       ENDIF
       !
       ! ... calculate all eigenvalues (overwritten to v)
       !
       v(:,:) = h(:,:)
       !
       CALL ZHEGV( 1, 'V', 'U', n, v, ldh, &
                   s, ldh, e, work, lwork, rwork, info )
       !
    ELSE
#endif
       !
       ! ... save the diagonal of input H (it will be overwritten)
       !
#ifdef USE_GPU
       if( cpu_path==0 ) then
         !call cufMemcpy2D( h_temp_d, size( h_temp_d, 1 ), h, ldh, ldh, n )
         call copy_diag(hdiag_d, h, ldh, n)
       else
#endif
         h_temp(1:ldh,1:n) = h(1:ldh,1:n)  
         !DO i = 1, n
         !  hdiag(i) = DBLE( h(i,i) )
         !END DO
#ifdef USE_GPU
       endif   
#endif

        !
        ! ... calculate only m lowest eigenvalues
        !
        abstol = 0.D0
       ! abstol = 2.D0*DLAMCH( 'S' )
        !
        ! ... the following commented lines calculate optimal lwork
        !
        !lwork = -1
        !
        !CALL ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
        !             0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
        !             work, lwork, rwork, iwork, ifail, info )
        !
        !lwork = INT( work(1) ) + 1
        !
        !IF( lwork > SIZE( work ) ) THEN
        !   DEALLOCATE( work )
        !   ALLOCATE( work( lwork ) )
        !END IF
        !
#ifdef USE_GPU
      if( cpu_path == 0 ) then

#ifdef USE_GPU_MPI
        call zhegvdx_gpu(n, h, ldh, s, ldh, v, ldh, 1, m, e, work_d,&
                         lwork_d, rwork_d, lrwork_d, &
                         work, lwork, rwork, lrwork, &
                         iwork, liwork, v_h, size(v_h, 1), e_h, info, .TRUE.)
#else
        call zhegvdx_gpu(n, h, ldh, s, ldh, v, ldh, 1, m, e, work_d,&
                         lwork_d, rwork_d, lrwork_d, &
                         work, lwork, rwork, lrwork, &
                         iwork, liwork, v_h, size(v_h, 1), e_h, info)
#endif

        ! Note: if zhegvdx_gpu fails, info = -1
        mm = m

     else

        CALL ZHEGVX( 1, 'V', 'I', 'U', n, h_temp, ldh, s_temp, ldh, &
                     0.D0, 0.D0, 1, m, abstol, mm, e_h, v_h, ldh, &
                     work, lwork, rwork, iwork, ifail, info )

     endif
#else
        CALL ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
                     0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
                     work, lwork, rwork, iwork, ifail, info )
#endif
        !
        !
        ! ... restore input H matrix from saved diagonal and lower triangle
        !
#ifdef USE_GPU
    if(cpu_path==0) then
      !call cufMemcpy2D( h, ldh, h_temp_d, size( h_temp_d, 1 ), ldh, n )
      call restore_upper_tri(h, ldh, hdiag_d, n)
    endif
#else
    h(1:ldh,1:n) = h_temp(1:ldh,1:n)
    !DO i = 1, n
    !  h_temp(i,i) = CMPLX( hdiag(i), 0.0_DP ,kind=DP)
    !  DO j = i + 1, n
    !    h_temp(i,j) = CONJG( h_temp(j,i) )
    !  END DO
    !  DO j = n + 1, ldh
    !    h_temp(j,i) = ( 0.0_DP, 0.0_DP )
    !  END DO
    !END DO
#endif

        !
#ifndef USE_GPU
     END IF
#endif
     !
     !
     IF ( info > n ) THEN
        CALL errore( 'cdiaghg', 'S matrix not positive definite', ABS( info ) )
     ELSE IF ( info > 0 ) THEN
        CALL errore( 'cdiaghg', 'eigenvectors failed to converge', ABS( info ) )
     ELSE IF ( info < 0 ) THEN
        CALL errore( 'cdiaghg', 'incorrect call to ZHEGV*', ABS( info ) )
     END IF
     !
     ! ... restore input S matrix from saved diagonal and lower triangle
     !
#ifdef USE_GPU
    if(cpu_path==0) then
     !call EpMemcpy2D( s, ldh, s_temp_d, size( s_temp_d, 1 ), ldh, n )
      call restore_upper_tri(s, ldh, sdiag_d, n)
    endif
#else
    s(1:ldh,1:n) = s_temp(1:ldh,1:n)
    !DO i = 1, n
    !  s(i,i) = CMPLX( sdiag(i), 0.0_DP ,kind=DP)
    !  DO j = i + 1, n
    !    s_temp(i,j) = CONJG( s_temp(j,i) )
    !  END DO
    !  DO j = n + 1, ldh
    !    s_temp(j,i) = ( 0.0_DP, 0.0_DP )
    !  END DO
    !END DO
#endif

    !
  END IF
  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
#ifdef USE_GPU
#ifdef USE_GPU_MPI
  istat = cudaDeviceSynchronize()
  CALL mp_bcast( e(1:n), root_bgrp, intra_bgrp_comm )
  CALL mp_bcast( v(1:ldh,1:m), root_bgrp, intra_bgrp_comm )
#else
  CALL mp_bcast( e_h(1:n), root_bgrp, intra_bgrp_comm )
  CALL mp_bcast( v_h(1:ldh,1:m), root_bgrp, intra_bgrp_comm )
  e(1:n) = e_h(1:n)
  istat = cudaMemcpy2D(v, ldh, v_h, size(v_h,1), n, m)
#endif
#else
  CALL mp_bcast( e, root_bgrp, intra_bgrp_comm )
  CALL mp_bcast( v, root_bgrp, intra_bgrp_comm )
#endif
  !
  IF( first_time_zhegvx ) first_time_zhegvx = 0
  !
  CALL stop_clock( 'cdiaghg' )
  !
  RETURN
  !
END SUBROUTINE MY_ROUTINE( cdiaghg )

#ifndef USE_GPU
!----------------------------------------------------------------------------
SUBROUTINE pcdiaghg( n, h, s, ldh, e, v, desc )
  !----------------------------------------------------------------------------
  !
  ! ... calculates eigenvalues and eigenvectors of the generalized problem
  ! ... Hv=eSv, with H hermitean matrix, S overlap matrix.
  ! ... On output both matrix are unchanged
  !
  ! ... Parallel version, with full data distribution
  !
  USE kinds,            ONLY : DP
  USE mp,               ONLY : mp_bcast
  USE mp_bands,         ONLY : root_bgrp, intra_bgrp_comm
  USE zhpev_module,     ONLY : pzhpev_drv, zhpev_drv
  USE descriptors,      ONLY : la_descriptor
  USE parallel_toolkit, ONLY : zsqmdst, zsqmcll
#if defined __SCALAPACK
  USE mp_diag,          ONLY : ortho_cntx, me_blacs, np_ortho, me_ortho, ortho_comm
  USE zhpev_module,     ONLY : pzheevd_drv
#endif

  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: n, ldh
    ! dimension of the matrix to be diagonalized
    ! leading dimension of h, as declared in the calling pgm unit
  COMPLEX(DP), INTENT(INOUT) :: h(ldh,ldh), s(ldh,ldh)
    ! actually intent(in) but compilers don't know and complain
    ! matrix to be diagonalized
    ! overlap matrix
  REAL(DP), INTENT(OUT) :: e(n)
    ! eigenvalues
  COMPLEX(DP), INTENT(OUT) :: v(ldh,ldh)
    ! eigenvectors (column-wise)
  TYPE(la_descriptor), INTENT(IN) :: desc
  !
  INTEGER             :: nx
#if defined __SCALAPACK
  INTEGER             :: descsca( 16 ), info
#endif
    ! local block size
  COMPLEX(DP), ALLOCATABLE :: ss(:,:), hh(:,:), tt(:,:)
    ! work space used only in parallel diagonalization
  !
  ! ... input s and h are copied so that they are not destroyed
  !
  CALL start_clock( 'cdiaghg' )
  !
  IF( desc%active_node > 0 ) THEN
     !
     nx   = desc%nrcx
     !
     IF( nx /= ldh ) &
        CALL errore(" pcdiaghg ", " inconsistent leading dimension ", ldh )
     !
     ALLOCATE( hh( nx, nx ) )
     ALLOCATE( ss( nx, nx ) )
     !
     hh(1:nx,1:nx) = h(1:nx,1:nx)
     ss(1:nx,1:nx) = s(1:nx,1:nx)
     !
  END IF

  CALL start_clock( 'cdiaghg:choldc' )
  !
  ! ... Cholesky decomposition of sl ( L is stored in sl )
  !
  IF( desc%active_node > 0 ) THEN
     !
#if defined __SCALAPACK
     CALL descinit( descsca, n, n, desc%nrcx, desc%nrcx, 0, 0, ortho_cntx, SIZE( ss, 1 ) , info )
     !
     IF( info /= 0 ) CALL errore( ' cdiaghg ', ' desckinit ', ABS( info ) )
#endif
     !
#if defined __SCALAPACK

     CALL pzpotrf( 'L', n, ss, 1, 1, descsca, info )

     IF( info /= 0 ) CALL errore( ' cdiaghg ', ' problems computing cholesky ', ABS( info ) )
#else
     CALL qe_pzpotrf( ss, nx, n, desc )
#endif
     !
  END IF
  !
  CALL stop_clock( 'cdiaghg:choldc' )
  !
  ! ... L is inverted ( sl = L^-1 )
  !
  CALL start_clock( 'cdiaghg:inversion' )
  !
  IF( desc%active_node > 0 ) THEN
     !
#if defined __SCALAPACK
     !CALL clear_upper_tr( ss )
     ! set to zero the upper triangle of ss
     !
     CALL sqr_zsetmat( 'U', n, ZERO, ss, size(ss,1), desc )
     !
     CALL pztrtri( 'L', 'N', n, ss, 1, 1, descsca, info )
     !
     IF( info /= 0 ) CALL errore( ' cdiaghg ', ' problems computing inverse ', ABS( info ) )
#else
     CALL qe_pztrtri( ss, nx, n, desc )
#endif
     !
  END IF
  !
  CALL stop_clock( 'cdiaghg:inversion' )
  !
  ! ... vl = L^-1*H
  !
  CALL start_clock( 'cdiaghg:paragemm' )
  !
  IF( desc%active_node > 0 ) THEN
     !
     CALL sqr_zmm_cannon( 'N', 'N', n, ONE, ss, nx, hh, nx, ZERO, v, nx, desc )
     !
  END IF
  !
  ! ... hl = ( L^-1*H )*(L^-1)^T
  !
  IF( desc%active_node > 0 ) THEN
     !
     CALL sqr_zmm_cannon( 'N', 'C', n, ONE, v, nx, ss, nx, ZERO, hh, nx, desc )
     !
     ! ensure that "hh" is really Hermitian, it is sufficient to set the diagonal
     ! properly, because only the lower triangle of hh will be used
     ! 
     CALL sqr_zsetmat( 'H', n, ZERO, hh, size(hh,1), desc )
     !
  END IF
  !
  CALL stop_clock( 'cdiaghg:paragemm' )
  !
  !
  IF ( desc%active_node > 0 ) THEN
     ! 
#ifdef TEST_DIAG
     CALL test_drv_begin()
#endif

#if defined(__SCALAPACK)
     !
     CALL pzheevd_drv( .true., n, desc%nrcx, hh, e, ortho_cntx, ortho_comm )
     !
#else
     !
     CALL qe_pzheevd( .true., n, desc, hh, SIZE( hh, 1 ), e )
     !
#endif
     !
#ifdef TEST_DIAG
     CALL test_drv_end()
#endif
     !
  END IF
  !
  ! ... v = (L^T)^-1 v
  !
  CALL start_clock( 'cdiaghg:paragemm' )
  !
  IF ( desc%active_node > 0 ) THEN
     !
     CALL sqr_zmm_cannon( 'C', 'N', n, ONE, ss, nx, hh, nx, ZERO, v, nx, desc )
     !
  END IF
  !
  CALL mp_bcast( e, root_bgrp, intra_bgrp_comm )
  !
  CALL stop_clock( 'cdiaghg:paragemm' )
  !
  IF ( desc%active_node > 0 ) THEN
     DEALLOCATE( ss, hh )
  END IF
  !
  CALL stop_clock( 'cdiaghg' )
  !
  RETURN
  !
CONTAINS
  !
  SUBROUTINE test_drv_begin()
     ALLOCATE( tt( n, n ) )
     CALL zsqmcll( n, hh, nx, tt, n, desc, desc%comm )
     RETURN
  END SUBROUTINE test_drv_begin
  !
  SUBROUTINE test_drv_end()
     !
     INTEGER :: i, j, k
     COMPLEX(DP), ALLOCATABLE :: diag(:,:)
     !
     IF( desc%myc == 0 .AND. desc%myr == 0 ) THEN

        write( 100, fmt="(A20,2D18.10)" ) ' e code = ', e( 1 ), e( n )
        ALLOCATE( diag( n*(n+1)/2, 1 ) )
        k = 1
        ! write( 100, fmt="(I5)" ) n
        DO j = 1, n
           DO i = j, n
              diag( k, 1 ) = tt( i, j )
              ! write( 100, fmt="(2I5,2D18.10)" ) i, j, tt( i, j )
              k = k + 1
           END DO
        END DO
        call zhpev_drv( 'V', 'L', N, diag(:,1), e, tt, n )
        write( 100, fmt="(A20,2D18.10)" ) ' e test = ', e( 1 ), e( n )
        ! write( 100, * ) 'eigenvalues and eigenvectors'
        DO j = 1, n
           ! write( 100, fmt="(1I5,1D18.10,A)" ) j, e( j )
           DO i = 1, n
              ! write( 100, fmt="(2I5,2D18.10)" ) i, j, tt( i, j )
           END DO
        END DO
        close(100)
        DEALLOCATE( diag )
     END IF
     CALL mp_bcast( tt, 0, desc%comm )
     CALL zsqmdst( n, tt, n, hh, nx, desc )
     DEALLOCATE( tt )
     CALL errore('cdiaghg','stop serial',1)
     RETURN
  END SUBROUTINE test_drv_end
  !
END SUBROUTINE pcdiaghg
!
#endif
