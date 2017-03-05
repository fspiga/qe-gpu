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
  INTEGER, PARAMETER :: jdr_min_size=1000
  COMPLEX(DP), PINNED, ALLOCATABLE :: Z(:,:)
  REAL(DP), ALLOCATABLE, PINNED :: e_h(:)
  COMPLEX(DP), ALLOCATABLE, PINNED :: v_h(:,:)
  REAL(DP),    DEVICE, ALLOCATABLE :: rwork_d(:), sdiag_d(:), hdiag_d(:)
  COMPLEX(DP), DEVICE, ALLOCATABLE :: work_d(:)
  COMPLEX(DP), DEVICE, ALLOCATABLE :: h_temp_d(:,:),s_temp_d(:,:)
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
  USE cublas
  USE zheevd_jdr,       ONLY : zheevd_gpu
  USE zhegvx_module,    ONLY : h_temp=>h_temp_d, s_temp=>s_temp_d, hdiag=>hdiag_d, sdiag=>sdiag_d
  USE zhegvx_module,    ONLY : e_h, v_h, rwork_d, work_d, Z, jdr_min_size
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
  INTEGER  :: cuf_i, ii, istat
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
     nb = 64
     lwork  = max(2*n + n*n, n + n*nb)
     lrwork = 1 + 5*n + 2*n*n
     liwork = 3 + 5*n
     max_lwork  = 2*(max(2*ldh + ldh*ldh, ldh + ldh*nb))
     max_lrwork = 2*(1 + 5*ldh + 2*ldh*ldh)
     max_liwork = 2*(3 + 5*ldh)
#else
     nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
     lwork = max(2*n, (nb+1)*n)
     lrwork = 7*n
     liwork = 5*n
     nb = max(nb, ILAENV( 1, 'ZHETRD', 'U', ldh, -1, -1, -1 ))
     max_lwork = 2*max(2*ldh, (nb+1)*ldh)
     max_lrwork = 2*7*ldh
     max_liwork = 2*5*ldh
#endif

  IF(first_time_zhegvx) THEN
     ALLOCATE( work( max_lwork ) )
     ALLOCATE( rwork( max_lrwork ) )
     ALLOCATE( iwork( max_liwork ) )
     ALLOCATE( ifail( ldh ) )
     ALLOCATE( sdiag( ldh ) )
     ALLOCATE( hdiag( ldh ) )
     ALLOCATE( h_temp(ldh,ldh) )
     ALLOCATE( s_temp(ldh,ldh) )
#ifdef USE_GPU
     ALLOCATE( work_d( max_lwork ) )
     ALLOCATE( rwork_d( max_lrwork ) )
     ALLOCATE( Z(ldh,ldh) )
#endif
  ELSE IF( size( work ) < 2*lwork .or. size( rwork ) < 2*lrwork .or. size( iwork ) < 2*liwork ) then
     DEALLOCATE( work, rwork, iwork, ifail, sdiag, hdiag, h_temp, s_temp )
#ifdef USE_GPU
     DEALLOCATE( work_d, rwork_d, Z )
#endif
     ALLOCATE( work( max_lwork ) )
     ALLOCATE( rwork( max_lrwork ) )
     ALLOCATE( iwork( max_liwork ) )
     ALLOCATE( ifail( ldh ) )
     ALLOCATE( sdiag( ldh ) )
     ALLOCATE( hdiag( ldh ) )
     ALLOCATE( h_temp(ldh,ldh) )
     ALLOCATE( s_temp(ldh,ldh) )
#ifdef USE_GPU
     ALLOCATE( work_d( max_lwork ) )
     ALLOCATE( rwork_d( max_lrwork ) )     
     ALLOCATE( Z(ldh,ldh) )
#endif
  END IF

#if 1
    s_temp = s
#else
#ifdef USE_GPU
!$cuf kernel do(1) <<<*,*>>>
#endif
     DO i = 1, n
        sdiag(i) = DBLE( s(i,i) )
     END DO
#endif

#ifndef USE_GPU
     !
     all_eigenvalues = ( m == n )
     !
     !
     IF ( all_eigenvalues ) THEN
        !
        ALLOCATE( rwork( 3*n - 2 ) )
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
#if 1
    h_temp = h
#else
#ifdef USE_GPU
!$cuf kernel do(1) <<<*,*>>>
#endif
        DO i = 1, n
           hdiag(i) = DBLE( h(i,i) ) !  <--- THIS LINE GIVES ERROR: INVALID DEVICE FUNCTION
        END DO
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
        call magmaf_zpotrf_gpu('U', n, s, ldh, info)
        call magmaf_zhegst_gpu( 1, 'U', n, h, ldh, s, ldh, info)
        IF(n > jdr_min_size) THEN
           call zheevd_gpu('V', 'U', 1, m, n, h, ldh, v, ldh, e, work_d, 2*lwork, rwork_d, 2*(1+5*n+2*n*n),  &
                                                                 work  , 2*lwork, rwork  , 2*(1+5*n+2*n*n), iwork  , 2*(3+5*n), &
                                                    v_h, ldh, e_h )
           mm = m
        ELSE
           call magmaf_zheevdx_gpu('V', 'I', 'U', n, h, ldh, 0.D0, 0.D0, 1, m, mm, e_h, Z, ldh, work, 2*lwork, rwork, 2*(1+5*n+2*n*n), iwork, 2*(3+5*n), info)
           v_h = h(:,1:m)
           v = v_h
           e = e_h
        END IF
        call cublasZtrsm( 'L', 'U', 'N', 'N', n, mm, ONE, s, ldh, v, ldh )
!        h = h_temp
!        s = s_temp
        v_h = v
        e_h = e
#else
        CALL ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
                     0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
                     work, lwork, rwork, iwork, ifail, info )
#endif
        !
        !
        ! ... restore input H matrix from saved diagonal and lower triangle
        !
#if 1
        h = h_temp
#else
#ifdef USE_GPU
!$cuf kernel do(1) <<<*,*>>>
#endif
        DO i = 1, n
           h(i,i) = CMPLX( hdiag(i), 0.0_DP ,kind=DP)
           DO j = i + 1, n
              h(i,j) = CONJG( h(j,i) )
           END DO
           DO j = n + 1, ldh
              h(j,i) = ( 0.0_DP, 0.0_DP )
           END DO
        END DO
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
#if 1
     s = s_temp
#else
#ifdef USE_GPU
!$cuf kernel do(1) <<<*,*>>>
#endif
     DO i = 1, n
        s(i,i) = CMPLX( sdiag(i), 0.0_DP ,kind=DP)
        DO j = i + 1, n
           s(i,j) = CONJG( s(j,i) )
        END DO
        DO j = n + 1, ldh
           s(j,i) = ( 0.0_DP, 0.0_DP )
        END DO
     END DO
#endif
     !
  END IF
  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
  CALL mp_bcast( e, root_bgrp, intra_bgrp_comm )
  CALL mp_bcast( v, root_bgrp, intra_bgrp_comm )
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
