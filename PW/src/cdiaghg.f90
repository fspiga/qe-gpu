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
!----------------------------------------------------------------------------
SUBROUTINE cdiaghg( n, m, h, s, ldh, e, v )
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
  !
  INTEGER                  :: lwork, nb, mm, info, i, j
    ! mm = number of calculated eigenvectors
  REAL(DP)                 :: abstol
  INTEGER,     ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP),    ALLOCATABLE :: rwork(:), sdiag(:), hdiag(:)
  COMPLEX(DP), ALLOCATABLE :: work(:)
    ! various work space
  LOGICAL                  :: all_eigenvalues
 ! REAL(DP), EXTERNAL       :: DLAMCH
  INTEGER,  EXTERNAL       :: ILAENV
    ! ILAENV returns optimal block size "nb"
  !
  !
  CALL start_clock( 'cdiaghg' )
  !
  ! ... only the first processor diagonalizes the matrix
  !
  IF ( me_bgrp == root_bgrp ) THEN
     !
     ! ... save the diagonal of input S (it will be overwritten)
     !
     ALLOCATE( sdiag( n ) )
     DO i = 1, n
        sdiag(i) = DBLE( s(i,i) )
     END DO
     !
     all_eigenvalues = ( m == n )
     !
     ! ... check for optimal block size
     !
     nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
     !
     IF ( nb < 1 .OR. nb >= n) THEN
        !
        lwork = 2*n
        !
     ELSE
        !
        lwork = ( nb + 1 )*n
        !
     END IF
     !
     ALLOCATE( work( lwork ) )
     !
     IF ( all_eigenvalues ) THEN
        !
        ALLOCATE( rwork( 3*n - 2 ) )
        !
        ! ... calculate all eigenvalues (overwritten to v)
        !
        v(:,:) = h(:,:)
        !
   print *,"CALLING CPU ZHEGV !!!!!!!!!"
        CALL ZHEGV( 1, 'V', 'U', n, v, ldh, &
                    s, ldh, e, work, lwork, rwork, info )
        !
     ELSE
        !
        ALLOCATE( rwork( 7*n ) )
        !
        ! ... save the diagonal of input H (it will be overwritten)
        !
        ALLOCATE( hdiag( n ) )
        DO i = 1, n
           hdiag(i) = DBLE( h(i,i) )
        END DO
        !
        ALLOCATE( iwork( 5*n ) )
        ALLOCATE( ifail( n ) )
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
   print *,"CALLING CPU ZHEGVX !!!!!!!!!"

        CALL ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
                     0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
                     work, lwork, rwork, iwork, ifail, info )
        !
        DEALLOCATE( ifail )
        DEALLOCATE( iwork )
        !
        ! ... restore input H matrix from saved diagonal and lower triangle
        !
        DO i = 1, n
           h(i,i) = CMPLX( hdiag(i), 0.0_DP ,kind=DP)
           DO j = i + 1, n
              h(i,j) = CONJG( h(j,i) )
           END DO
           DO j = n + 1, ldh
              h(j,i) = ( 0.0_DP, 0.0_DP )
           END DO
        END DO
        !
        DEALLOCATE( hdiag )
        !
     END IF
     !
     !
     DEALLOCATE( rwork )
     DEALLOCATE( work )
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
     DO i = 1, n
        s(i,i) = CMPLX( sdiag(i), 0.0_DP ,kind=DP)
        DO j = i + 1, n
           s(i,j) = CONJG( s(j,i) )
        END DO
        DO j = n + 1, ldh
           s(j,i) = ( 0.0_DP, 0.0_DP )
        END DO
     END DO
     !
     DEALLOCATE( sdiag )
     !
  END IF
  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
  CALL mp_bcast( e, root_bgrp, intra_bgrp_comm )
  CALL mp_bcast( v, root_bgrp, intra_bgrp_comm )
  !
  CALL stop_clock( 'cdiaghg' )
  !
  RETURN
  !
END SUBROUTINE cdiaghg
!

#ifdef USE_CUDA


module cdiaghg_compute_gpu_module
  use cudafor
  use kinds, ONLY : DP
IMPLICIT NONE
  INTEGER,     pinned, save, ALLOCATABLE :: iwork(:), ifail(:)
  REAL(DP),    pinned, save, ALLOCATABLE :: rwork(:), sdiag(:), hdiag(:)
  COMPLEX(DP), pinned, save, ALLOCATABLE :: work(:)
  INTEGER,     DEVICE, ALLOCATABLE :: iwork_d(:), ifail_d(:)
  REAL(DP),    DEVICE, ALLOCATABLE :: rwork_d(:)!, sdiag_d(:), hdiag_d(:)
  COMPLEX(DP), DEVICE, ALLOCATABLE :: work_d(:)
  COMPLEX(DP), DEVICE, ALLOCATABLE :: h_d_temp(:,:),s_d_temp(:,:)
  COMPLEX(DP), pinned, save, ALLOCATABLE :: A(:,:), Z(:,:)

contains
!----------------------------------------------------------------------------
SUBROUTINE cdiaghg_gpu(n, m, h, h_d, s, s_d, ldh, e, e_d, v, v_d )
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
  USE cublas
  USE zheevd_jdr,       ONLY : zheevd_gpu
  USE nvtx
!  USE cudafor
  !USE magma
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
  !


  COMPLEX(DP), DEVICE, INTENT(INOUT) :: h_d(ldh,n), s_d(ldh,n)
    ! actually intent(in) but compilers don't know and complain
    ! matrix to be diagonalized
    ! overlap matrix
  REAL(DP), DEVICE, INTENT(OUT) :: e_d(n)
    ! eigenvalues
  COMPLEX(DP), DEVICE, INTENT(OUT) :: v_d(ldh,m)
    ! eigenvectors (column-wise)


  INTEGER                  :: lwork, lrwork, liwork,  nb, mm, info, i, j, cuf_i, ii, istat
    ! mm = number of calculated eigenvectors
  REAL(DP)                 :: abstol
  COMPLEX(DP) :: temp_v
    ! various work space
  LOGICAL                  :: all_eigenvalues
 ! REAL(DP), EXTERNAL       :: DLAMCH
  INTEGER,  EXTERNAL       :: ILAENV
    ! ILAENV returns optimal block size "nb"
  !
  REAL(DP) :: timer
#if 0
  INTEGER, SAVE :: counter=0
#endif
  INTEGER, PARAMETER :: jdr_min_size=1000
  INTEGER, save :: first_time_cdiaghg_gpu=1
  !
  CALL start_clock( 'cdiaghg' )
  !
  ! ... only the first processor diagonalizes the matrix
  !
  IF ( me_bgrp == root_bgrp ) THEN !(1) THEN !( me_bgrp == root_bgrp ) THEN
     !
     ! ... save the diagonal of input S (it will be overwritten)
     !
#if 0
     ALLOCATE( sdiag( n ) )
     DO i = 1, n
        sdiag(i) = DBLE( s(i,i) )
     END DO
#endif

#if 0     
     if(first_time_cdiaghg_gpu) ALLOCATE( sdiag_d( n ) )
!$cuf kernel do(1) <<<*,*>>>
     DO i = 1, n
        sdiag_d(i) = DBLE( s_d(i,i) )
     END DO
     !
     all_eigenvalues = ( m == n )
     !
     ! ... check for optimal block size
     !
!     nb = ILAENV( 1, 'ZHETRD', 'U', n, -1, -1, -1 )
     !
     IF ( nb < 64 ) nb = 64
     IF ( nb < 1 .OR. nb >= n) THEN
        !
        lwork = 2*n
        !
     ELSE
        !
        lwork = ( nb + 1 )*n
        !
     END IF
#endif
     nb = 64
     lwork  = max(2*n + n*n, n + n*nb)
     lrwork = 1 + 5*n + 2*n*n
     liwork = 3 + 5*n
     !
     if(first_time_cdiaghg_gpu) then 
        ALLOCATE( work  ( 2*( max(2*ldh + ldh*ldh,ldh + ldh*nb) ) ) )
        ALLOCATE( work_d( 2*( 2*ldh + ldh*ldh ) ) )
        ALLOCATE( rwork  ( 2*( 1 + 5*ldh + 2*ldh*ldh) ) ) !7*n ) )
        ALLOCATE( rwork_d( 2*( 1 + 5*ldh + 2*ldh*ldh) ) )
        ALLOCATE( iwork  ( 2*( 3 + 5*ldh ) ) )
!        ALLOCATE( hdiag( ldh ) )
!        ALLOCATE( hdiag_d( ldh ) )
        ALLOCATE( ifail( ldh ) )
        ALLOCATE( h_d_temp(ldh,ldh) )
        ALLOCATE( s_d_temp(ldh,ldh) )
!        ALLOCATE( A(ldh,ldh) )
        ALLOCATE( Z(ldh,ldh) )
     else
        if( size( work ) < 2*lwork .or. size( rwork ) < 2*lrwork .or. size( iwork ) < 2*liwork ) then
        print *,"buffers too small, re-allocating cdiaghg_gpu buffers..."
        DEALLOCATE( work, work_d, rwork, rwork_d, iwork, ifail, h_d_temp, s_d_temp, Z )
        ALLOCATE( work  ( 2*( max( 2*ldh + ldh*ldh, ldh + ldh*nb ) ) ) )
        ALLOCATE( work_d( 2*( 2*ldh + ldh*ldh ) ) )
        ALLOCATE( rwork  ( 2*( 1 + 5*ldh + 2*ldh*ldh) ) ) !7*n ) )
        ALLOCATE( rwork_d( 2*( 1 + 5*ldh + 2*ldh*ldh) ) )
        ALLOCATE( iwork  ( 2*( 3 + 5*ldh ) ) )
!        ALLOCATE( hdiag( ldh ) )
!        ALLOCATE( hdiag_d( ldh ) )
        ALLOCATE( ifail( ldh ) )
        ALLOCATE( h_d_temp(ldh,ldh) )
        ALLOCATE( s_d_temp(ldh,ldh) )
!        ALLOCATE( A(ldh,ldh) )
        ALLOCATE( Z(ldh,ldh) )
        endif
     endif
     !
     !TODO: consider doing ZHEGV for GPU version
     !      For now we force the ZHEGVX path
     IF (0) THEN !( all_eigenvalues ) THEN
        !
        !
        ! ... calculate all eigenvalues (overwritten to v)
        !
        v(:,:) = h_d(:,:)
        s = s_d
        !
        CALL ZHEGV( 1, 'V', 'U', n, v, ldh, &
                    s, ldh, e, work, lwork, rwork, info )
        e_d = e
        v_d = v
        !
     ELSE
        !
        !
        ! ... save the diagonal of input H (it will be overwritten)
        !
#if 0
        !ALLOCATE( hdiag( n ) )
        DO i = 1, n
           hdiag(i) = DBLE( h(i,i) )
        END DO
        !
     !ALLOCATE( hdiag_d( n ) )
!$cuf kernel do(1) <<<*,*>>>
     DO i = 1, n
        hdiag_d(i) = DBLE( h_d(i,i) )
     END DO

        !
        !ALLOCATE( iwork_d( 3+5*n ) )
        !ALLOCATE( ifail( n ) )
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
!        call compare(h, h_d, "h cdghg")
!        call compare(s, s_d, "s cdghg")
!        call compare(sdiag, sdiag_d, "sdiag")
!        call compare(hdiag, hdiag_d, "hdiag")
!        h_d = h
!        s_d = s

!        CALL ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
!                     0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
!                     work, lwork, rwork, iwork, ifail, info )

!         CALL magmaf_zhegvdx( 

!        v_d = v
!        e_d = e

!        ALLOCATE( h_d_temp(ldh,n) )
!        ALLOCATE( A(ldh,n) )
!        ALLOCATE( Z(ldh,n) )
#endif
        h_d_temp = h_d
        s_d_temp = s_d
#if 0
        counter = counter + 1 
        print *,"counter = ",counter

      print *,"n,m,ldh",n,m,ldh

    if(counter>=17) then
      !print *,"n,m,ldh",n,m,ldh
      h = h_d
      open(unit=10, status='replace', file="H_input", form='unformatted')
      write(10) h
      close(10)
      s = s_d
      open(unit=10, status='replace', file="S_input", form='unformatted')
      write(10) s
      close(10)
    endif
#endif
!        h = h_d
!        s = s_d

!        CALL magmaf_ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
!                     0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
!                     work, lwork, rwork, iwork, ifail, info )

!   e_d = e

!        call compare(e, e_d, "e zhegvx")
!        call compare(v, v_d, "v zhegvx")
#if 1

!#define DUMP

#ifdef DUMP
    s = s_d
    open(unit=10, status='replace', file="B_orig", form='unformatted')
    write(10) s
    close(10)
#endif

        !s = s_d
        !call zpotrf('U',n,s,ldh,info)
  CALL start_clock( 'zpotrf' )
        call magmaf_zpotrf_gpu('U', n, s_d, ldh, info)
!  istat = cudaDeviceSynchronize()
  CALL stop_clock(  'zpotrf' )

  CALL start_clock( 'zhegst' )

!        call compare(s,s_d,"s zpotrf")
        !h = h_d
        !call zhegst(1, 'U', n, h, ldh, s, ldh, info)
#ifdef DUMP
    print *,"n,ldh",n,ldh
    h = h_d
    open(unit=10, status='replace', file="A_input", form='unformatted')
    write(10) h
    close(10)
    s = s_d
    open(unit=10, status='replace', file="B_input", form='unformatted')
    write(10) s
    close(10)
#endif
        call magmaf_zhegst_gpu( 1, 'U', n, h_d, ldh, s_d, ldh, info)
#ifdef DUMP
    h = h_d
    open(unit=10, status='replace', file="A_output", form='unformatted')
    write(10) h
    close(10)
    STOP
#endif
        !call compare(h,h_d,"h zhegst")
        !call compare(s,s_d,"s zhegst")
!  istat = cudaDeviceSynchronize()
  CALL stop_clock( 'zhegst' )
!        v = ZERO
!        v_d = ZERO
!#define VERBOSE
#ifdef VERBOSE
 print *,"before zheevx..."
 print *,n,m,ldh
#endif
        !call ZHEEVX('V', 'I', 'U', n, h, ldh, 0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, work, lwork, rwork, iwork, ifail, info)
! e_d = e
! e = 0.d0
! print *,"after zheevx..."
! print *,"info = ", info
! do ii = 1, 10
!  print *,"v ",ii,v(ii,1),e(ii)
! end do
!  e = 0.0
!  e_d = 0.0
!  v_d = h_d_temp
  CALL start_clock( 'zheevx' )
if(n > jdr_min_size) then
#ifdef VERBOSE
timer = MPI_Wtime()
#endif
        call zheevd_gpu('V', 'U', 1, m, n, h_d, ldh, v_d, ldh, e_d, work_d, 2*lwork, rwork_d, 2*(1+5*n+2*n*n),  &
                                                                    work  , 2*lwork, rwork  , 2*(1+5*n+2*n*n), iwork  , 2*(3+5*n), &
                                                     v  , ldh, e )
mm = m
!!        call magmaf_zheevx_gpu('V', 'I', 'U', n, h_d, ldh, 0.D0, 0.D0, 1, m, abstol, mm, e, v_d, ldh, A, ldh, Z, ldh, work, lwork, rwork, iwork, ifail, info)
#ifdef VERBOSE
istat = cudaDeviceSynchronize()
timer = MPI_Wtime() - timer
print *,"*** zheevdx_jdr in ",timer, "secs",info
#endif
!  istat = cudaDeviceSynchronize()
  CALL stop_clock( 'zheevx' )

else
!  e_d = e
!  v = v_d
!  h_d = h_d_temp
  !v_d = h_d_temp
!  e = 0.0
!  h = h_d
#if 1
#ifdef VERBOSE
timer = MPI_Wtime()
#endif
!  call zheev('V', 'U', n, h, ldh, e, A, n*ldh, rwork, info)
        call magmaf_zheevdx_gpu('V', 'I', 'U', n, h_d, ldh, 0.D0, 0.D0, 1, m, mm, e, Z, ldh, work, 2*lwork, rwork, 2*(1+5*n+2*n*n), iwork, 2*(3+5*n), info)
#ifdef VERBOSE
timer = MPI_Wtime() - timer
  
print *,"*** zheevdx_magma in ",timer, "secs",info
#endif

#endif

  istat = cudaDeviceSynchronize()
  CALL stop_clock( 'zheevx' )

!   v = h_d(:,1:m)
!   v_d = v
!   e_d = e
#if 1
#if 0
  istat = cudaMemcpy2D( v_d(1,1), ldh, h_d(1,1), ldh, n, m , cudaMemcpyDeviceToDevice)
  if(istat) then
    print *,"ERROR in cudaMemcpy2D in cdiaghg.c"
    stop
  endif 
#else
  v = h_d(:,1:m)
  v_d = v
#endif
#endif
  !v = h_d(:,1:m)
  !h_d = h_d_temp
   e_d = e
endif

  !call compare(e(1:m), e_d(1:m), "e vdx")
  !call compare(v, v_d, "v vdx")
  !v_d = v 
  !e_d = e

! print *,"info_gpu = ",info
 !v = v_d
! do ii = 1, 10
!  temp_v = v_d(ii,1)
!  print *,"v ",ii,temp_v,e(ii)
! end do
!        call compare(e, e_d, "e b4 ztrsm")
!        call compare(v, v_d, "v b4 ztrsm")
  CALL start_clock( 'ztrsm' )
        call cublasZtrsm( 'L', 'U', 'N', 'N', n, mm, ONE, s_d, ldh, v_d, ldh )
!  istat = cudaDeviceSynchronize()
  CALL stop_clock( 'ztrsm' )

        !call ZTRSM( 'L', 'U', 'N', 'N', n, mm, ONE, s, ldh, v, ldh )

        h_d = h_d_temp
        s_d = s_d_temp
        v = v_d
        e = e_d
#if 0
      print *,"mm= ",mm

    if(counter>=17) then
      !print *,"mm= ",mm
      open(unit=10, status='replace', file="V_output", form='unformatted')
      write(10) v
      close(10)
      open(unit=10, status='replace', file="e_output", form='unformatted')
      write(10) e
      close(10)

      open(unit=10, status='old', file="H_input", form='unformatted')
      read(10) h
      close(10)

      open(unit=10, status='old', file="S_input", form='unformatted')
      read(10) s
      close(10)

        CALL ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
                     0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
                     work, lwork, rwork, iwork, ifail, info )

       call compare(e, e_d, "e comp")
       call compare(v, v_d, "v comp")

!      CALL magmaf_ZHEGVX( 1, 'V', 'I', 'U', n, h, ldh, s, ldh, &
!                          0.D0, 0.D0, 1, m, abstol, mm, e, v, ldh, &
!                          work, lwork, rwork, iwork, ifail, info )


    endif
#endif
!        call compare( e, e_d, "e eigenvalues")
!        call compare( v, v_d, "v eigenvectors")
        !
#endif
        !v_d = v
        !e = e_d
        !h_d = h
        !s_d = s

!        DEALLOCATE( work_d, rwork_d )
!        DEALLOCATE( A, Z, h_d_temp )
!        DEALLOCATE( ifail )
!        DEALLOCATE( iwork )

#if 0
        !h = h_d
        !
        ! ... restore input H matrix from saved diagonal and lower triangle
        !
        DO i = 1, n
           h(i,i) = CMPLX( hdiag(i), 0.0_DP ,kind=DP)
           DO j = i + 1, n
              h(i,j) = CONJG( h(j,i) )
           END DO
           DO j = n + 1, ldh
              h(j,i) = ( 0.0_DP, 0.0_DP )
           END DO
        END DO
        !
        DEALLOCATE( hdiag )
#endif
#if 0
!$cuf kernel do(1) <<<*,*>>>
     DO i = 1, n
           h_d(i,i) = CMPLX( hdiag_d(i), 0.0_DP ,kind=DP)
           DO j = i + 1, n
              h_d(i,j) = CONJG( h_d(j,i) )
           END DO
           DO j = n + 1, ldh
              h_d(j,i) = ( 0.0_DP, 0.0_DP )
           END DO
     END DO
#endif
!     call compare( h, h_d, "restore h")
!     DEALLOCATE( hdiag_d )
        !
     END IF
     !
     !
 !    DEALLOCATE( rwork )
 !    DEALLOCATE( work )
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
     !s = s_d
#if 0
     DO i = 1, n
        s(i,i) = CMPLX( sdiag(i), 0.0_DP ,kind=DP)
        DO j = i + 1, n
           s(i,j) = CONJG( s(j,i) )
        END DO
        DO j = n + 1, ldh
           s(j,i) = ( 0.0_DP, 0.0_DP )
        END DO
     END DO

     !
     DEALLOCATE( sdiag )
#endif
     !
#if 0
!$cuf kernel do(1) <<<*,*>>>
     DO i = 1, n
           s_d(i,i) = CMPLX( sdiag_d(i), 0.0_DP ,kind=DP)
           DO j = i + 1, n
              s_d(i,j) = CONJG( s_d(j,i) )
           END DO
           DO j = n + 1, ldh
              s_d(j,i) = ( 0.0_DP, 0.0_DP )
           END DO
     END DO
!     call compare( s, s_d, "restore s")
  !   DEALLOCATE( sdiag_d )
     !
#endif
  END IF

  !
  ! ... broadcast eigenvectors and eigenvalues to all other processors
  !
call nvtxStartRange("bcast",2)
  CALL mp_bcast( e, root_bgrp, intra_bgrp_comm )
  CALL mp_bcast( v, root_bgrp, intra_bgrp_comm )
call nvtxEndRange
  !
  e_d = e
  v_d = v
  CALL stop_clock( 'cdiaghg' )
  !
  first_time_cdiaghg_gpu = 0
  RETURN
  !
END SUBROUTINE cdiaghg_gpu
!
END MODULE cdiaghg_compute_gpu_module

#endif


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
