!
! Copyright (C) Quantum ESPRESSO group
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

! This mode is set by default
#define __FFTW_ALL_XY_PLANES

!=----------------------------------------------------------------------=!
   MODULE fft_scalar_essl
!=----------------------------------------------------------------------=!

       USE fft_param

       USE, intrinsic ::  iso_c_binding
#ifdef USE_CUDA
       USE cudafor
       USE cufft
#endif
       IMPLICIT NONE
       SAVE

#if defined(__LINUX_ESSL)

        PRIVATE
        PUBLIC :: cft_1z, cft_2xy, cfft3d, cfft3ds

       INTERFACE cft_1z
          MODULE PROCEDURE cft_1z_cpu
#ifdef USE_CUDA
          MODULE PROCEDURE cft_1z_gpu
#endif
       END INTERFACE

       INTERFACE cft_2xy
          MODULE PROCEDURE cft_2xy_cpu
#ifdef USE_CUDA
          MODULE PROCEDURE cft_2xy_gpu
#endif
       END INTERFACE

       INTERFACE cfft3d
          MODULE PROCEDURE cfft3d_cpu
#ifdef USE_CUDA
          MODULE PROCEDURE cfft3d_gpu
#endif
       END INTERFACE

! ...   Local Parameter


        !   Workspace that is statically allocated is defined here
        !   in order to avoid multiple copies of the same workspace
        !   lwork:   Dimension of the work space array (if any)

        !   ESSL IBM library: see the ESSL manual for DCFT

        INTEGER, PARAMETER :: lwork = 20000 + ( 2*nfftx + 256 ) * 64 + 3*nfftx
        REAL (DP) :: work( lwork )


!=----------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------=!

!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT along "z"
!
!
!
!=----------------------------------------------------------------------=!
!

   SUBROUTINE cft_1z_cpu(c, nsl, nz, ldz, isign, cout)

!     driver routine for nsl 1d complex fft's of length nz
!     ldz >= nz is the distance between sequences to be transformed
!     (ldz>nz is used on some architectures to reduce memory conflicts)
!     input  :  c(ldz*nsl)   (complex)
!     output : cout(ldz*nsl) (complex - NOTA BENE: transform is not in-place!)
!     isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
!     Up to "ndims" initializations (for different combinations of input
!     parameters nz, nsl, ldz) are stored and re-used if available

     INTEGER, INTENT(IN) :: isign
     INTEGER, INTENT(IN) :: nsl, nz, ldz

     COMPLEX (DP) :: c(:), cout(:)

     REAL (DP)  :: tscale
     INTEGER    :: i, err, idir, ip, void
     INTEGER, SAVE :: zdims( 3, ndims ) = -1
     INTEGER, SAVE :: icurrent = 1
     LOGICAL :: done

     INTEGER :: tid

     ! ...   Machine-Dependent parameters, work arrays and tables of factors

     !   ltabl   Dimension of the tables of factors calculated at the
     !           initialization stage

#if defined(__OPENMP)
     INTEGER :: offset, ldz_t
     INTEGER :: omp_get_max_threads
     EXTERNAL :: omp_get_max_threads
#endif


     !   ESSL IBM library: see the ESSL manual for DCFT

     INTEGER, PARAMETER :: ltabl = 20000 + 3 * nfftx
     REAL (DP), SAVE :: fw_tablez( ltabl, ndims )
     REAL (DP), SAVE :: bw_tablez( ltabl, ndims )


     IF( nsl < 0 ) THEN
       CALL fftx_error__(" fft_scalar: cft_1z ", " nsl out of range ", nsl)
     END IF

     !
     !   Here initialize table only if necessary
     !

     DO ip = 1, ndims

        !   first check if there is already a table initialized
        !   for this combination of parameters

        done = ( nz == zdims(1,ip) )

        !   The initialization in ESSL and FFTW v.3 depends on all three parameters

        done = done .AND. ( nsl == zdims(2,ip) ) .AND. ( ldz == zdims(3,ip) )
        IF (done) EXIT
     END DO

     IF( .NOT. done ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       ! WRITE( stdout, fmt="('DEBUG cft_1z, reinitializing tables ', I3)" ) icurrent

       tscale = 1.0_DP / nz

       CALL DCFT ( 1, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl,  1, &
          tscale, fw_tablez(1, icurrent), ltabl, work(1), lwork)
       CALL DCFT ( 1, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, -1, &
          1.0_DP, bw_tablez(1, icurrent), ltabl, work(1), lwork)

       zdims(1,icurrent) = nz; zdims(2,icurrent) = nsl; zdims(3,icurrent) = ldz;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined(__FFT_CLOCKS)
     CALL start_clock( 'cft_1z' )
#endif


     ! essl uses a different convention for forward/backward transforms
     ! wrt most other implementations: notice the sign of "idir"

     IF( isign < 0 ) THEN
        idir   =+1
        tscale = 1.0_DP / nz
        CALL DCFT (0, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, idir, &
             tscale, fw_tablez(1, ip), ltabl, work, lwork)
     ELSE IF( isign > 0 ) THEN
        idir   =-1
        tscale = 1.0_DP
        CALL DCFT (0, c(1), 1, ldz, cout(1), 1, ldz, nz, nsl, idir, &
             tscale, bw_tablez(1, ip), ltabl, work, lwork)
     END IF

#if defined(__FFT_CLOCKS)
     CALL stop_clock( 'cft_1z' )
#endif

     RETURN

   END SUBROUTINE cft_1z_cpu

#ifdef USE_CUDA
   SUBROUTINE cft_1z_gpu(c, nsl, nz, ldz, isign, cout)

!     driver routine for nsl 1d complex fft's of length nz
!     ldz >= nz is the distance between sequences to be transformed
!     (ldz>nz is used on some architectures to reduce memory conflicts)
!     input  :  c(ldz*nsl)   (complex)
!     ### GPU VERION IN PLACE!!! #### output : cout(ldz*nsl) (complex - NOTA
!     BENE: transform is not in-place!)
!     isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
!     Up to "ndims" initializations (for different combinations of input
!     parameters nz, nsl, ldz) are stored and re-used if available
#ifdef TRACK_FLOPS
     USE flops_tracker, ONLY : fft_ops
#endif
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: isign
     INTEGER, INTENT(IN) :: nsl, nz, ldz

     COMPLEX (DP), DEVICE :: c(:), cout(:)

     REAL (DP)  :: tscale
     INTEGER    :: i, err, idir, ip, void, istat
#ifdef TRACK_FLOPS
     REAL (DP), SAVE :: zflops( ndims ) = 0.d0
#endif
     INTEGER, SAVE :: zdims( 3, ndims ) = -1
     INTEGER, SAVE :: icurrent = 1
     LOGICAL :: found

     INTEGER :: tid

     !   Pointers to the "C" structures containing FFT factors ( PLAN )
     !   C_POINTER is defined in include/fft_defs.h
     !   for 32bit executables, C_POINTER is integer(4)
     !   for 64bit executables, C_POINTER is integer(8)

     INTEGER, SAVE :: cufft_planz( ndims ) = 0
     !C_POINTER, SAVE :: fw_planz( ndims ) = 0
     !C_POINTER, SAVE :: bw_planz( ndims ) = 0

     IF( nsl < 0 ) THEN
       CALL fftx_error__(" fft_scalar: cft_1z ", " nsl out of range ", nsl)
     END IF

     !
     !   Here initialize table only if necessary
     !
     CALL lookup()

     IF( .NOT. found ) THEN

       !   no table exist for these parameters
       !   initialize a new one
#if defined(__GPU_DEBUG)
       print *,"INIT CUFFT Z PLAN: ",nz,"x",nsl
#endif
       CALL init_plan()

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined(__FFT_CLOCKS)
     CALL start_clock( 'GPU_cft_1z' )
#endif

     IF (isign < 0) THEN
        !print *,"exec cufft FWD",nz,ldz,nsl
        !call flush(6)
        istat = cufftExecZ2Z( cufft_planz( ip), c(1), c(1), CUFFT_FORWARD )
        tscale = 1.0_DP / nz
!$cuf kernel do(1) <<<*,*>>>
        DO i=1, ldz * nsl
           cout( i ) = c( i ) * tscale
        END DO
     ELSE IF (isign > 0) THEN
        !print *,"exec cufft INV",nz,ldz,nsl
        !call flush(6)
        istat = cufftExecZ2Z( cufft_planz( ip), c(1), cout(1), CUFFT_INVERSE )
!CUFFT_FORWARD )
     END IF

#if defined(__FFT_CLOCKS)
     CALL stop_clock( 'GPU_cft_1z' )
#endif

#ifdef TRACK_FLOPS
     fft_ops = fft_ops + zflops( ip )
#endif

     RETURN

     CONTAINS 

     SUBROUTINE lookup()
     DO ip = 1, ndims
        !   first check if there is already a table initialized
        !   for this combination of parameters
        found = ( nz == zdims(1,ip) ) .AND. ( nsl == zdims(2,ip) ) .AND. ( ldz == zdims(3,ip) )
        IF (found) EXIT
     END DO
     END SUBROUTINE lookup

     SUBROUTINE init_plan()
       IMPLICIT NONE
       INTEGER, PARAMETER :: RANK=1
       INTEGER :: FFT_DIM(RANK), DATA_DIM(RANK)
       INTEGER :: STRIDE, DIST, BATCH

        FFT_DIM(1) = nz
       DATA_DIM(1) = ldz
            STRIDE = 1
              DIST = ldz
             BATCH = nsl

       IF( cufft_planz( icurrent) /= 0 ) istat = cufftDestroy( cufft_planz( icurrent) )

       istat = cufftPlanMany( cufft_planz( icurrent), RANK, FFT_DIM, & 
                              DATA_DIM, STRIDE, DIST, &
                              DATA_DIM, STRIDE, DIST, &
                              CUFFT_Z2Z, BATCH )

#ifdef TRACK_FLOPS
       zflops( icurrent ) = 5.0d0 * REAL( nz ) * log( REAL( nz ) )/log( 2.d0 )
#endif

       zdims(1,icurrent) = nz; zdims(2,icurrent) = nsl; zdims(3,icurrent) = ldz;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan

   END SUBROUTINE cft_1z_gpu

#endif

!
!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT along "x" and "y" direction
!
!
!
!=----------------------------------------------------------------------=!
!
!

   SUBROUTINE cft_2xy_cpu(r, nzl, nx, ny, ldx, ldy, isign, pl2ix)

!     driver routine for nzl 2d complex fft's of lengths nx and ny
!     input : r(ldx*ldy)  complex, transform is in-place
!     ldx >= nx, ldy >= ny are the physical dimensions of the equivalent
!     2d array: r2d(ldx, ldy) (x first dimension, y second dimension)
!     (ldx>nx, ldy>ny used on some architectures to reduce memory conflicts)
!     pl2ix(nx) (optional) is 1 for columns along y to be transformed
!     isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
!     Up to "ndims" initializations (for different combinations of input
!     parameters nx,ny,nzl,ldx) are stored and re-used if available

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: isign, ldx, ldy, nx, ny, nzl
     INTEGER, OPTIONAL, INTENT(IN) :: pl2ix(:)
     COMPLEX (DP) :: r( : )
     INTEGER :: i, k, j, err, idir, ip, kk, void
     REAL(DP) :: tscale
     INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: dims( 4, ndims) = -1
     LOGICAL :: dofft( nfftx ), done
     INTEGER, PARAMETER  :: stdout = 6

#if defined(__OPENMP)
     INTEGER :: offset
     INTEGER :: nx_t, ny_t, nzl_t, ldx_t, ldy_t
     INTEGER  :: itid, mytid, ntids
     INTEGER  :: omp_get_thread_num, omp_get_num_threads,omp_get_max_threads
     EXTERNAL :: omp_get_thread_num, omp_get_num_threads, omp_get_max_threads
#endif

     INTEGER, PARAMETER :: ltabl = 20000 + 3 * nfftx
     REAL (DP), SAVE :: fw_tablex( ltabl, ndims ), fw_tabley( ltabl, ndims )
     REAL (DP), SAVE :: bw_tablex( ltabl, ndims ), bw_tabley( ltabl, ndims )


     dofft( 1 : nx ) = .TRUE.
     IF( PRESENT( pl2ix ) ) THEN
       IF( SIZE( pl2ix ) < nx ) &
         CALL fftx_error__( ' cft_2xy ', ' wrong dimension for arg no. 8 ', 1 )
       DO i = 1, nx
         IF( pl2ix(i) < 1 ) dofft( i ) = .FALSE.
       END DO
     END IF

     ! WRITE( stdout,*) 'DEBUG: ', COUNT( dofft )

     !
     !   Here initialize table only if necessary
     !

     DO ip = 1, ndims

       !   first check if there is already a table initialized
       !   for this combination of parameters

       done = ( ny == dims(1,ip) ) .AND. ( nx == dims(3,ip) )
       done = done .AND. ( ldx == dims(2,ip) ) .AND.  ( nzl == dims(4,ip) )
       IF (done) EXIT

     END DO

     IF( .NOT. done ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       ! WRITE( stdout, fmt="('DEBUG cft_2xy, reinitializing tables ', I3)" ) icurrent

#if defined(__OPENMP)

       tscale = 1.0_DP / ( nx * ny )
       CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, nx,  1, 1.0_DP, &
          fw_tabley( 1, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, nx, -1, 1.0_DP, &
          bw_tabley(1, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny,  1, &
          tscale, fw_tablex( 1, icurrent), ltabl, work(1), lwork)
       CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny, -1, &
          1.0_DP, bw_tablex(1, icurrent), ltabl, work(1), lwork)

#else

       tscale = 1.0_DP / ( nx * ny )
       CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, 1,  1, 1.0_DP, &
          fw_tabley( 1, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, 1, -1, 1.0_DP, &
          bw_tabley(1, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny,  1, &
          tscale, fw_tablex( 1, icurrent), ltabl, work(1), lwork)
       CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny, -1, &
          1.0_DP, bw_tablex(1, icurrent), ltabl, work(1), lwork)

#endif

       dims(1,icurrent) = ny; dims(2,icurrent) = ldx;
       dims(3,icurrent) = nx; dims(4,icurrent) = nzl;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined(__FFT_CLOCKS)
     CALL start_clock( 'cft_2xy' )
#endif


#if defined(__OPENMP)

   IF( isign < 0 ) THEN
      tscale = 1.0_DP / ( nx * ny )
      do k = 1, nzl
         kk = 1 + ( k - 1 ) * ldx * ldy
         CALL DCFT ( 0, r( kk ), 1, ldx, r( kk ), 1, ldx, nx, ny, &
              1, tscale, fw_tablex( 1, ip ), ltabl, work( 1 ), lwork)
         CALL DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, nx, &
              1, 1.0_DP, fw_tabley(1, ip), ltabl, work( 1 ), lwork)
      end do
   ELSE IF( isign > 0 ) THEN
      DO k = 1, nzl
         kk = 1 + ( k - 1 ) * ldx * ldy
         CALL DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, nx, &
                   -1, 1.0_DP, bw_tabley(1, ip), ltabl, work( 1 ), lwork)
         CALL DCFT ( 0, r( kk ), 1, ldx, r( kk ), 1, ldx, nx, ny, &
                   -1, 1.0_DP, bw_tablex(1, ip), ltabl, work( 1 ), lwork)
      END DO
   END IF

#else

   IF( isign < 0 ) THEN
      idir = 1
      tscale = 1.0_DP / ( nx * ny )
      do k = 1, nzl
         kk = 1 + ( k - 1 ) * ldx * ldy
         CALL DCFT ( 0, r(kk), 1, ldx, r(kk), 1, ldx, nx, ny, idir, &
              tscale, fw_tablex( 1, ip ), ltabl, work( 1 ), lwork)
         do i = 1, nx
            IF( dofft( i ) ) THEN
               kk = i + ( k - 1 ) * ldx * ldy
               call DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, 1, &
                    idir, 1.0_DP, fw_tabley(1, ip), ltabl, work( 1 ), lwork)
            END IF
         end do
      end do
   ELSE IF( isign > 0 ) THEN
      idir = -1
      DO k = 1, nzl
         do i = 1, nx
            IF( dofft( i ) ) THEN
               kk = i + ( k - 1 ) * ldx * ldy
               call DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, 1, &
                    idir, 1.0_DP, bw_tabley(1, ip), ltabl, work( 1 ), lwork)
            END IF
         end do
         kk = 1 + ( k - 1 ) * ldx * ldy
         CALL DCFT ( 0, r( kk ), 1, ldx, r( kk ), 1, ldx, nx, ny, idir, &
              1.0_DP, bw_tablex(1, ip), ltabl, work( 1 ), lwork)
      END DO
   END IF
#endif

#if defined(__FFT_CLOCKS)
     CALL stop_clock( 'cft_2xy' )
#endif

     RETURN

   END SUBROUTINE cft_2xy_cpu

#ifdef USE_CUDA
   SUBROUTINE cft_2xy_gpu(r, temp, nzl, nx, ny, ldx, ldy, isign, pl2ix)

!     driver routine for nzl 2d complex fft's of lengths nx and ny
!     input : r(ldx*ldy)  complex, transform is in-place
!     ldx >= nx, ldy >= ny are the physical dimensions of the equivalent
!     2d array: r2d(ldx, ldy) (x first dimension, y second dimension)
!     (ldx>nx, ldy>ny used on some architectures to reduce memory conflicts)
!     pl2ix(nx) (optional) is 1 for columns along y to be transformed
!     isign > 0 : forward (f(G)=>f(R)), isign <0 backward (f(R) => f(G))
!     Up to "ndims" initializations (for different combinations of input
!     parameters nx,ny,nzl,ldx) are stored and re-used if available
#ifdef TRACK_FLOPS
     USE flops_tracker, ONLY : fft_ops
#endif
     IMPLICIT NONE

     INTEGER, INTENT(IN) :: isign, ldx, ldy, nx, ny, nzl
     INTEGER, OPTIONAL, INTENT(IN) :: pl2ix(:)
!pgi$ ignore_tkr r, temp
     COMPLEX (DP), DEVICE :: r(ldx,ldy,nzl), temp(ldy,nzl,ldx)
     INTEGER :: i, k, j, err, idir, ip, kk, void, istat
     REAL(DP) :: tscale
     INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: dims( 6, ndims) = -1
     ! dims(5,:) = batch_1
     ! dims(6,:) = batch_2
     LOGICAL :: dofft( nfftx ), found
     INTEGER, PARAMETER  :: stdout = 6
#ifdef TRACK_FLOPS
     REAL (DP), SAVE :: xyflops( ndims ) = 0.d0
#endif

#if defined(__FFTW_ALL_XY_PLANES)
     INTEGER, SAVE :: cufft_plan_2d( ndims ) = 0
#else
     INTEGER, SAVE :: cufft_plan_x( ndims ) = 0
     INTEGER, SAVE :: cufft_plan_y( 2, ndims ) = 0
#endif
     INTEGER :: batch_1, batch_2
     !C_POINTER, SAVE :: fw_plan_2d( ndims ) = 0
     !C_POINTER, SAVE :: bw_plan_2d( ndims ) = 0

     dofft( 1 : nx ) = .TRUE.
     batch_1 = nx
     batch_2 = 0
     IF( PRESENT( pl2ix ) ) THEN
       IF( SIZE( pl2ix ) < nx ) &
         CALL fftx_error__( ' cft_2xy ', ' wrong dimension for arg no. 8 ', 1 )
       DO i = 1, nx
         IF( pl2ix(i) < 1 ) dofft( i ) = .FALSE.
       END DO

       i=1
       do while(pl2ix(i) >= 1 .and. i<=nx); i=i+1; END DO
       batch_1 = i-1
       do while(pl2ix(i) < 1 .and. i<=nx); i=i+1; END DO
       batch_2 = nx-i+1
#if 0
       !batch_2_start = i
       !do while(pl2ix(i) >= 1 .and. i<nx); i=i+1; END DO
      
       do while( i<=nx ) 
         do while(pl2ix(i) < 1 .and. i<=nx); i=i+1; END DO
         batch_start = i
         do while(pl2ix(i) >= 1 .and. i<=nx); i=i+1; END DO
         batch_end = i-1
         batch_count = batch_end - batch_start + 1
!         print *,"batch: ",batch_start,batch_end,batch_count
       enddo
#endif
     END IF

     !
     !   Here initialize table only if necessary
     !

     CALL lookup()

     IF( .NOT. found ) THEN

       !   no table exist for these parameters
       !   initialize a new one
       CALL init_plan()

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined(__FFT_CLOCKS)
     CALL start_clock( 'GPU_cft_2xy' )
#endif
! TODO: consider case where not all XY planes are processed

     IF( isign < 0 ) THEN
        !
        !print *,"exec cufft FWD",nx,ny,ldx,ldy,nzl

        tscale = 1.0_DP / ( nx * ny )
        !
#if defined(__FFTW_ALL_XY_PLANES)
        istat = cufftExecZ2Z( cufft_plan_2d(ip), r(1,1,1), r(1,1,1),
CUFFT_FORWARD )
#else
        istat = cufftExecZ2Z( cufft_plan_x(ip), r(1,1,1), r(1,1,1), CUFFT_FORWARD )
        if(istat) print *,"error in fftxy fftx istat = ",istat

!$cuf kernel do(3) <<<*,(16,16,1)>>>
        DO k=1, nzl
           DO i=1, ldx 
              DO j=1, ldy
                temp(j,k,i) = r(i,j,k)                
              END DO
           END DO
        END DO
        
        
        if(batch_1>0) then
           istat = cufftExecZ2Z( cufft_plan_y(1,ip), temp(1,1,1), temp(1,1,1), CUFFT_FORWARD )
           if(istat) print *,"error in fftxy ffty batch_1 istat = ",istat
        end if

        if(batch_2>0) then     
           istat = cufftExecZ2Z( cufft_plan_y(2,ip), temp(1,1,nx-batch_2+1), temp(1,1,nx-batch_2+1), CUFFT_FORWARD )
           if(istat) print *,"error in fftxy ffty batch_2 istat = ",istat
        end if

!$cuf kernel do(3) <<<*,(16,16,1)>>>
        DO k=1, nzl
           DO j=1, ldy
             DO i=1, ldx
                r(i,j,k) = temp(j,k,i) * tscale
              END DO
           END DO
        END DO
#endif

#if 0
!$cuf kernel do(1) <<<*,*>>>
        DO i=1, ldx * ldy * nzl
           r( i ) = r( i ) * tscale
        END DO
#endif
        !CALL ZDSCAL( ldx * ldy * nzl, tscale, r(1), 1)
        !
     ELSE IF( isign > 0 ) THEN
        !
        !print *,"exec cufft INV",nx,ny,ldx,ldy,nzl
#if defined(__FFTW_ALL_XY_PLANES)
        istat = cufftExecZ2Z( cufft_plan_2d(ip), r(1,1,1), r(1,1,1), CUFFT_INVERSE )
#else
!$cuf kernel do(3) <<<*,(16,16,1)>>>
        DO k=1, nzl
           DO i=1, ldx
              DO j=1, ldy
                temp(j,k,i) = r(i,j,k)
              END DO
           END DO
        END DO

        if(batch_1>0) then
           istat = cufftExecZ2Z( cufft_plan_y(1,ip), temp(1,1,1), temp(1,1,1), CUFFT_INVERSE )
           if(istat) print *,"error in fftxy ffty batch_1 istat = ",istat
        end if

        if(batch_2>0) then
           istat = cufftExecZ2Z( cufft_plan_y(2,ip), temp(1,1,nx-batch_2+1), temp(1,1,nx-batch_2+1), CUFFT_INVERSE )
           if(istat) print *,"error in fftxy ffty batch_2 istat = ",istat
        end if

!$cuf kernel do(3) <<<*,(16,16,1)>>>
        DO k=1, nzl
           DO j=1, ldy
             DO i=1, ldx
                r(i,j,k) = temp(j,k,i)
              END DO
           END DO
        END DO

!        do i = 1, nx
!           IF( dofft( i ) ) THEN
!             istat = cufftExecZ2Z( cufft_plan_y(ip), r(i), r(i), CUFFT_INVERSE
!             )
!             if(istat) print *,"error in fftxy ffty istat = ",istat,i
!           END IF
!        end do

        istat = cufftExecZ2Z( cufft_plan_x(ip), r(1,1,1), r(1,1,1), CUFFT_INVERSE )
        if(istat) print *,"error in fftxy fftx istat = ",istat

#endif
        !
     END IF


#if defined(__FFT_CLOCKS)
     CALL stop_clock( 'GPU_cft_2xy' )
#endif

#ifdef TRACK_FLOPS
     fft_ops = fft_ops + xyflops( ip )
#endif

     RETURN

   CONTAINS

     SUBROUTINE lookup()
     DO ip = 1, ndims
       !   first check if there is already a table initialized
       !   for this combination of parameters
       found = ( ny == dims(1,ip) ) .AND. ( nx == dims(3,ip) )
       found = found .AND. ( ldx == dims(2,ip) ) .AND.  ( nzl == dims(4,ip) )
       found = found .AND. ( batch_1 == dims(5,ip) ) .AND. (batch_2 == dims(6,ip) )
       IF (found) EXIT
     END DO
     END SUBROUTINE lookup

     SUBROUTINE init_plan()
       IMPLICIT NONE
#if defined(__FFTW_ALL_XY_PLANES)
       INTEGER, PARAMETER :: RANK=2
       INTEGER :: FFT_DIM(RANK), DATA_DIM(RANK)
       INTEGER :: STRIDE, DIST, BATCH

        FFT_DIM(1) = ny
        FFT_DIM(2) = nx
       DATA_DIM(1) = ldy
       DATA_DIM(2) = ldx
            STRIDE = 1
              DIST = ldx*ldy
             BATCH = nzl

       IF( cufft_plan_2d( icurrent) /= 0 )  istat = cufftDestroy( cufft_plan_2d(icurrent) )

       istat = cufftPlanMany( cufft_plan_2d( icurrent), RANK, FFT_DIM, &
                              DATA_DIM, STRIDE, DIST, &
                              DATA_DIM, STRIDE, DIST, &
                              CUFFT_Z2Z, BATCH )
#else
       INTEGER, PARAMETER :: RANK=1
       INTEGER :: FFT_DIM_X(RANK), DATA_DIM_X(RANK), FFT_DIM_Y(RANK), DATA_DIM_Y(RANK)
       INTEGER :: STRIDE_X, STRIDE_Y, DIST_X, DIST_Y, BATCH_X, BATCH_Y1, BATCH_Y2

        FFT_DIM_X(1) = nx
       DATA_DIM_X(1) = ldx
            STRIDE_X = 1
              DIST_X = ldx
             BATCH_X = ny*nzl

        FFT_DIM_Y(1) = ny
       DATA_DIM_Y(1) = ldy
            STRIDE_Y = 1
              DIST_Y = ldy
            BATCH_Y1 = nzl*BATCH_1
            BATCH_Y2 = nzl*BATCH_2


       IF( cufft_plan_x( icurrent) /= 0 )  istat = cufftDestroy( cufft_plan_x(icurrent) ) 
       IF( cufft_plan_y( 1, icurrent) /= 0 )  istat = cufftDestroy( cufft_plan_y(1,icurrent) )
       IF( cufft_plan_y( 2, icurrent) /= 0 )  istat = cufftDestroy( cufft_plan_y(2,icurrent) )

#if defined(__GPU_DEBUG)
       print *,"INIT CUFFT XY PLAN: ",nx,"x",ny,"x",nzl
#endif

       istat = cufftPlanMany( cufft_plan_x( icurrent), RANK, FFT_DIM_X, &
                              DATA_DIM_X, STRIDE_X, DIST_X, &
                              DATA_DIM_X, STRIDE_X, DIST_X, &
                              CUFFT_Z2Z, BATCH_X )

       istat = cufftPlanMany( cufft_plan_y( 1, icurrent), RANK, FFT_DIM_Y, &
                              DATA_DIM_Y, STRIDE_Y, DIST_Y, &
                              DATA_DIM_Y, STRIDE_Y, DIST_Y, &
                              CUFFT_Z2Z, BATCH_Y1 )

       istat = cufftPlanMany( cufft_plan_y( 2, icurrent), RANK, FFT_DIM_Y, &
                              DATA_DIM_Y, STRIDE_Y, DIST_Y, &
                              DATA_DIM_Y, STRIDE_Y, DIST_Y, &
                              CUFFT_Z2Z, BATCH_Y2 )


#endif

#ifdef TRACK_FLOPS
       xyflops( icurrent ) = REAL( ny*nzl )                    * 5.0d0 * REAL( nx ) * log( REAL( nx )  )/log( 2.d0 ) &
                           + REAL( nzl*BATCH_1 + nzl*BATCH_2 ) * 5.0d0 * REAL( ny ) * log( REAL( ny )  )/log( 2.d0 )

#endif

       dims(1,icurrent) = ny; dims(2,icurrent) = ldx;
       dims(3,icurrent) = nx; dims(4,icurrent) = nzl;
       dims(5,icurrent) = BATCH_1; dims(6,icurrent) = BATCH_2;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan

   END SUBROUTINE cft_2xy_gpu
#endif


!
!=----------------------------------------------------------------------=!
!
!
!
!         3D scalar FFTs
!
!
!
!=----------------------------------------------------------------------=!
!

   SUBROUTINE cfft3d_cpu( f, nx, ny, nz, ldx, ldy, ldz, howmany, isign )

  !     driver routine for 3d complex fft of lengths nx, ny, nz
  !     input  :  f(ldx*ldy*ldz)  complex, transform is in-place
  !     ldx >= nx, ldy >= ny, ldz >= nz are the physical dimensions
  !     of the equivalent 3d array: f3d(ldx,ldy,ldz)
  !     (ldx>nx, ldy>ny, ldz>nz may be used on some architectures
  !      to reduce memory conflicts - not implemented for FFTW)
  !     isign > 0 : f(G) => f(R)   ; isign < 0 : f(R) => f(G)
  !
  !     Up to "ndims" initializations (for different combinations of input
  !     parameters nx,ny,nz) are stored and re-used if available

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nx, ny, nz, ldx, ldy, ldz, howmany, isign
     COMPLEX (DP) :: f(:)
     INTEGER :: i, k, j, err, idir, ip
     REAL(DP) :: tscale
     INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: dims(3,ndims) = -1

     IF ( nx < 1 ) &
         call fftx_error__('cfft3d',' nx is less than 1 ', 1)
     IF ( ny < 1 ) &
         call fftx_error__('cfft3d',' ny is less than 1 ', 1)
     IF ( nz < 1 ) &
         call fftx_error__('cfft3d',' nz is less than 1 ', 1)
     IF ( howmany /= 1 ) &
         call fftx_error__('cfft3d',' howmany different from 1, not yetimplemented for ESSL ', 1)
     !
     !   Here initialize table only if necessary
     !
     ip = -1
     DO i = 1, ndims

       !   first check if there is already a table initialized
       !   for this combination of parameters

       IF ( ( nx == dims(1,i) ) .and. &
            ( ny == dims(2,i) ) .and. &
            ( nz == dims(3,i) ) ) THEN
         ip = i
         EXIT
       END IF
     END DO

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       ! no initialization for 3d FFT's from ESSL

       dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1

     END IF

     !
     !   Now perform the 3D FFT using the machine specific driver
     !

     IF ( isign < 0 ) THEN
       tscale = 1.0_DP / ( nx * ny * nz )
       idir = +1
     ELSE IF( isign > 0 ) THEN
       tscale = 1.0_DP
       idir = -1
     END IF

     IF( isign /= 0 ) CALL dcft3( f(1), ldx,ldx*ldy, f(1), ldx,ldx*ldy, &
          nx,ny,nz, idir, tscale, work(1), lwork)

     RETURN
   END SUBROUTINE cfft3d_cpu

#ifdef USE_CUDA
SUBROUTINE cfft3d_gpu( f, nx, ny, nz, ldx, ldy, ldz, howmany, isign )

  !     driver routine for 3d complex fft of lengths nx, ny, nz
  !     input  :  f(ldx*ldy*ldz)  complex, transform is in-place
  !     ldx >= nx, ldy >= ny, ldz >= nz are the physical dimensions
  !     of the equivalent 3d array: f3d(ldx,ldy,ldz)
  !     (ldx>nx, ldy>ny, ldz>nz may be used on some architectures
  !      to reduce memory conflicts - not implemented for FFTW)
  !     isign > 0 : f(G) => f(R)   ; isign < 0 : f(R) => f(G)
  !
  !     Up to "ndims" initializations (for different combinations of input
  !     parameters nx,ny,nz) are stored and re-used if available

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nx, ny, nz, ldx, ldy, ldz, howmany, isign
     COMPLEX (DP), device :: f(:)
     INTEGER :: i, k, j, err, idir, ip, istat
     REAL(DP) :: tscale
     INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: dims(3,ndims) = -1

!     C_POINTER, save :: fw_plan(ndims) = 0
!     C_POINTER, save :: bw_plan(ndims) = 0
     INTEGER, SAVE :: cufft_plan_3d( ndims ) = 0


     IF ( nx < 1 ) &
         call fftx_error__('cfft3d',' nx is less than 1 ', 1)
     IF ( ny < 1 ) &
         call fftx_error__('cfft3d',' ny is less than 1 ', 1)
     IF ( nz < 1 ) &
         call fftx_error__('cfft3',' nz is less than 1 ', 1)

     !
     !   Here initialize table only if necessary
     !
     CALL lookup()

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       CALL init_plan()

     END IF

     !
     !   Now perform the 3D FFT using the machine specific driver
     !

     IF( isign < 0 ) THEN

        istat = cufftExecZ2Z( cufft_plan_3d(ip), f(1), f(1), CUFFT_FORWARD )

       tscale = 1.0_DP / DBLE( nx * ny * nz )
!$cuf kernel do(1) <<<*,*>>>
        DO i=1, nx*ny*nz
           f( i ) = f( i ) * tscale
        END DO
!       call ZDSCAL( nx * ny * nz, tscale, f(1), 1)

     ELSE IF( isign > 0 ) THEN

!       call FFTW_INPLACE_DRV_3D( bw_plan(ip), 1, f(1), 1, 1 )
        istat = cufftExecZ2Z( cufft_plan_3d(ip), f(1), f(1), CUFFT_INVERSE )

     END IF

     RETURN

   CONTAINS 

     SUBROUTINE lookup()
     ip = -1
     DO i = 1, ndims
       !   first check if there is already a table initialized
       !   for this combination of parameters
       IF ( ( nx == dims(1,i) ) .and. &
            ( ny == dims(2,i) ) .and. &
            ( nz == dims(3,i) ) ) THEN
         ip = i
         EXIT
       END IF
     END DO
     END SUBROUTINE lookup

     SUBROUTINE init_plan()
       INTEGER, PARAMETER :: RANK=3
       INTEGER :: FFT_DIM(RANK), DATA_DIM(RANK)
       INTEGER :: STRIDE, DIST, BATCH

        FFT_DIM(1) = nz
        FFT_DIM(2) = ny
        FFT_DIM(3) = nx
       DATA_DIM(1) = ldz
       DATA_DIM(2) = ldy
       DATA_DIM(3) = ldx
            STRIDE = 1
              DIST = ldx*ldy*ldz
             BATCH = 1

       IF( cufft_plan_3d( icurrent) /= 0 )  istat = cufftDestroy( cufft_plan_3d(icurrent) )

       istat = cufftPlanMany( cufft_plan_3d( icurrent), RANK, FFT_DIM, &
                              DATA_DIM, STRIDE, DIST, &
                              DATA_DIM, STRIDE, DIST, &
                              CUFFT_Z2Z, BATCH )

       !IF ( nx /= ldx .or. ny /= ldy .or. nz /= ldz ) &
       !  call fftx_error__('cfft3','not implemented',1)
       !IF( fw_plan(icurrent) /= 0 ) CALL DESTROY_PLAN_3D( fw_plan(icurrent) )
       !IF( bw_plan(icurrent) /= 0 ) CALL DESTROY_PLAN_3D( bw_plan(icurrent) )
       !idir = -1; CALL CREATE_PLAN_3D( fw_plan(icurrent), nx, ny, nz, idir)
       !idir =  1; CALL CREATE_PLAN_3D( bw_plan(icurrent), nx, ny, nz, idir)
       dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1
     END SUBROUTINE init_plan

   END SUBROUTINE cfft3d_gpu
#endif

!
!=----------------------------------------------------------------------=!
!
!
!
!         3D scalar FFTs,  but using sticks!
!
!
!
!=----------------------------------------------------------------------=!
!

SUBROUTINE cfft3ds (f, nx, ny, nz, ldx, ldy, ldz, howmany, isign, &
     do_fft_z, do_fft_y)
  !
  !     driver routine for 3d complex "reduced" fft - see cfft3d
  !     The 3D fft are computed only on lines and planes which have
  !     non zero elements. These lines and planes are defined by
  !     the two integer vectors do_fft_y(nx) and do_fft_z(ldx*ldy)
  !     (1 = perform fft, 0 = do not perform fft)
  !     This routine is implemented only for fftw, essl, acml
  !     If not implemented, cfft3d is called instead
  !
  !----------------------------------------------------------------------
  !
  implicit none

  integer :: nx, ny, nz, ldx, ldy, ldz, howmany, isign
  !
  !   logical dimensions of the fft
  !   physical dimensions of the f array
  !   sign of the transformation

  complex(DP) :: f ( ldx * ldy * ldz )
  integer :: do_fft_y(:), do_fft_z(:)
  !
  integer :: m, incx1, incx2
  INTEGER :: i, k, j, err, idir, ip,  ii, jj
  REAL(DP) :: tscale
  INTEGER, SAVE :: icurrent = 1
  INTEGER, SAVE :: dims(3,ndims) = -1

  INTEGER, PARAMETER :: ltabl = 20000 + 3 * nfftx
  REAL (DP), SAVE :: fw_table( ltabl, 3, ndims )
  REAL (DP), SAVE :: bw_table( ltabl, 3, ndims )

  tscale = 1.0_DP

  ! WRITE( stdout, fmt="('DEBUG cfft3ds :',6I6)") nx, ny, nz, ldx, ldy, ldz
  ! WRITE( stdout, fmt="('DEBUG cfft3ds :',24I2)") do_fft_y
  ! WRITE( stdout, fmt="('DEBUG cfft3ds :',24I2)") do_fft_z


  IF( ny /= ldy ) &
    CALL fftx_error__(' cfft3ds ', ' wrong dimensions: ny /= ldy ', 1 )

     IF ( howmany /= 1 ) &
         call fftx_error__('cfft3ds',' howmany different from 1, not yetimplemented for ESSL ', 1)

     ip = -1
     DO i = 1, ndims

       !   first check if there is already a table initialized
       !   for this combination of parameters

       IF( ( nx == dims(1,i) ) .and. ( ny == dims(2,i) ) .and. &
           ( nz == dims(3,i) ) ) THEN
         ip = i
         EXIT
       END IF

     END DO

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       !
       ! ESSL sign convention for fft's is the opposite of the "usual" one
       !
       tscale = 1.0_DP
       !  x - direction
       incx1 = 1; incx2 = ldx; m = ldy*nz
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nx, m,  1, 1.0_DP, &
          fw_table( 1, 1, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nx, m, -1, 1.0_DP, &
          bw_table(1, 1, icurrent), ltabl, work(1), lwork )
       !  y - direction
       incx1 = ldx; incx2 = ldx*ldy; m = nz;
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, ny, m,  1, 1.0_DP, &
          fw_table( 1, 2, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, ny, m, -1, 1.0_DP, &
          bw_table(1, 2, icurrent), ltabl, work(1), lwork )
       !  z - direction
       incx1 = ldx * ldy; incx2 = 1; m = 1
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nz, m,  1, 1.0_DP, &
          fw_table(1, 3, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nz, m, -1, 1.0_DP, &
          bw_table(1, 3, icurrent), ltabl, work(1), lwork )


       dims(1,icurrent) = nx; dims(2,icurrent) = ny; dims(3,icurrent) = nz
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1

     END IF


     IF ( isign > 0 ) THEN

        !
        !  k-direction ...
        !

        incx1 = ldx * ldy;  incx2 = 1;  m = 1
        do i =1,nx
           do j =1,ny
              ii = i + ldx *(j -1)
              if ( do_fft_z(ii) == 1 ) then
                 call dcft (0, f( ii ), incx1, incx2, f( ii ), incx1, incx2, nz, m, &
                   -isign, 1.0_DP, bw_table ( 1, 3, ip ), ltabl, work( 1 ), lwork)
              end if
           end do
        end do

        !
        !  ... j-direction ...
        !

        incx1 = ldx;  incx2 = ldx*ldy;  m = nz

        do i = 1, nx
           if ( do_fft_y( i ) == 1 ) then
              call dcft (0, f (i), incx1, incx2, f (i), incx1, incx2, nx, m, &
                -isign, 1.0_DP, bw_table ( 1, 2,  ip ), ltabl, work( 1 ), lwork)
           endif
        enddo

        !
        !  ... i - direction
        !

        incx1 = 1;  incx2 = ldx;  m = ldy * nz

        call dcft (0, f (1), incx1,incx2, f (1), incx1,incx2, nx, m, &
          -isign, 1.0_DP, bw_table ( 1, 1,  ip ), ltabl, work( 1 ), lwork)


     ELSE

        !
        !  i - direction ...
        !

        incx1 = 1;  incx2 = ldx;  m = ldy*nz

        call dcft (0, f (1), incx1,incx2, f (1), incx1,incx2, nx, m, &
          -isign, 1.0_DP, fw_table ( 1, 1, ip ), ltabl, work( 1 ), lwork)

        !
        !  ... j-direction ...
        !

        incx1 = ldx;  incx2 = ldx*ldy;  m = nz

        do i = 1, nx
           if ( do_fft_y ( i ) == 1 ) then
              call dcft (0, f (i), incx1, incx2, f (i), incx1, incx2, ny, m, &
                -isign, 1.0_DP, fw_table ( 1, 2, ip ), ltabl, work( 1 ), lwork)
           endif
        enddo

        !
        !  ... k-direction
        !

        incx1 = ldx * ny;  incx2 = 1;  m = 1

        do i = 1, nx
           do j = 1, ny
              ii = i + ldx * (j-1)
              if (do_fft_z(ii) == 1 ) then
                 call dcft (0, f( ii ), incx1, incx2, f( ii ), incx1, incx2, nz, m, &
                  -isign, 1.0_DP, fw_table ( 1, 3, ip ), ltabl, work( 1 ), lwork)
              end if
           end do
        end do

        call DSCAL (2 * ldx * ldy * nz, 1.0_DP/(nx * ny * nz), f(1), 1)

     END IF
     RETURN
   END SUBROUTINE cfft3ds
#endif
!=----------------------------------------------------------------------=!
 END MODULE fft_scalar_essl
!=----------------------------------------------------------------------=!

