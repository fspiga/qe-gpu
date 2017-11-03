!
! Copyright (C) 2001-2008 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE gradcorr( rho, rhog, rho_core, rhog_core, etxc, vtxc, v )
  !----------------------------------------------------------------------------
  !
  USE constants,            ONLY : e2
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : nl, ngm, g
  USE lsda_mod,             ONLY : nspin
  USE cell_base,            ONLY : omega, alat
  USE funct,                ONLY : gcxc, gcx_spin, gcc_spin, igcc_is_lyp, &
                                   gcc_spin_more, dft_is_gradient, get_igcc
  USE spin_orb,             ONLY : domag
  USE noncollin_module,     ONLY : ux
  USE wavefunctions_module, ONLY : psic
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft

  !
  IMPLICIT NONE
  !
  REAL(DP),    INTENT(IN)    :: rho(dfftp%nnr,nspin), rho_core(dfftp%nnr)
  COMPLEX(DP), INTENT(IN)    :: rhog(ngm,nspin), rhog_core(ngm)
  REAL(DP),    INTENT(INOUT) :: v(dfftp%nnr,nspin)
  REAL(DP),    INTENT(INOUT) :: vtxc, etxc
  !
  INTEGER :: k, ipol, is, nspin0, ir, jpol
  !
  REAL(DP),    ALLOCATABLE :: grho(:,:,:), h(:,:,:), dh(:)
  REAL(DP),    ALLOCATABLE :: rhoout(:,:), segni(:), vgg(:,:), vsave(:,:)
  REAL(DP),    ALLOCATABLE :: gmag(:,:,:)

  COMPLEX(DP), ALLOCATABLE :: rhogsum(:,:)
  !
  REAL(DP) :: grho2(2), sx, sc, v1x, v2x, v1c, v2c, &
              v1xup, v1xdw, v2xup, v2xdw, v1cup, v1cdw , &
              etxcgc, vtxcgc, segno, arho, fac, zeta, rh, grh2, amag 
  REAL(DP) :: v2cup, v2cdw,  v2cud, rup, rdw, &
              grhoup, grhodw, grhoud, grup, grdw, seg
  !
  REAL(DP), PARAMETER :: epsr = 1.D-6, epsg = 1.D-10
  !
  !
  IF ( .NOT. dft_is_gradient() ) RETURN
  !
  call start_clock( 'gradcorr' )
  !
  etxcgc = 0.D0
  vtxcgc = 0.D0
  !
  nspin0=nspin
  if (nspin==4) nspin0=1
  if (nspin==4.and.domag) nspin0=2
  fac = 1.D0 / DBLE( nspin0 )
  !
  ALLOCATE(    h( 3, dfftp%nnr, nspin0) )
  ALLOCATE( grho( 3, dfftp%nnr, nspin0) )
  ALLOCATE( rhoout( dfftp%nnr, nspin0) )
  IF (nspin==4.AND.domag) THEN
     ALLOCATE( vgg( dfftp%nnr, nspin0 ) )
     ALLOCATE( vsave( dfftp%nnr, nspin ) )
     ALLOCATE( segni( dfftp%nnr ) )
     vsave=v
     v=0.d0
  ENDIF
  !
  ALLOCATE( rhogsum( ngm, nspin0 ) )
  !
  ! ... calculate the gradient of rho + rho_core in real space
  !
  IF ( nspin == 4 .AND. domag ) THEN
     !
     CALL compute_rho(rho,rhoout,segni,dfftp%nnr)
     !
     ! ... bring starting rhoout to G-space
     !
     DO is = 1, nspin0
        !
        psic(:) = rhoout(:,is)
        !
        CALL fwfft ('Dense', psic, dfftp)
        !
        rhogsum(:,is) = psic(nl(:))
        !
     END DO
  ELSE
     !
     rhoout(:,1:nspin0)  = rho(:,1:nspin0)
     rhogsum(:,1:nspin0) = rhog(:,1:nspin0)
     !
  ENDIF

  DO is = 1, nspin0
     !
     rhoout(:,is)  = fac * rho_core(:)  + rhoout(:,is)
     rhogsum(:,is) = fac * rhog_core(:) + rhogsum(:,is)
     !
call start_clock( 'gradrho' )
     CALL gradrho( dfftp%nnr, rhogsum(1,is), ngm, g, nl, grho(1,1,is) )
call stop_clock( 'gradrho' )
     !
  END DO

  !
  DEALLOCATE( rhogsum )
  !
  IF ( nspin0 == 1 ) THEN
     !
     ! ... This is the spin-unpolarised case
     !
     DO k = 1, dfftp%nnr
        !
        arho = ABS( rhoout(k,1) )
        !
        IF ( arho > epsr ) THEN
           !
           grho2(1) = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
           !
           IF ( grho2(1) > epsg ) THEN
              !
              segno = SIGN( 1.D0, rhoout(k,1) )
              !
              CALL gcxc( arho, grho2(1), sx, sc, v1x, v2x, v1c, v2c )
              !
              ! ... first term of the gradient correction : D(rho*Exc)/D(rho)
              !
              v(k,1) = v(k,1) + e2 * ( v1x + v1c )
              !
              ! ... h contains :
              !
              ! ...    D(rho*Exc) / D(|grad rho|) * (grad rho) / |grad rho|
              !
              h(:,k,1) = e2 * ( v2x + v2c ) * grho(:,k,1)
              !
              vtxcgc = vtxcgc+e2*( v1x + v1c ) * ( rhoout(k,1) - rho_core(k) )
              etxcgc = etxcgc+e2*( sx + sc ) * segno
              !
           ELSE
              h(:,k,1)=0.D0
           END IF
           !
        ELSE
           !
           h(:,k,1) = 0.D0
           !
        END IF
        !
     END DO
     !
  ELSE
     !
     ! ... spin-polarised case
     !
!$omp parallel do private( rh, grho2, sx, v1xup, v1xdw, v2xup, v2xdw, rup, rdw, &
!$omp             grhoup, grhodw, grhoud, sc, v1cup, v1cdw, v2cup, v2cdw, v2cud, &
!$omp             zeta, grh2, v2c, grup, grdw  ), &
!$omp             reduction(+:etxcgc,vtxcgc)
     DO k = 1, dfftp%nnr
        !
        rh = rhoout(k,1) + rhoout(k,2)
        !
        grho2(:) = grho(1,k,:)**2 + grho(2,k,:)**2 + grho(3,k,:)**2
        !
        CALL gcx_spin( rhoout(k,1), rhoout(k,2), grho2(1), &
                       grho2(2), sx, v1xup, v1xdw, v2xup, v2xdw )
        !
        IF ( rh > epsr ) THEN
           !
           IF ( igcc_is_lyp() ) THEN
              !
              rup = rhoout(k,1)
              rdw = rhoout(k,2)
              !
              grhoup = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
              grhodw = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
              !
              grhoud = grho(1,k,1) * grho(1,k,2) + &
                       grho(2,k,1) * grho(2,k,2) + &
                       grho(3,k,1) * grho(3,k,2)
              !
              CALL gcc_spin_more( rup, rdw, grhoup, grhodw, grhoud, &
                                  sc, v1cup, v1cdw, v2cup, v2cdw, v2cud )
              !
           ELSE
              !
              zeta = ( rhoout(k,1) - rhoout(k,2) ) / rh
              if (nspin.eq.4.and.domag) zeta=abs(zeta)*segni(k)
              !
              grh2 = ( grho(1,k,1) + grho(1,k,2) )**2 + &
                     ( grho(2,k,1) + grho(2,k,2) )**2 + &
                     ( grho(3,k,1) + grho(3,k,2) )**2
              !
              CALL gcc_spin( rh, zeta, grh2, sc, v1cup, v1cdw, v2c )
              !
              v2cup = v2c
              v2cdw = v2c
              v2cud = v2c
              !
           END IF
           !
        ELSE
           !
           sc    = 0.D0
           v1cup = 0.D0
           v1cdw = 0.D0
           v2c   = 0.D0
           v2cup = 0.D0
           v2cdw = 0.D0
           v2cud = 0.D0
           !
        ENDIF
        !
        ! ... first term of the gradient correction : D(rho*Exc)/D(rho)
        !
        v(k,1) = v(k,1) + e2 * ( v1xup + v1cup )
        v(k,2) = v(k,2) + e2 * ( v1xdw + v1cdw )
        !
        ! ... h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
        !
        DO ipol = 1, 3
           !
           grup = grho(ipol,k,1)
           grdw = grho(ipol,k,2)
           h(ipol,k,1) = e2 * ( ( v2xup + v2cup ) * grup + v2cud * grdw )
           h(ipol,k,2) = e2 * ( ( v2xdw + v2cdw ) * grdw + v2cud * grup )
           !
        END DO
        !
        vtxcgc = vtxcgc + &
                 e2 * ( v1xup + v1cup ) * ( rhoout(k,1) - rho_core(k) * fac )
        vtxcgc = vtxcgc + &
                 e2 * ( v1xdw + v1cdw ) * ( rhoout(k,2) - rho_core(k) * fac )
        etxcgc = etxcgc + e2 * ( sx + sc )
        !
     END DO
!$omp end parallel do
     !
  END IF
  !
  DO is = 1, nspin0
     !
     rhoout(:,is) = rhoout(:,is) - fac * rho_core(:)
     !
  END DO
  !
  DEALLOCATE( grho )
  !
  ALLOCATE( dh( dfftp%nnr ) )    
  !
  ! ... second term of the gradient correction :
  ! ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
  !
  DO is = 1, nspin0
     !
call start_clock( 'graddot' )
     CALL grad_dot( dfftp%nnr, h(1,1,is), ngm, g, nl, alat, dh )
call stop_clock( 'graddot' )
     !
     v(:,is) = v(:,is) - dh(:)
     !
     vtxcgc = vtxcgc - SUM( dh(:) * rhoout(:,is) )
     !
  END DO
  !
  vtxc = vtxc + omega * vtxcgc / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  etxc = etxc + omega * etxcgc / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )

  IF (nspin==4.AND.domag) THEN
     DO is=1,nspin0
        vgg(:,is)=v(:,is)
     ENDDO
     v=vsave
     DO k=1,dfftp%nnr
        v(k,1)=v(k,1)+0.5d0*(vgg(k,1)+vgg(k,2))
        amag=sqrt(rho(k,2)**2+rho(k,3)**2+rho(k,4)**2)
        IF (amag.GT.1.d-12) THEN
           v(k,2)=v(k,2)+segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,2)/amag
           v(k,3)=v(k,3)+segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,3)/amag
           v(k,4)=v(k,4)+segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,4)/amag
        ENDIF
     ENDDO
  ENDIF
  !
  DEALLOCATE( dh )
  DEALLOCATE( h )
  DEALLOCATE( rhoout )
  IF (nspin==4.and.domag) THEN
     DEALLOCATE( vgg )
     DEALLOCATE( vsave )
     DEALLOCATE( segni )
  ENDIF
  !

  call stop_clock( 'gradcorr' )
  RETURN
  !
END SUBROUTINE gradcorr
!
!----------------------------------------------------------------------------
SUBROUTINE gradrho( nrxx, a, ngm, g, nl, ga )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a in R-space (a is in G-space)
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE gvect,     ONLY : nlm
  USE control_flags, ONLY : gamma_only
  USE fft_base,      ONLY : dfftp
  USE fft_interfaces,ONLY : invfft

  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN)  :: nrxx
  INTEGER,     INTENT(IN)  :: ngm, nl(ngm)
  COMPLEX(DP), INTENT(IN)  :: a(ngm)
  REAL(DP),    INTENT(IN)  :: g(3,ngm)
  REAL(DP),    INTENT(OUT) :: ga(3,nrxx)
  !
  INTEGER                  :: ipol
  COMPLEX(DP), ALLOCATABLE :: gaux(:)
  !
  !
  ALLOCATE( gaux( nrxx ) )
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  ga(:,:) = 0.D0
  !
  DO ipol = 1, 3
     !
     gaux(:) = CMPLX(0.d0,0.d0,kind=dp)
     !
     gaux(nl(:)) = g(ipol,:) * CMPLX( -AIMAG( a(:) ), REAL( a(:) ) ,kind=DP)
     !
     IF ( gamma_only ) THEN
        !
        gaux(nlm(:)) = CMPLX( REAL( gaux(nl(:)) ), -AIMAG( gaux(nl(:)) ) ,kind=DP)
        !
     END IF
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL invfft ('Dense', gaux, dfftp)
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     ga(ipol,:) = ga(ipol,:) + tpiba * REAL( gaux(:) )
     !
  END DO
  !
  DEALLOCATE( gaux )
  !
  RETURN
  !
END SUBROUTINE gradrho
!
!----------------------------------------------------------------------------
SUBROUTINE gradient( nrxx, a, ngm, g, nl, ga )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a in R-space (a is also in R-space)
  !
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : nlm
  USE control_flags, ONLY : gamma_only
  USE fft_base,      ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nrxx
  INTEGER,  INTENT(IN)  :: ngm, nl(ngm)
  REAL(DP), INTENT(IN)  :: a(nrxx), g(3,ngm)
  REAL(DP), INTENT(OUT) :: ga(3,nrxx)
  !
  INTEGER                  :: ipol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)
  !
  !
  ALLOCATE(  aux( nrxx ) )
  ALLOCATE( gaux( nrxx ) )
  !
  aux = CMPLX( a(:), 0.D0 ,kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Dense', aux, dfftp)
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  DO ipol = 1, 3
     !
     gaux(:) = CMPLX(0.d0,0.d0, kind=dp)
     !
     gaux(nl(:)) = g(ipol,:) * &
                   CMPLX( -AIMAG( aux(nl(:)) ), REAL( aux(nl(:)) ) ,kind=DP)
     !
     IF ( gamma_only ) THEN
        !
        gaux(nlm(:)) = CMPLX( REAL( gaux(nl(:)) ), -AIMAG( gaux(nl(:)) ) ,kind=DP)
        !
     END IF
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL invfft ('Dense', gaux, dfftp)
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     ga(ipol,:) = tpiba * DBLE( gaux(:) )
     !
  END DO
  !
  DEALLOCATE( gaux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE gradient
!
!----------------------------------------------------------------------------
SUBROUTINE grad_dot( nrxx, a, ngm, g, nl, alat, da )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates da = \sum_i \grad_i a_i in R-space
  !
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : nlm
  USE control_flags, ONLY : gamma_only
  USE fft_base,      ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)     :: nrxx, ngm, nl(ngm)
  REAL(DP), INTENT(IN)     :: a(3,nrxx), g(3,ngm), alat
  REAL(DP), INTENT(OUT)    :: da(nrxx)
  !
  INTEGER                  :: n, ipol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)
  !
  !
  ALLOCATE( aux( nrxx ), gaux( nrxx ) )
  !
  gaux(:) = CMPLX(0.d0,0.d0, kind=dp)
  !
  DO ipol = 1, 3
     !
     aux = CMPLX( a(ipol,:), 0.D0 ,kind=DP)
     !
     ! ... bring a(ipol,r) to G-space, a(G) ...
     !
     CALL fwfft ('Dense', aux, dfftp)
     !
     DO n = 1, ngm
        !
        gaux(nl(n)) = gaux(nl(n)) + g(ipol,n) * &
                      CMPLX( -AIMAG( aux(nl(n)) ), REAL( aux(nl(n)) ) ,kind=DP)
        !
     END DO
    !
  END DO
  !
  IF ( gamma_only ) THEN
     !
     DO n = 1, ngm
        !
        gaux(nlm(n)) = CONJG( gaux(nl(n)) )
        !
     END DO
     !
  END IF
  !
  ! ... bring back to R-space, (\grad_ipol a)(r) ...
  !
  CALL invfft ('Dense', gaux, dfftp)
  !
  ! ... add the factor 2\pi/a  missing in the definition of G and sum
  !
  da(:) = tpiba * REAL( gaux(:) )
  !
  DEALLOCATE( aux, gaux )
  !
  RETURN
  !
END SUBROUTINE grad_dot
!--------------------------------------------------------------------
SUBROUTINE hessian( nrxx, a, ngm, g, nl, ga, ha )
!--------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a in R-space 
  ! ... and ha = \hessian a in R-space (a is also in R-space)
  !
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : nlm
  USE control_flags, ONLY : gamma_only
  USE fft_base,      ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nrxx
  INTEGER,  INTENT(IN)  :: ngm, nl(ngm)
  REAL(DP), INTENT(IN)  :: a(nrxx), g(3,ngm)
  REAL(DP), INTENT(OUT) :: ga( 3, nrxx )
  REAL(DP), INTENT(OUT) :: ha( 3, 3, nrxx )
  !
  INTEGER                  :: ipol, jpol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:), haux(:)
  !
  !
  ALLOCATE(  aux( nrxx ) )
  ALLOCATE( gaux( nrxx ) )
  ALLOCATE( haux( nrxx ) )
  !
  aux = CMPLX( a(:), 0.D0 ,kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Dense', aux, dfftp)
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  DO ipol = 1, 3
     !
     gaux(:) = CMPLX(0.d0,0.d0, kind=dp)
     !
     gaux(nl(:)) = g(ipol,:) * &
                   CMPLX( -AIMAG( aux(nl(:)) ), REAL( aux(nl(:)) ) ,kind=DP)
     !
     IF ( gamma_only ) THEN
        !
        gaux(nlm(:)) = CMPLX( REAL( gaux(nl(:)) ), -AIMAG( gaux(nl(:)) ) ,kind=DP)
        !
     END IF
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL invfft ('Dense', gaux, dfftp)
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     ga(ipol,:) = tpiba * DBLE( gaux(:) )
     !
     ! ... compute the second derivatives
     !
     DO jpol = 1, ipol
        !
        haux(:) = CMPLX(0.d0,0.d0, kind=dp)
        !
        haux(nl(:)) = - g(ipol,:) * g(jpol,:) * &
                       CMPLX( REAL( aux(nl(:)) ), AIMAG( aux(nl(:)) ) ,kind=DP)
        !
        IF ( gamma_only ) THEN
           !
           haux(nlm(:)) = CMPLX( REAL( haux(nl(:)) ), -AIMAG( haux(nl(:)) ) ,kind=DP)
           !
        END IF
        !
        ! ... bring back to R-space, (\grad_ipol a)(r) ...
        !
        CALL invfft ('Dense', haux, dfftp)
        !
        ! ...and add the factor 2\pi/a  missing in the definition of G
        !
        ha(ipol, jpol, :) = tpiba * tpiba * DBLE( haux(:) )
        !
        ha(jpol, ipol, :) = ha(ipol, jpol, :) 
        !
     END DO
     !
  END DO
  !
  DEALLOCATE( haux )
  DEALLOCATE( gaux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE hessian

!--------------------------------------------------------------------
SUBROUTINE ggradient( nrxx, a, ngm, g, nl, ga, gga )
!--------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a in R-space 
  ! ... and gga = \grad \grad a in R-space (a is also in R-space)
  !
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : nlm
  USE control_flags, ONLY : gamma_only
  USE fft_base,      ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nrxx
  INTEGER,  INTENT(IN)  :: ngm, nl(ngm)
  REAL(DP), INTENT(IN)  :: a(nrxx), g(3,ngm)
  REAL(DP), INTENT(OUT) :: ga( 3, nrxx )
  REAL(DP), INTENT(OUT) :: gga( 3, 3, nrxx )
  !
  INTEGER                  :: ipol, jpol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:), ggaux(:)
  !
  !
  ALLOCATE(  aux( nrxx ) )
  ALLOCATE( gaux( nrxx ) )
  ALLOCATE( ggaux( nrxx ) )
  !
  aux = CMPLX( a(:), 0.D0 ,kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Dense', aux, dfftp)
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  DO ipol = 1, 3
     !
     gaux(:) = CMPLX(0.d0,0.d0, kind=dp)
     !
     gaux(nl(:)) = g(ipol,:) * &
                   CMPLX( -AIMAG( aux(nl(:)) ), REAL( aux(nl(:)) ) ,kind=DP)
     !
     IF ( gamma_only ) THEN
        !
        gaux(nlm(:)) = CMPLX( REAL( gaux(nl(:)) ), -AIMAG( gaux(nl(:)) ) ,kind=DP)
        !
     END IF
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL invfft ('Dense', gaux, dfftp)
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     ga(ipol,:) = tpiba * DBLE( gaux(:) )
     !
     ! ... compute the second derivatives
     !
     DO jpol = 1, ipol
        !
        ggaux(:) = CMPLX(0.d0,0.d0, kind=dp)
        !
        ggaux(nl(:)) = - g(ipol,:) * g(jpol,:) * &
                       CMPLX( REAL( aux(nl(:)) ), AIMAG( aux(nl(:)) ) ,kind=DP)
        !
        IF ( gamma_only ) THEN
           !
           ggaux(nlm(:)) = CMPLX( REAL( ggaux(nl(:)) ), -AIMAG( ggaux(nl(:)) ) ,kind=DP)
           !
        END IF
        !
        ! ... bring back to R-space, (\grad_ipol a)(r) ...
        !
        CALL invfft ('Dense', ggaux, dfftp)
        !
        ! ...and add the factor 2\pi/a  missing in the definition of G
        !
        gga(ipol, jpol, :) = tpiba * tpiba * DBLE( ggaux(:) )
        !
        gga(jpol, ipol, :) = gga(ipol, jpol, :) 
        !
     END DO
     !
  END DO
  !
  DEALLOCATE( ggaux )
  DEALLOCATE( gaux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE ggradient

!--------------------------------------------------------------------
SUBROUTINE laplacian( nrxx, a, ngm, gg, nl, lapla )
!--------------------------------------------------------------------
  !
  ! ... Calculates lapla = \laplace a in R-space (a is also in R-space)
  !
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba2
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : nlm, gstart
  USE control_flags, ONLY : gamma_only
  USE fft_base,      ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nrxx
  INTEGER,  INTENT(IN)  :: ngm, nl(ngm)
  REAL(DP), INTENT(IN)  :: a(nrxx), gg(ngm)
  REAL(DP), INTENT(OUT) :: lapla( nrxx )
  !
  INTEGER                  :: ig
  COMPLEX(DP), ALLOCATABLE :: aux(:), laux(:)
  !
  !
  ALLOCATE(  aux( nrxx ) )
  ALLOCATE( laux( nrxx ) )
  !
  aux = CMPLX( a(:), 0.D0 ,kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Dense', aux, dfftp)
  !
  ! ... Compute the laplacian
  !
  laux(:) = CMPLX(0.d0,0.d0, kind=dp)
  !
  DO ig = gstart, ngm
     !
     laux(nl(ig)) = -gg(ig)*aux(nl(ig))
     !
  END DO
  !
  IF ( gamma_only ) THEN
     !
     laux(nlm(:)) = CMPLX( REAL(laux(nl(:)) ), -AIMAG(laux(nl(:)) ) ,kind=DP)
     !
  ENDIF
  !
  ! ... bring back to R-space, (\lapl a)(r) ...
  !
  CALL invfft ('Dense', laux, dfftp)
  !
  ! ... add the missing factor (2\pi/a)^2 in G
  !
  lapla = tpiba2 * DBLE( laux )   
  !
  DEALLOCATE( laux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE laplacian

!--------------------------------------------------------------------
SUBROUTINE external_gradient( a, grada )
!--------------------------------------------------------------------
  ! 
  ! Interface for computing gradients in real space, to be called by
  ! an external module
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm, nl, g
  !
  IMPLICIT NONE
  !
  REAL( DP ), INTENT(IN)   :: a( dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: grada( 3, dfftp%nnr )

! A in real space, grad(A) in real space
  CALL gradient( dfftp%nnr, a, ngm, g, nl, grada )

  RETURN

END SUBROUTINE external_gradient

!--------------------------------------------------------------------
SUBROUTINE external_ggradient( a, grada, ggrada )
!--------------------------------------------------------------------
  ! 
  ! Interface for computing gradient and hessian in real 
  ! space, to be called by an external module
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm, nl, g
  !
  IMPLICIT NONE
  !
  REAL( DP ), INTENT(IN)   :: a( dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: grada( 3, dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: ggrada( 3, 3, dfftp%nnr )

! A in real space, grad(A) in real space
  CALL ggradient( dfftp%nnr, a, ngm, g, nl, grada, ggrada )

  RETURN

END SUBROUTINE external_ggradient
!--------------------------------------------------------------------
SUBROUTINE external_hessian( a, grada, hessa )
!--------------------------------------------------------------------
  ! 
  ! Interface for computing hessian in real space, to be called by
  ! an external module
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm, nl, g
  !
  IMPLICIT NONE
  !
  REAL( DP ), INTENT(IN)   :: a( dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: grada( 3, dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: hessa( 3, 3, dfftp%nnr )

! A in real space, grad(A) and hess(A) in real space
  CALL hessian( dfftp%nnr, a, ngm, g, nl, grada, hessa )

  RETURN

END SUBROUTINE external_hessian
!--------------------------------------------------------------------
!--------------------------------------------------------------------
SUBROUTINE external_laplacian( a, lapla )
!--------------------------------------------------------------------
  ! 
  ! Interface for computing laplacian in real space, to be called by 
  ! an external module
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm, nl, gg
  !
  IMPLICIT NONE
  !
  REAL( DP ), INTENT(IN)   :: a( dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: lapla( dfftp%nnr )

! A in real space, lapl(A) in real space
  CALL laplacian( dfftp%nnr, a, ngm, gg, nl, lapla )

  RETURN

END SUBROUTINE external_laplacian
!--------------------------------------------------------------------

#ifdef USE_CUDA
module gradcorr_gpu_kernels_m
  use cudafor
  use funct_gpu
  contains
    attributes(global) subroutine gradcorr_gpu_upolar(nnr, igcx, igcc, e2, epsr, epsg, rhoout, rho_core, grho, v, h, vtxcgcs, etxcgcs)
      use kinds, only: DP
      implicit none
      integer, value :: nnr, igcx, igcc
      real(DP), value :: e2, epsr, epsg
      real(DP), device, intent(in) :: rhoout(nnr, 2), rho_core(nnr), grho(3, nnr, 1)
      real(DP), device, intent(out) :: v(nnr, 2), h(3, nnr, 1), vtxcgcs(nnr), etxcgcs(nnr)

      integer :: k
      real(DP) :: arho, segno, grho1, sx, sc, v1x, v2x, v1c, v2c
      real(DP) :: g1, g2, g3, r1

      k = (blockIdx%x - 1) * blockDim%x + threadIdx%x

      if (k > nnr) return

      r1 = rhoout(k,1)
      arho = ABS( r1 )
      !
      IF ( arho > epsr ) THEN
         !
         g1 = grho(1,k,1)
         g2 = grho(2,k,1)
         g3 = grho(3,k,1)

         grho1 = g1*g1 + g2*g2 + g3*g3
         !
         IF ( grho1 > epsg ) THEN
            !
            segno = SIGN( 1.D0, r1 )
            !
            CALL gcxc_dev( igcx, igcc, arho, grho1, sx, sc, v1x, v2x, v1c, v2c)
            !
            ! ... first term of the gradient correction : D(rho*Exc)/D(rho)
            !
            v(k,1) = v(k,1) + e2 * ( v1x + v1c )
            !
            ! ... h contains :
            !
            ! ...    D(rho*Exc) / D(|grad rho|) * (grad rho) / |grad rho|
            !
            h(1,k,1) = e2 * ( v2x + v2c ) * g1
            h(2,k,1) = e2 * ( v2x + v2c ) * g2
            h(3,k,1) = e2 * ( v2x + v2c ) * g3
            !
            vtxcgcs(k) = e2*( v1x + v1c ) * ( r1 - rho_core(k) )
            etxcgcs(k) = e2*( sx + sc ) * segno
            !
         ELSE
            h(1,k,1)=0.D0
            h(2,k,1)=0.D0
            h(3,k,1)=0.D0
            vtxcgcs(k) = 0.D0
            etxcgcs(k) = 0.D0
         END IF
         !
      ELSE
         !
         h(1,k,1) = 0.D0
         h(2,k,1) = 0.D0
         h(3,k,1) = 0.D0
         vtxcgcs(k) = 0.D0
         etxcgcs(k) = 0.D0
         !
      END IF
    end subroutine gradcorr_gpu_upolar

    attributes(global) subroutine gradcorr_gpu_polar(nnr, igcx, igcc, e2, epsr, fac, rhoout, rho_core, grho, v, h, vtxcgcs, etxcgcs)
      use kinds, only: DP
      implicit none
      integer, value :: nnr, igcx, igcc
      real(DP), value :: e2, epsr, fac
      real(DP), device, intent(in) :: rhoout(nnr, 2), rho_core(nnr), grho(3, nnr, 2)
      real(DP), device, intent(out) :: v(nnr, 2), h(3, nnr, 2), vtxcgcs(nnr), etxcgcs(nnr)

      integer :: k, ipol
      real(DP) :: zeta, arho, segno, grho1, grho2, sx, sc, v1xup, v1xdw, v2xup, v2xdw
      real(DP) :: v1cup, v1cdw, v2cup, v2cdw, v2cud, v2c
      real(DP) :: rh, grh2, grup, grdw, r1, r2
      real(DP) :: g1(3), g2(3)

      k = (blockIdx%x - 1) * blockDim%x + threadIdx%x

      if (k > nnr) return

      r1 = rhoout(k,1)
      r2 = rhoout(k,2)

      rh = r1 + r2

      g1(1) = grho(1,k,1)
      g1(2) = grho(2,k,1)
      g1(3) = grho(3,k,1)
      g2(1) = grho(1,k,2)
      g2(2) = grho(2,k,2)
      g2(3) = grho(3,k,2)
      !
      !grho2(:) = grho(1,k,:)**2 + grho(2,k,:)**2 + grho(3,k,:)**2
      grho1 = g1(1)**2 + g1(2)**2 + g1(3)**2
      grho2 = g2(1)**2 + g2(2)**2 + g2(3)**2
      !
      CALL gcx_spin_dev( igcx, r1, r2, grho1, &
                     grho2, sx, v1xup, v1xdw, v2xup, v2xdw )
      !
      IF ( rh > epsr ) THEN
         !
         !IF ( igcc_is_lyp() ) THEN
         !   print*, "gradcorr_gpu: igcc_is_lyp() not implemented!"; flush(6); stop
            !
            !rup = rhoout(k,1)
            !rdw = rhoout(k,2)
            !!
            !grhoup = grho(1,k,1)**2 + grho(2,k,1)**2 + grho(3,k,1)**2
            !grhodw = grho(1,k,2)**2 + grho(2,k,2)**2 + grho(3,k,2)**2
            !!
            !grhoud = grho(1,k,1) * grho(1,k,2) + &
            !         grho(2,k,1) * grho(2,k,2) + &
            !         grho(3,k,1) * grho(3,k,2)
            !!
            !CALL gcc_spin_more( rup, rdw, grhoup, grhodw, grhoud, &
            !                    sc, v1cup, v1cdw, v2cup, v2cdw, v2cud )
            !
         !ELSE
            !
            zeta = ( r1 - r2 ) / rh
            !if (nspin.eq.4.and.domag) zeta=abs(zeta)*segni(k)
            !
            grh2 = ( g1(1) + g2(1) )**2 + &
                   ( g1(2) + g2(2) )**2 + &
                   ( g1(3) + g2(3) )**2
            !
            CALL gcc_spin_dev( igcc, rh, zeta, grh2, sc, v1cup, v1cdw, v2c )
            !
            v2cup = v2c
            v2cdw = v2c
            v2cud = v2c
            !
         !END IF
         !
      ELSE
         !
         sc    = 0.D0
         v1cup = 0.D0
         v1cdw = 0.D0
         v2c   = 0.D0
         v2cup = 0.D0
         v2cdw = 0.D0
         v2cud = 0.D0
         !
      ENDIF
      !
      ! ... first term of the gradient correction : D(rho*Exc)/D(rho)
      !
      v(k,1) = v(k,1) + e2 * ( v1xup + v1cup )
      v(k,2) = v(k,2) + e2 * ( v1xdw + v1cdw )
      !
      ! ... h contains D(rho*Exc)/D(|grad rho|) * (grad rho) / |grad rho|
      !
      DO ipol = 1, 3
         !
         grup = g1(ipol)
         grdw = g2(ipol)
         h(ipol,k,1) = e2 * ( ( v2xup + v2cup ) * grup + v2cud * grdw )
         h(ipol,k,2) = e2 * ( ( v2xdw + v2cdw ) * grdw + v2cud * grup )
         !
      END DO
      !
      vtxcgcs(k) = e2 * ( v1xup + v1cup ) * ( r1 - rho_core(k) * fac ) + &
                    e2 * ( v1xdw + v1cdw ) * ( r2 - rho_core(k) * fac )
      etxcgcs(k) =  e2 * ( sx + sc)
      !

    end subroutine gradcorr_gpu_polar

end module gradcorr_gpu_kernels_m

SUBROUTINE gradcorr_gpu(rho_d, rhog_d, rho_core_d, rhog_core_d, etxc, vtxc, v_d )
  !----------------------------------------------------------------------------
  !
  USE constants,            ONLY : e2, pi
  USE kinds,                ONLY : DP
  USE gvect,                ONLY : nl, nl_d, ngm, g, g_d
  USE lsda_mod,             ONLY : nspin
  USE cell_base,            ONLY : omega, alat
  USE funct,                ONLY : gcxc, gcx_spin, gcc_spin, igcc_is_lyp, &
                                   gcc_spin_more, dft_is_gradient, get_igcc, get_igcx
  USE spin_orb,             ONLY : domag
  USE noncollin_module,     ONLY : ux
  USE wavefunctions_module, ONLY : psic
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : fwfft
  USE funct_gpu
  USE gradcorr_gpu_kernels_m

  !
  IMPLICIT NONE
  !
  REAL(DP),    INTENT(IN), device    :: rho_d(dfftp%nnr,nspin), rho_core_d(dfftp%nnr)
  COMPLEX(DP), INTENT(IN), device    :: rhog_d(ngm,nspin), rhog_core_d(ngm)
  REAL(DP),    INTENT(INOUT), device :: v_d(dfftp%nnr,nspin)
  REAL(DP),    INTENT(INOUT) :: vtxc, etxc

  !
  INTEGER :: k, ipol, is, nspin0, ir, jpol
  !
  REAL(DP),    ALLOCATABLE :: segni(:), vgg(:,:), vsave(:,:)
  REAL(DP),    ALLOCATABLE :: gmag(:,:,:)

  REAL(DP),    ALLOCATABLE, device :: grho_d(:,:,:), h_d(:,:,:), dh_d(:), rhoout_d(:,:)

  COMPLEX(DP), ALLOCATABLE, device :: rhogsum_d(:,:)
  !
  REAL(DP) :: grho1, grho2(2), etxcgc, vtxcgc, segno, fac, amag 

  REAL(DP), device, allocatable :: etxcs_d(:), vtxcs_d(:)
  !
  REAL(DP), PARAMETER :: epsr = 1.D-6, epsg = 1.D-10
  integer :: i, j, igcx, igcc
  integer :: threads, blocks
  !
  !
  IF ( .NOT. dft_is_gradient() ) RETURN
  !
  call start_clock( 'gradcorr' )
  !
  etxcgc = 0.D0
  vtxcgc = 0.D0
  !
  nspin0=nspin
  if (nspin==4) nspin0=1
  if (nspin==4.and.domag) nspin0=2
  fac = 1.D0 / DBLE( nspin0 )
  !
  ALLOCATE(    h_d( 3, dfftp%nnr, nspin0) )
  ALLOCATE( grho_d( 3, dfftp%nnr, nspin0) )
  ALLOCATE( rhoout_d( dfftp%nnr, nspin0) )
!  IF (nspin==4.AND.domag) THEN
!     ALLOCATE( vgg( dfftp%nnr, nspin0 ) )
!     ALLOCATE( vsave( dfftp%nnr, nspin ) )
!     ALLOCATE( segni( dfftp%nnr ) )
!     vsave=v
!     v=0.d0
!  ENDIF
  !
  ALLOCATE( rhogsum_d( ngm, nspin0 ) )
  !
  ! ... calculate the gradient of rho + rho_core in real space
  !
  IF ( nspin == 4 .AND. domag ) THEN
     print*, "gradcorr_gpu: nspin == 4 .and. domag not supported!"; flush(6); stop
!     !
!     CALL compute_rho(rho,rhoout,segni,dfftp%nnr)
!     !
!     ! ... bring starting rhoout to G-space
!     !
!     DO is = 1, nspin0
!        !
!        psic(:) = rhoout(:,is)
!        !
!        CALL fwfft ('Dense', psic, dfftp)
!        !
!        rhogsum(:,is) = psic(nl(:))
!        !
!     END DO
  ELSE
     !
     !$cuf kernel do(2) <<<*,*>>>
     do j = 1, nspin0
       do i = lbound(rhoout_d,1), ubound(rhoout_d,1)
         rhoout_d(i,j)  = rho_d(i,j)
       end do
     end do

     !$cuf kernel do(2) <<<*,*>>>
     do j = 1, nspin0
       do i = lbound(rhogsum_d,1), ubound(rhogsum_d,1)
         rhogsum_d(i,j)  = rhog_d(i,j)
       end do
     end do
     !
  ENDIF

  DO is = 1, nspin0
     !
     !$cuf kernel do(1) <<<*,*>>>
     do i = lbound(rhoout_d,1), ubound(rhoout_d,1)
       rhoout_d(i,is)  = fac * rho_core_d(i)  + rhoout_d(i,is)
     end do
     !$cuf kernel do(1) <<<*,*>>>
     do i = lbound(rhogsum_d,1), ubound(rhogsum_d,1)
       rhogsum_d(i,is) = fac * rhog_core_d(i) + rhogsum_d(i,is)
     end do
     !
call start_clock( 'gradrho' )
     CALL gradrho_gpu( dfftp%nnr, rhogsum_d(1,is), ngm, g_d, nl_d, grho_d(1,1,is) )
call stop_clock( 'gradrho' )
     
     !
  END DO
  !
  DEALLOCATE( rhogsum_d )
  !
  igcx = get_igcx
  igcc = get_igcc

  IF ( nspin0 == 1 ) THEN
    !
    !
    ! ... This is the spin-unpolarised case
    !
    allocate(etxcs_d(dfftp%nnr))
    allocate(vtxcs_d(dfftp%nnr))

    threads = 128
    blocks = ceiling(real(dfftp%nnr)/threads)
    call gradcorr_gpu_upolar<<<blocks, threads>>>(dfftp%nnr, igcx, igcc, e2, epsr, epsg, rhoout_d, &
      rho_core_d, grho_d, v_d, h_d, vtxcs_d, etxcs_d)

    ! Complete reductions after kernel
    !$cuf kernel do(1) <<<*, *>>>
    DO k = 1, dfftp%nnr
      vtxcgc = vtxcgc + vtxcs_d(k)
      etxcgc = etxcgc + etxcs_d(k)
    ENDDO

    deallocate(etxcs_d)
    deallocate(vtxcs_d)

  ELSE
     !
     ! ... spin-polarised case
     !
      allocate(etxcs_d(dfftp%nnr))
      allocate(vtxcs_d(dfftp%nnr))

      threads = 128
      blocks = ceiling(real(dfftp%nnr)/threads)
      call gradcorr_gpu_polar<<<blocks, threads>>>(dfftp%nnr, igcx, igcc, e2, epsr, fac, &
        rhoout_d, rho_core_d, grho_d, v_d, h_d, vtxcs_d, etxcs_d)

      ! Complete reductions after kernel
      !$cuf kernel do(1) <<<*, *>>>
      DO k = 1, dfftp%nnr
        vtxcgc = vtxcgc + vtxcs_d(k)
        etxcgc = etxcgc + etxcs_d(k)
      ENDDO

      deallocate(etxcs_d)
      deallocate(vtxcs_d)
     
  END IF
  !
  !$cuf kernel do(2)
  do is = 1, nspin0
    do i = lbound(rhoout_d, 1), ubound(rhoout_d,1)
      rhoout_d(i,is) = rhoout_d(i,is) - fac * rho_core_d(i)
    end do
  end do
  !
  DEALLOCATE( grho_d )
  !
  ALLOCATE( dh_d( dfftp%nnr ) )    
  !
  ! ... second term of the gradient correction :
  ! ... \sum_alpha (D / D r_alpha) ( D(rho*Exc)/D(grad_alpha rho) )
  !
  DO is = 1, nspin0
     !
call start_clock( 'graddot' )
     CALL grad_dot_gpu( dfftp%nnr, h_d(1,1,is), ngm, g_d, nl_d, alat, dh_d )
call stop_clock( 'graddot' )
     !
     !$cuf kernel do(1)
     do i = lbound(v_d, 1), ubound(v_d,1)
       v_d(i,is) = v_d(i,is) - dh_d(i)
     end do
     !
     !$cuf kernel do(1)
     do i = lbound(dh_d, 1), ubound(dh_d,1)
       vtxcgc = vtxcgc - dh_d(i) * rhoout_d(i,is)
     end do
     !
  END DO

  !
  vtxc = vtxc + omega * vtxcgc / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  etxc = etxc + omega * etxcgc / ( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )

!  IF (nspin==4.AND.domag) THEN
!     DO is=1,nspin0
!        vgg(:,is)=v(:,is)
!     ENDDO
!     v=vsave
!     DO k=1,dfftp%nnr
!        v(k,1)=v(k,1)+0.5d0*(vgg(k,1)+vgg(k,2))
!        amag=sqrt(rho(k,2)**2+rho(k,3)**2+rho(k,4)**2)
!        IF (amag.GT.1.d-12) THEN
!           v(k,2)=v(k,2)+segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,2)/amag
!           v(k,3)=v(k,3)+segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,3)/amag
!           v(k,4)=v(k,4)+segni(k)*0.5d0*(vgg(k,1)-vgg(k,2))*rho(k,4)/amag
!        ENDIF
!     ENDDO
!  ENDIF
  !
  DEALLOCATE( dh_d )
  DEALLOCATE( h_d )
  DEALLOCATE( rhoout_d )
  IF (nspin==4.and.domag) THEN
     DEALLOCATE( vgg )
     DEALLOCATE( vsave )
     DEALLOCATE( segni )
  ENDIF
  !
  call stop_clock( 'gradcorr' )
  RETURN
  !
END SUBROUTINE gradcorr_gpu
!
!----------------------------------------------------------------------------
SUBROUTINE gradrho_gpu( nrxx, a, ngm, g, nl, ga )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a in R-space (a is in G-space)
  !
  USE kinds,     ONLY : DP
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE gvect,     ONLY : nlm
  USE control_flags, ONLY : gamma_only
  USE fft_base,      ONLY : dfftp
  USE fft_interfaces,ONLY : invfft

  !
  IMPLICIT NONE
  !
  INTEGER,     INTENT(IN)  :: nrxx
  INTEGER,     INTENT(IN)  :: ngm
  INTEGER,     INTENT(IN), device  :: nl(ngm)
  COMPLEX(DP), INTENT(IN), device  :: a(ngm)
  REAL(DP),    INTENT(IN), device  :: g(3,ngm)
  REAL(DP),    INTENT(OUT), device :: ga(3,nrxx)
  !
  INTEGER                  :: i, ipol
  COMPLEX(DP), ALLOCATABLE, device :: gaux_d(:)
  !
  !
  ALLOCATE( gaux_d (nrxx) )
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  ga(:,:) = 0.D0
  !
  DO ipol = 1, 3
     !
     gaux_d(:) = CMPLX(0.d0,0.d0,kind=dp)
     !
     !$cuf kernel do(1) <<<*,*>>>
     do i = lbound(a,1), ubound(a,1)
       gaux_d(nl(i)) = g(ipol,i) * CMPLX( -AIMAG( a(i) ), REAL( a(i) ) ,kind=DP)
     end do
     !
     IF ( gamma_only ) THEN
        print*, "gradrho_gpu: GAMMA_ONLY not implemented!"
        flush(6); stop
        !
        !gaux(nlm(:)) = CMPLX( REAL( gaux(nl(:)) ), -AIMAG( gaux(nl(:)) ) ,kind=DP)
        !
     END IF
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL invfft ('Dense', gaux_d, dfftp)
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     !$cuf kernel do(1) <<<*,*>>>
     do i = lbound(ga,2), ubound(ga,2)
       ga(ipol,i) = ga(ipol,i) + tpiba * REAL( gaux_d(i) )
     end do
     !
  END DO
  !
  DEALLOCATE( gaux_d )
  !
  RETURN
  !
END SUBROUTINE gradrho_gpu
!
!----------------------------------------------------------------------------
SUBROUTINE gradient_gpu( nrxx, a, ngm, g, nl, ga )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a in R-space (a is also in R-space)
  !
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : nlm
  USE control_flags, ONLY : gamma_only
  USE fft_base,      ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nrxx
  INTEGER,  INTENT(IN)  :: ngm, nl(ngm)
  REAL(DP), INTENT(IN)  :: a(nrxx), g(3,ngm)
  REAL(DP), INTENT(OUT) :: ga(3,nrxx)
  !
  INTEGER                  :: ipol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:)
  !
  !
  ALLOCATE(  aux( nrxx ) )
  ALLOCATE( gaux( nrxx ) )
  !
  aux = CMPLX( a(:), 0.D0 ,kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Dense', aux, dfftp)
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  DO ipol = 1, 3
     !
     gaux(:) = CMPLX(0.d0,0.d0, kind=dp)
     !
     gaux(nl(:)) = g(ipol,:) * &
                   CMPLX( -AIMAG( aux(nl(:)) ), REAL( aux(nl(:)) ) ,kind=DP)
     !
     IF ( gamma_only ) THEN
        !
        gaux(nlm(:)) = CMPLX( REAL( gaux(nl(:)) ), -AIMAG( gaux(nl(:)) ) ,kind=DP)
        !
     END IF
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL invfft ('Dense', gaux, dfftp)
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     ga(ipol,:) = tpiba * DBLE( gaux(:) )
     !
  END DO
  !
  DEALLOCATE( gaux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE gradient_gpu
!
!----------------------------------------------------------------------------
SUBROUTINE grad_dot_gpu( nrxx, a, ngm, g, nl, alat, da )
  !----------------------------------------------------------------------------
  !
  ! ... Calculates da = \sum_i \grad_i a_i in R-space
  !
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : nlm
  USE control_flags, ONLY : gamma_only
  USE fft_base,      ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)     :: nrxx, ngm
  INTEGER,  INTENT(IN), device     :: nl(ngm)
  REAL(DP), INTENT(IN)     ::  alat
  REAL(DP), INTENT(IN), device     :: a(3,nrxx), g(3,ngm)
  REAL(DP), INTENT(OUT), device    :: da(nrxx)
  !
  INTEGER                  :: n, ipol, i
  COMPLEX(DP), ALLOCATABLE, device :: aux_d(:), gaux_d(:)
  !
  !
  ALLOCATE( aux_d( nrxx ), gaux_d( nrxx ) )
  !
  gaux_d(:) = CMPLX(0.d0,0.d0, kind=dp)
  !
  DO ipol = 1, 3
     !
     !$cuf kernel do(1) <<<*,*>>>
     do i = lbound(a,2), ubound(a,2)
       aux_d(i) = CMPLX( a(ipol,i), 0.D0 ,kind=DP)
     end do
     !
     ! ... bring a(ipol,r) to G-space, a(G) ...
     !
     CALL fwfft ('Dense', aux_d, dfftp)
     !
     !$cuf kernel do(1) <<<*,*>>>
     DO n = 1, ngm
        !
        gaux_d(nl(n)) = gaux_d(nl(n)) + g(ipol,n) * &
                      CMPLX( -AIMAG( aux_d(nl(n)) ), REAL( aux_d(nl(n)) ) ,kind=DP)
        !
     END DO
    !
  END DO
  !
  IF ( gamma_only ) THEN
     print*, "grad_dot_gpu: GAMMA_ONLY not implemented!!"
     flush(6); stop
     !
     !DO n = 1, ngm
     !   !
     !   gaux(nlm(n)) = CONJG( gaux(nl(n)) )
     !   !
     !END DO
     !
  END IF
  !
  ! ... bring back to R-space, (\grad_ipol a)(r) ...
  !
  CALL invfft ('Dense', gaux_d, dfftp)
  !
  ! ... add the factor 2\pi/a  missing in the definition of G and sum
  !
  !$cuf kernel do(1) <<<*,*>>>
  do i = lbound(da,1), ubound(da,1)
    da(i) = tpiba * REAL( gaux_d(i) )
  end do
  !
  DEALLOCATE( aux_d, gaux_d )
  !
  RETURN
  !
END SUBROUTINE grad_dot_gpu
!--------------------------------------------------------------------
SUBROUTINE hessian_gpu( nrxx, a, ngm, g, nl, ga, ha )
!--------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a in R-space 
  ! ... and ha = \hessian a in R-space (a is also in R-space)
  !
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : nlm
  USE control_flags, ONLY : gamma_only
  USE fft_base,      ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nrxx
  INTEGER,  INTENT(IN)  :: ngm, nl(ngm)
  REAL(DP), INTENT(IN)  :: a(nrxx), g(3,ngm)
  REAL(DP), INTENT(OUT) :: ga( 3, nrxx )
  REAL(DP), INTENT(OUT) :: ha( 3, 3, nrxx )
  !
  INTEGER                  :: ipol, jpol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:), haux(:)
  !
  !
  ALLOCATE(  aux( nrxx ) )
  ALLOCATE( gaux( nrxx ) )
  ALLOCATE( haux( nrxx ) )
  !
  aux = CMPLX( a(:), 0.D0 ,kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Dense', aux, dfftp)
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  DO ipol = 1, 3
     !
     gaux(:) = CMPLX(0.d0,0.d0, kind=dp)
     !
     gaux(nl(:)) = g(ipol,:) * &
                   CMPLX( -AIMAG( aux(nl(:)) ), REAL( aux(nl(:)) ) ,kind=DP)
     !
     IF ( gamma_only ) THEN
        !
        gaux(nlm(:)) = CMPLX( REAL( gaux(nl(:)) ), -AIMAG( gaux(nl(:)) ) ,kind=DP)
        !
     END IF
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL invfft ('Dense', gaux, dfftp)
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     ga(ipol,:) = tpiba * DBLE( gaux(:) )
     !
     ! ... compute the second derivatives
     !
     DO jpol = 1, ipol
        !
        haux(:) = CMPLX(0.d0,0.d0, kind=dp)
        !
        haux(nl(:)) = - g(ipol,:) * g(jpol,:) * &
                       CMPLX( REAL( aux(nl(:)) ), AIMAG( aux(nl(:)) ) ,kind=DP)
        !
        IF ( gamma_only ) THEN
           !
           haux(nlm(:)) = CMPLX( REAL( haux(nl(:)) ), -AIMAG( haux(nl(:)) ) ,kind=DP)
           !
        END IF
        !
        ! ... bring back to R-space, (\grad_ipol a)(r) ...
        !
        CALL invfft ('Dense', haux, dfftp)
        !
        ! ...and add the factor 2\pi/a  missing in the definition of G
        !
        ha(ipol, jpol, :) = tpiba * tpiba * DBLE( haux(:) )
        !
        ha(jpol, ipol, :) = ha(ipol, jpol, :) 
        !
     END DO
     !
  END DO
  !
  DEALLOCATE( haux )
  DEALLOCATE( gaux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE hessian_gpu

!--------------------------------------------------------------------
SUBROUTINE ggradient_gpu( nrxx, a, ngm, g, nl, ga, gga )
!--------------------------------------------------------------------
  !
  ! ... Calculates ga = \grad a in R-space 
  ! ... and gga = \grad \grad a in R-space (a is also in R-space)
  !
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : nlm
  USE control_flags, ONLY : gamma_only
  USE fft_base,      ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nrxx
  INTEGER,  INTENT(IN)  :: ngm, nl(ngm)
  REAL(DP), INTENT(IN)  :: a(nrxx), g(3,ngm)
  REAL(DP), INTENT(OUT) :: ga( 3, nrxx )
  REAL(DP), INTENT(OUT) :: gga( 3, 3, nrxx )
  !
  INTEGER                  :: ipol, jpol
  COMPLEX(DP), ALLOCATABLE :: aux(:), gaux(:), ggaux(:)
  !
  !
  ALLOCATE(  aux( nrxx ) )
  ALLOCATE( gaux( nrxx ) )
  ALLOCATE( ggaux( nrxx ) )
  !
  aux = CMPLX( a(:), 0.D0 ,kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Dense', aux, dfftp)
  !
  ! ... multiply by (iG) to get (\grad_ipol a)(G) ...
  !
  DO ipol = 1, 3
     !
     gaux(:) = CMPLX(0.d0,0.d0, kind=dp)
     !
     gaux(nl(:)) = g(ipol,:) * &
                   CMPLX( -AIMAG( aux(nl(:)) ), REAL( aux(nl(:)) ) ,kind=DP)
     !
     IF ( gamma_only ) THEN
        !
        gaux(nlm(:)) = CMPLX( REAL( gaux(nl(:)) ), -AIMAG( gaux(nl(:)) ) ,kind=DP)
        !
     END IF
     !
     ! ... bring back to R-space, (\grad_ipol a)(r) ...
     !
     CALL invfft ('Dense', gaux, dfftp)
     !
     ! ...and add the factor 2\pi/a  missing in the definition of G
     !
     ga(ipol,:) = tpiba * DBLE( gaux(:) )
     !
     ! ... compute the second derivatives
     !
     DO jpol = 1, ipol
        !
        ggaux(:) = CMPLX(0.d0,0.d0, kind=dp)
        !
        ggaux(nl(:)) = - g(ipol,:) * g(jpol,:) * &
                       CMPLX( REAL( aux(nl(:)) ), AIMAG( aux(nl(:)) ) ,kind=DP)
        !
        IF ( gamma_only ) THEN
           !
           ggaux(nlm(:)) = CMPLX( REAL( ggaux(nl(:)) ), -AIMAG( ggaux(nl(:)) ) ,kind=DP)
           !
        END IF
        !
        ! ... bring back to R-space, (\grad_ipol a)(r) ...
        !
        CALL invfft ('Dense', ggaux, dfftp)
        !
        ! ...and add the factor 2\pi/a  missing in the definition of G
        !
        gga(ipol, jpol, :) = tpiba * tpiba * DBLE( ggaux(:) )
        !
        gga(jpol, ipol, :) = gga(ipol, jpol, :) 
        !
     END DO
     !
  END DO
  !
  DEALLOCATE( ggaux )
  DEALLOCATE( gaux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE ggradient_gpu

!--------------------------------------------------------------------
SUBROUTINE laplacian_gpu( nrxx, a, ngm, gg, nl, lapla )
!--------------------------------------------------------------------
  !
  ! ... Calculates lapla = \laplace a in R-space (a is also in R-space)
  !
  USE constants, ONLY : tpi
  USE cell_base, ONLY : tpiba2
  USE kinds,     ONLY : DP
  USE gvect,     ONLY : nlm, gstart
  USE control_flags, ONLY : gamma_only
  USE fft_base,      ONLY : dfftp
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  !
  INTEGER,  INTENT(IN)  :: nrxx
  INTEGER,  INTENT(IN)  :: ngm, nl(ngm)
  REAL(DP), INTENT(IN)  :: a(nrxx), gg(ngm)
  REAL(DP), INTENT(OUT) :: lapla( nrxx )
  !
  INTEGER                  :: ig
  COMPLEX(DP), ALLOCATABLE :: aux(:), laux(:)
  !
  !
  ALLOCATE(  aux( nrxx ) )
  ALLOCATE( laux( nrxx ) )
  !
  aux = CMPLX( a(:), 0.D0 ,kind=DP)
  !
  ! ... bring a(r) to G-space, a(G) ...
  !
  CALL fwfft ('Dense', aux, dfftp)
  !
  ! ... Compute the laplacian
  !
  laux(:) = CMPLX(0.d0,0.d0, kind=dp)
  !
  DO ig = gstart, ngm
     !
     laux(nl(ig)) = -gg(ig)*aux(nl(ig))
     !
  END DO
  !
  IF ( gamma_only ) THEN
     !
     laux(nlm(:)) = CMPLX( REAL(laux(nl(:)) ), -AIMAG(laux(nl(:)) ) ,kind=DP)
     !
  ENDIF
  !
  ! ... bring back to R-space, (\lapl a)(r) ...
  !
  CALL invfft ('Dense', laux, dfftp)
  !
  ! ... add the missing factor (2\pi/a)^2 in G
  !
  lapla = tpiba2 * DBLE( laux )   
  !
  DEALLOCATE( laux )
  DEALLOCATE( aux )
  !
  RETURN
  !
END SUBROUTINE laplacian_gpu

!--------------------------------------------------------------------
SUBROUTINE external_gradient_gpu( a, grada )
!--------------------------------------------------------------------
  ! 
  ! Interface for computing gradients in real space, to be called by
  ! an external module
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm, nl, g
  !
  IMPLICIT NONE
  !
  REAL( DP ), INTENT(IN)   :: a( dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: grada( 3, dfftp%nnr )

! A in real space, grad(A) in real space
  CALL gradient( dfftp%nnr, a, ngm, g, nl, grada )

  RETURN

END SUBROUTINE external_gradient_gpu

!--------------------------------------------------------------------
SUBROUTINE external_ggradient_gpu( a, grada, ggrada )
!--------------------------------------------------------------------
  ! 
  ! Interface for computing gradient and hessian in real 
  ! space, to be called by an external module
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm, nl, g
  !
  IMPLICIT NONE
  !
  REAL( DP ), INTENT(IN)   :: a( dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: grada( 3, dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: ggrada( 3, 3, dfftp%nnr )

! A in real space, grad(A) in real space
  CALL ggradient( dfftp%nnr, a, ngm, g, nl, grada, ggrada )

  RETURN

END SUBROUTINE external_ggradient_gpu
!--------------------------------------------------------------------
SUBROUTINE external_hessian_gpu( a, grada, hessa )
!--------------------------------------------------------------------
  ! 
  ! Interface for computing hessian in real space, to be called by
  ! an external module
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm, nl, g
  !
  IMPLICIT NONE
  !
  REAL( DP ), INTENT(IN)   :: a( dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: grada( 3, dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: hessa( 3, 3, dfftp%nnr )

! A in real space, grad(A) and hess(A) in real space
  CALL hessian( dfftp%nnr, a, ngm, g, nl, grada, hessa )

  RETURN

END SUBROUTINE external_hessian_gpu
!--------------------------------------------------------------------
!--------------------------------------------------------------------
SUBROUTINE external_laplacian_gpu( a, lapla )
!--------------------------------------------------------------------
  ! 
  ! Interface for computing laplacian in real space, to be called by 
  ! an external module
  !
  USE kinds,            ONLY : DP
  USE fft_base,         ONLY : dfftp
  USE gvect,            ONLY : ngm, nl, gg
  !
  IMPLICIT NONE
  !
  REAL( DP ), INTENT(IN)   :: a( dfftp%nnr )
  REAL( DP ), INTENT(OUT)  :: lapla( dfftp%nnr )

! A in real space, lapl(A) in real space
  CALL laplacian( dfftp%nnr, a, ngm, gg, nl, lapla )

  RETURN

END SUBROUTINE external_laplacian_gpu
!--------------------------------------------------------------------

#endif
