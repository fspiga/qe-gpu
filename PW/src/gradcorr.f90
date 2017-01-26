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

subroutine gcxc_gpu (rho, grho, sx, sc, v1x, v2x, v1c, v2c)
  use kinds, ONLY : DP
  use constants, ONLY : pi
  implicit none
  real(DP) :: rho, grho, sx, sc, v1x, v2x, v1c, v2c

  !locals
  real(DP), parameter:: small = 1.E-10_DP
  real(DP) :: kf, agrho, s1, s2, ds, dsg, exunif, fx
  real(DP) :: dxunif, dfx, f1, f2, f3, dfx1
  real(DP) :: p,  amu, ab, c, dfxdp, dfxds, upbe, uge, s, ak, aa
  real(DP), parameter :: third = 1._DP / 3._DP, c1 = 0.75_DP / pi , &
       c2 = 3.093667726280136_DP, c5 = 4._DP * third, &
       c6 = c2*2.51984210, c7=5._DP/6._DP, c8=0.8_DP, & ! (3pi^2)^(1/3)*2^(4/3)
       k1 = 0.804_DP, mu1 = 0.21951_DP 
!  real(DP) :: k (6), mu(6), ev(6)
  !           pbe        rpbe        pbesol   pbeq2d      optB88  optB86b
!  data k / 0.804_DP,   1.2450D0,   0.804_DP , 0.804_DP ,    0.0 ,  0.0 /, &
!       mu/ 0.21951_DP, 0.21951_DP, 0.12345679012345679012_DP,             &
!                                   0.12345679012345679,     0.22 , 0.1234/, &
!       ev / 1.647127_DP, 0.980118_DP, 0.017399_DP, 1.523671_DP, 0.367229_DP, &
!                                   0.011282_DP /  ! a and b parameters of Engel and Vosko
  !
  real(DP), parameter :: ga = 0.031091d0, be1 = 0.066725d0
!  real(DP) :: be (3)
!  data be / 0.066725d0, 0.046d0,     0.066725d0/
  real(DP), parameter :: pi34 = 0.6203504908994d0
  real(DP), parameter :: xkf = 1.919158292677513d0, xks = 1.128379167095513d0
  real(DP) :: ks, rs, ec, vc, t, expe, af, bf, y, xy, qy
  real(DP) :: h0, dh0, ddh0, sc2D, v1c2D, v2c2D
  !
  real(DP) :: a, b1, b2, c0, d0, d1
  parameter (a = 0.031091d0, b1 = 7.5957d0, b2 = 3.5876d0, c0 = a)
  real(DP) :: lnrs, rs12, rs32, rs2, om, dom, olog
  real(DP), parameter :: a11 = 0.21370d0, b31 = 1.6382d0, b41 = 0.49294d0
!  real(DP) :: a1 (2), b3 (2), b4 (2)
!  data a1 / 0.21370d0, 0.026481d0 /, b3 / 1.6382d0, -0.46647d0 /, &
!       b4 / 0.49294d0, 0.13354d0 /
  !
  if (rho <= small) then
     sx = 0.0_DP
     v1x = 0.0_DP
     v2x = 0.0_DP
  else
!     call pbex (rho, grho, 1, sx, v1x, v2x)
     agrho = sqrt (grho)
     kf = c2 * rho**third
     dsg = 0.5_DP / kf
     s1 = agrho * dsg / rho
     s2 = s1 * s1
     ds = - c5 * s1
     f1 = s2 * mu1 / k1
     f2 = 1._DP + f1
     f3 = k1 / f2
     fx = k1 - f3
     exunif = - c1 * kf
     sx = exunif * fx
     dxunif = exunif * third
     dfx1 = f2 * f2
     dfx = 2._DP * mu1 * s1 / dfx1
     v1x = sx + dxunif * fx + exunif * dfx * ds
     v2x = exunif * dfx * dsg / agrho
     sx = sx * rho
  endif

  if (rho.le.small) then
     sc = 0.0_DP
     v1c = 0.0_DP
     v2c = 0.0_DP
  else
     !call pbec (rho, grho, 1, sc, v1c, v2c)
     rs = pi34 / rho**third
     !
     !call pw_gpu (rs, ec, vc)
     rs12 = sqrt (rs)
     rs32 = rs * rs12
     rs2 = rs**2
     om = 2.d0 * a * (b1 * rs12 + b2 * rs + b31 * rs32 + b41 * rs2)
     dom = 2.d0 * a * (0.5d0 * b1 * rs12 + b2 * rs + 1.5d0 * b31 * rs32 + 2.d0 * b41 * rs2)
     olog = log (1.d0 + 1.0d0 / om)
     ec = - 2.d0 * a * (1.d0 + a11 * rs) * olog
     vc = - 2.d0 * a * (1.d0 + 2.d0 / 3.d0 * a11 * rs) &
          * olog - 2.d0 / 3.d0 * a * (1.d0 + a11 * rs) * dom / (om * (om + 1.d0) )
     !
     kf = xkf / rs
     ks = xks * sqrt (kf)
     t = sqrt (grho) / (2.d0 * ks * rho)
     expe = exp ( - ec / ga)
     af = be1 / ga * (1.d0 / (expe-1.d0) )
     bf = expe * (vc - ec)
     y = af * t * t
     xy = (1.d0 + y) / (1.d0 + y + y * y)
     qy = y * y * (2.d0 + y) / (1.d0 + y + y * y) **2
     s1 = 1.d0 + be1 / ga * t * t * xy
     h0 = ga * log (s1)
     dh0 = be1 * t * t / s1 * ( - 7.d0 / 3.d0 * xy - qy * (af * bf / &
          be1-7.d0 / 3.d0) )
     ddh0 = be1 / (2.d0 * ks * ks * rho) * (xy - qy) / s1
     sc = rho * h0
     v1c = h0 + dh0
     v2c = ddh0
  endif

end subroutine gcxc_gpu



subroutine pbex_gpu (rho, grho, sx, v1x, v2x)
  !---------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE constants, ONLY : pi
  implicit none
  real(dp), intent(in) :: rho, grho
  real(dp), intent(out):: sx, v1x, v2x
  ! locals
  real(DP) :: kf, agrho, s1, s2, ds, dsg, exunif, fx
  real(DP) :: dxunif, dfx, f1, f2, f3, dfx1
  real(DP) :: p,  amu, ab, c, dfxdp, dfxds, upbe, uge, s, ak, aa
  real(DP), parameter :: third = 1._DP / 3._DP, c1 = 0.75_DP / pi , &
       c2 = 3.093667726280136_DP, c5 = 4._DP * third, &
       c6 = c2*2.51984210, c7=5._DP/6._DP, c8=0.8_DP ! (3pi^2)^(1/3)*2^(4/3)
  real(DP) :: k (6), mu(6), ev(6)
  !           pbe        rpbe        pbesol   pbeq2d      optB88  optB86b
  data k / 0.804_DP,   1.2450D0,   0.804_DP , 0.804_DP ,    0.0 ,  0.0 /, &
       mu/ 0.21951_DP, 0.21951_DP, 0.12345679012345679012_DP,             &
                                   0.12345679012345679,     0.22 , 0.1234/, &
       ev / 1.647127_DP, 0.980118_DP, 0.017399_DP, 1.523671_DP, 0.367229_DP, &
                                   0.011282_DP /  ! a and b parameters of Engel and Vosko
  !
  agrho = sqrt (grho)
  kf = c2 * rho**third
  dsg = 0.5_DP / kf
  s1 = agrho * dsg / rho
  s2 = s1 * s1
  ds = - c5 * s1
  f1 = s2 * mu(1) / k (1)
  f2 = 1._DP + f1
  f3 = k (1) / f2
  fx = k (1) - f3
  exunif = - c1 * kf
  sx = exunif * fx
  dxunif = exunif * third
  dfx1 = f2 * f2
  dfx = 2._DP * mu(1) * s1 / dfx1
  v1x = sx + dxunif * fx + exunif * dfx * ds
  v2x = exunif * dfx * dsg / agrho
  sx = sx * rho
  return
end subroutine pbex_gpu

subroutine pbec_gpu (rho, grho, sc, v1c, v2c)
  USE kinds, ONLY : DP
  implicit none
  real(DP), intent(in) :: rho, grho
  real(DP), intent(out):: sc, v1c, v2c
  !locals
  real(DP), parameter :: ga = 0.031091d0
  real(DP) :: be (3)
  data be / 0.066725d0, 0.046d0,     0.066725d0/
  real(DP), parameter :: third = 1.d0 / 3.d0, pi34 = 0.6203504908994d0
  real(DP), parameter :: xkf = 1.919158292677513d0, xks = 1.128379167095513d0
  real(DP) :: kf, ks, rs, ec, vc, t, expe, af, bf, y, xy, qy
  real(DP) :: s1, h0, dh0, ddh0, sc2D, v1c2D, v2c2D
  !
  rs = pi34 / rho**third
  call pw_gpu (rs, ec, vc)
  kf = xkf / rs
  ks = xks * sqrt (kf)
  t = sqrt (grho) / (2.d0 * ks * rho)
  expe = exp ( - ec / ga)
  af = be(1) / ga * (1.d0 / (expe-1.d0) )
  bf = expe * (vc - ec)
  y = af * t * t
  xy = (1.d0 + y) / (1.d0 + y + y * y)
  qy = y * y * (2.d0 + y) / (1.d0 + y + y * y) **2
  s1 = 1.d0 + be(1) / ga * t * t * xy
  h0 = ga * log (s1)
  dh0 = be(1) * t * t / s1 * ( - 7.d0 / 3.d0 * xy - qy * (af * bf / &
       be(1)-7.d0 / 3.d0) )
  ddh0 = be(1) / (2.d0 * ks * ks * rho) * (xy - qy) / s1
  sc = rho * h0
  v1c = h0 + dh0
  v2c = ddh0
  return
end subroutine pbec_gpu

subroutine pw_gpu (rs, ec, vc)
  USE kinds, ONLY : DP
  implicit none
  real(dp), intent(in) :: rs
  real(dp), intent(out):: ec, vc
  !locals
  real(DP) :: a, b1, b2, c0, c1, c2, c3, d0, d1
  parameter (a = 0.031091d0, b1 = 7.5957d0, b2 = 3.5876d0, c0 = a, &
       c1 = 0.046644d0, c2 = 0.00664d0, c3 = 0.01043d0, d0 = 0.4335d0, &
       d1 = 1.4408d0)
  real(DP) :: lnrs, rs12, rs32, rs2, om, dom, olog
  real(DP) :: a1 (2), b3 (2), b4 (2)
  data a1 / 0.21370d0, 0.026481d0 /, b3 / 1.6382d0, -0.46647d0 /, &
       b4 / 0.49294d0, 0.13354d0 /
  !
  rs12 = sqrt (rs)
  rs32 = rs * rs12
  rs2 = rs**2
  om = 2.d0 * a * (b1 * rs12 + b2 * rs + b3 (1) * rs32 + b4 (1) * rs2)
  dom = 2.d0 * a * (0.5d0 * b1 * rs12 + b2 * rs + 1.5d0 * b3 (1) * rs32 + 2.d0 * b4 (1) * rs2)
  olog = log (1.d0 + 1.0d0 / om)
  ec = - 2.d0 * a * (1.d0 + a1 (1) * rs) * olog
  vc = - 2.d0 * a * (1.d0 + 2.d0 / 3.d0 * a1 (1) * rs) &
       * olog - 2.d0 / 3.d0 * a * (1.d0 + a1 (1) * rs) * dom / (om * (om + 1.d0) )
  return
end subroutine pw_gpu


SUBROUTINE gradcorr_gpu( rho, rho_d, rhog, rhog_d, rho_core, rho_core_d, rhog_core, rhog_core_d, etxc, vtxc, v, v_d )
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
  REAL(DP),    DEVICE, INTENT(IN)    :: rho_d(dfftp%nnr,nspin), rho_core_d(dfftp%nnr)
  COMPLEX(DP), DEVICE, INTENT(IN)    :: rhog_d(ngm,nspin), rhog_core_d(ngm)
  REAL(DP),    DEVICE, INTENT(INOUT) :: v_d(dfftp%nnr,nspin)
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
              CALL gcxc_gpu( arho, grho2(1), sx, sc, v1x, v2x, v1c, v2c )
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
