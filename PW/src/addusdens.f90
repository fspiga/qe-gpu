!
! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
SUBROUTINE addusdens(rho)
  !----------------------------------------------------------------------
  !
  USE realus,               ONLY : addusdens_r
  USE control_flags,        ONLY : tqr
  USE noncollin_module,     ONLY : nspin_mag
  USE fft_base,             ONLY : dfftp
  USE kinds,                ONLY : DP
  !
  IMPLICIT NONE
  !
  !
  REAL(kind=dp), INTENT(inout) :: rho(dfftp%nnr,nspin_mag)
  !
  IF ( tqr ) THEN
     CALL addusdens_r(rho,.true.)
  ELSE
     CALL addusdens_g(rho)
  ENDIF
  !
  RETURN
  !
END SUBROUTINE addusdens
!
!----------------------------------------------------------------------
SUBROUTINE addusdens_g(rho)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the charge density the part which is due to
  !  the US augmentation.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft
  USE gvect,                ONLY : ngm, nl, nlm, gg, g, &
                                   eigts1, eigts2, eigts3, mill
  USE noncollin_module,     ONLY : noncolin, nspin_mag
  USE uspp,                 ONLY : becsum, okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  !
  IMPLICIT NONE
  !
  REAL(kind=dp), INTENT(inout) :: rho(dfftp%nnr,nspin_mag)
  !
  !     here the local variables
  !

  INTEGER :: ig, na, nt, ih, jh, ijh, is, nab, nb, nij
  ! counters

  REAL(DP), ALLOCATABLE :: tbecsum(:,:,:)
  ! \sum_kv <\psi_kv|\beta_l><beta_m|\psi_kv> for each species of atoms
  REAL(DP), ALLOCATABLE :: qmod (:), ylmk0 (:,:)
  ! modulus of G, spherical harmonics
  COMPLEX(DP), ALLOCATABLE :: skk(:,:), aux2(:,:)
  ! structure factors, US contribution to rho
  COMPLEX(DP), ALLOCATABLE ::  aux (:,:), qgm(:)
  ! work space for rho(G,nspin), Fourier transform of q

  IF (.not.okvan) RETURN

  CALL start_clock ('addusdens')

  ALLOCATE (aux ( ngm, nspin_mag) )
  ALLOCATE (qmod( ngm), qgm( ngm) )
  ALLOCATE (ylmk0( ngm, lmaxq * lmaxq) )

  aux (:,:) = (0.d0, 0.d0)
  CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  DO ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  ENDDO
  !
  DO nt = 1, ntyp
     IF ( upf(nt)%tvanp ) THEN
        !
        ! nij = max number of (ih,jh) pairs per atom type nt
        !
        nij = nh(nt)*(nh(nt)+1)/2
        !
        ! count max number of atoms of type nt
        !
        nab = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) nab = nab + 1
        ENDDO
        !
        ALLOCATE ( skk(ngm,nab), tbecsum(nij,nab,nspin_mag), aux2(ngm,nij) )
        !
        nb = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) THEN
              nb = nb + 1
              tbecsum(:,nb,:) = becsum(1:nij,na,1:nspin_mag)
!$omp parallel do default(shared) private(ig)
              DO ig = 1, ngm
                 skk(ig,nb) = eigts1 (mill (1,ig), na) * &
                              eigts2 (mill (2,ig), na) * &
                              eigts3 (mill (3,ig), na)
              ENDDO
!$omp end parallel do
           ENDIF
        ENDDO

        DO is = 1, nspin_mag
           ! sum over atoms
           CALL dgemm( 'N', 'T', 2*ngm, nij, nab, 1.0_dp, skk, 2*ngm,&
                tbecsum(1,1,is), nij, 0.0_dp, aux2, 2*ngm )
           ! sum over lm indices of Q_{lm}
           ijh = 0
           DO ih = 1, nh (nt)
              DO jh = ih, nh (nt)
                 ijh = ijh + 1
                 CALL qvan2 (ngm, ih, jh, nt, qmod, qgm, ylmk0)
!$omp parallel do default(shared) private(ig)
                 DO ig = 1, ngm
                    aux(ig,is) = aux(ig,is) + aux2(ig,ijh)*qgm(ig)
                 ENDDO
!$omp end parallel do
             ENDDO
           ENDDO
        ENDDO
        DEALLOCATE (aux2, tbecsum, skk )
     ENDIF
  ENDDO
  !
  DEALLOCATE (ylmk0)
  DEALLOCATE (qgm, qmod)
  !
  !     convert aux to real space and add to the charge density
  !
#ifdef DEBUG_ADDUSDENS
  CALL start_clock ('addus:fft')
#endif
  DO is = 1, nspin_mag
     psic(:) = (0.d0, 0.d0)
     psic( nl(:) ) = aux(:,is)
     IF (gamma_only) psic( nlm(:) ) = CONJG (aux(:,is))
     CALL invfft ('Dense', psic, dfftp)
     rho(:, is) = rho(:, is) +  DBLE (psic (:) )
  ENDDO
#ifdef DEBUG_ADDUSDENS
  CALL stop_clock ('addus:fft')
#endif
  DEALLOCATE (aux)

  CALL stop_clock ('addusdens')
  RETURN
END SUBROUTINE addusdens_g

#ifdef USE_CUDA
module ylmr2_gpu
contains
attributes(global) subroutine ylmr2_gpu_kernel (lmax,lmax2, ng, g, gg, ylm)
  implicit none
  !
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  REAL(DP), PARAMETER :: pi     = 3.14159265358979323846_DP
  REAL(DP), PARAMETER :: fpi    = 4.0_DP * pi
  integer, intent(in) :: lmax2, ng
  real(DP), intent(in) :: g (3, ng), gg (ng)
  real(DP), intent(out) :: ylm (ng,lmax2)
  !
  ! local variables
  !
  real(DP), parameter :: eps = 1.0d-9
  real(DP) ::  Q(0:4,0:4)  !Allocate Q for the maximum supported size 

  real(DP) :: cost , sent, phi
  real(DP) :: c, gmod
  integer :: lmax, ig, l, m, lm

  attributes(value)::lmax,lmax2,ng
  attributes(device):: g,gg,Q,ylm

  ig= threadIdx%x+BlockDim%x*(BlockIdx%x-1)

  if (ig <= ng) then

  !
  if (lmax == 0) then
     ylm(ig,1) =  sqrt (1.d0 / fpi)
     return
  end if
  !
  !  theta and phi are polar angles, cost = cos(theta)
  !

     gmod = sqrt (gg (ig) )
     if (gmod < eps) then
        cost = 0.d0
     else
        cost = g(3,ig)/gmod
     endif
     !
     !  beware the arc tan, it is defined modulo pi
     !
     if (g(1,ig) > eps) then
        phi  = atan( g(2,ig)/g(1,ig) )
     else if (g(1,ig) < -eps) then
        phi  = atan( g(2,ig)/g(1,ig) ) + pi
     else
        phi  = sign( pi/2.d0,g(2,ig) )
     end if
     sent = sqrt(max(0d0,1.d0-cost*cost))
  !
  !  Q(:,l,m) are defined as sqrt ((l-m)!/(l+m)!) * P(:,l,m) where
  !  P(:,l,m) are the Legendre Polynomials (0 <= m <= l)
  !
    Q(0,0) = 1.d0
    Q(1,0) = cost
    Q(1,1) =-sent/sqrt(2.d0)
    c = sqrt (3.d0 / fpi)
    ylm(ig, 1) = sqrt (1.d0 / fpi)* Q(0,0)
    ylm(ig, 2) = c* Q(1,0)
    ylm(ig, 3) = c*sqrt (2.d0)* Q(1,1) * cos (phi)
    ylm(ig, 4) = c*sqrt (2.d0)* Q(1,1) * sin (phi)
    lm = 4
    do l = 2, lmax
       c = sqrt (DBLE(2*l+1) / fpi)
        !
        !  recursion on l for Q(:,l,m)
        !
        do m = 0, l - 2
           Q(l,m) = cost*(2*l-1)/sqrt(DBLE(l*l-m*m))*Q(l-1,m) &
                       - sqrt(DBLE((l-1)*(l-1)-m*m))/sqrt(DBLE(l*l-m*m))*Q(l-2,m)
        end do
        Q(l,l-1) = cost * sqrt(DBLE(2*l-1)) * Q(l-1,l-1)
        Q(l,l)   = - sqrt(DBLE(2*l-1))/sqrt(DBLE(2*l))*sent*Q(l-1,l-1)

       !
       ! Y_lm, m = 0
       !
       lm = lm + 1
       ylm(ig, lm) = c * Q(l,0)
       !
       do m = 1, l
        !
        ! Y_lm, m > 0
        !
          ylm(ig, lm+2*m-1) = c * sqrt(2.d0) * Q(l,m) * cos (m*phi)
        !
        ! Y_lm, m < 0
        !
          ylm(ig, lm+2*m  ) = c * sqrt(2.d0) * Q(l,m) * sin (m*phi)
       end do
       lm=lm+2*l
    end do
  end if
  return
  end subroutine ylmr2_gpu_kernel

  subroutine ylmr2_d(lmax2, ng, g, gg, ylm)
  !-----------------------------------------------------------------------
  !
  !     Real spherical harmonics ylm(G) up to l=lmax
  !     lmax2 = (lmax+1)^2 is the total number of spherical harmonics
  !     Numerical recursive algorithm based on the one given in Numerical 
  !     Recipes but avoiding the calculation of factorials that generate 
  !     overflow for lmax > 11
  !
  use cudafor
  implicit none
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  REAL(DP), PARAMETER :: pi     = 3.14159265358979323846_DP
  REAL(DP), PARAMETER :: fpi    = 4.0_DP * pi
  integer, intent(in) :: lmax2, ng
  real(DP), intent(in) :: g (3, ng), gg (ng)
  real(DP), intent(out) :: ylm (ng,lmax2)
  attributes(device):: g,gg,ylm
  integer:: lmax
  type(dim3):: grid,tBlock

  !
  ! BEWARE: gg = g(1)^2 + g(2)^2 +g(3)^2  is not checked on input
  !         incorrect results will ensue if the above does not hold
  !

  if (ng < 1 .or. lmax2 < 1) return
  do lmax = 0, 25
     if ((lmax+1)**2 == lmax2) go to 10
  end do
  !call errore (' ylmr', 'l > 25 or wrong number of Ylm required',lmax2)
  print *,'l > 25 or wrong number of Ylm required'
10 continue

  tBlock = dim3(256,1,1)
  grid = dim3(ceiling(real(ng)/tBlock%x),1,1)
  call ylmr2_gpu_kernel<<<grid,tBlock>>>(lmax,lmax2, ng, g, gg, ylm)


  end subroutine ylmr2_d

subroutine dylmr2_gpu (nylm, ngy, g, gg, dylm, ipol)
  !-----------------------------------------------------------------------
  !
  !     compute \partial Y_lm(G) \over \partial (G)_ipol
  !     using simple numerical derivation (SdG)
  !     The spherical harmonics are calculated in ylmr2
  !
  USE kinds, ONLY : DP
  implicit none
  !
  !    here the I/O variables
  !
  integer :: nylm, ngy, ipol
  ! input: number of spherical harmonics
  ! input: the number of g vectors to compute
  ! input: desired polarization
  real(DP),device :: g (3, ngy), gg (ngy), dylm (ngy, nylm)
  ! input: the coordinates of g vectors
  ! input: the moduli of g vectors
  ! output: the spherical harmonics derivatives
  !
  !    and here the local variables
  !
  integer :: ig, lm
  ! counter on g vectors
  ! counter on l,m component

  real(DP), parameter :: delta = 1.d-6
  real(DP), allocatable, device :: dg (:), dgi (:), gx (:,:), ggx (:), ylmaux (:,:)
  ! dg is the finite increment for numerical derivation:
  ! dg = delta |G| = delta * sqrt(gg)
  ! dgi= 1 /(delta * sqrt(gg))
  ! gx = g +/- dg
  ! ggx = gx^2
  !
  allocate ( gx(3,ngy), ggx(ngy), dg(ngy), dgi(ngy), ylmaux(ngy,nylm) )

  !$cuf kernel do(1) <<<*,*>>>
  do ig = 1, ngy
     dg (ig) = delta * sqrt (gg (ig) )
     if (gg (ig) .gt. 1.d-9) then
        dgi (ig) = 1.d0 / dg (ig)
     else
        dgi (ig) = 0.d0
     endif
     gx (1, ig) = g (1, ig)
     gx (2, ig) = g (2, ig)
     gx (3, ig) = g (3, ig)
     gx (ipol, ig) = g (ipol, ig) + dg (ig)
     ggx (ig) = gx (1, ig) * gx (1, ig) + &
                gx (2, ig) * gx (2, ig) + &
                gx (3, ig) * gx (3, ig)
  enddo

  call ylmr2_d (nylm, ngy, gx, ggx, dylm)

  !$cuf kernel do(1) <<<*,*>>>
  do ig = 1, ngy
     gx (ipol, ig) = g (ipol, ig) - dg (ig)
     ggx (ig) = gx (1, ig) * gx (1, ig) + &
                gx (2, ig) * gx (2, ig) + &
                gx (3, ig) * gx (3, ig)
  enddo

  call ylmr2_d (nylm, ngy, gx, ggx, ylmaux)


  do lm = 1, nylm
    !$cuf kernel do(1) <<<*,*>>>
     do ig = 1, ngy
        dylm (ig, lm) = (dylm (ig, lm) - ylmaux (ig,lm)) * 0.5d0 * dgi (ig)
     enddo
  enddo

  deallocate ( gx, ggx, dg, dgi, ylmaux )

  return
end subroutine dylmr2_gpu
end module ylmr2_gpu

module qvan2_gpu_m
contains

attributes(global) subroutine qvan2_kernel(ngy, ih, jh, np, qmod, qg, ylmk0, lmaxq, nbetam, nlx)
  !-----------------------------------------------------------------------
  !
  !    This routine computes the fourier transform of the Q functions
  !    The interpolation table for the radial fourier trasform is stored 
  !    in qrad.
  !
  !    The formula implemented here is
  !
  !     q(g,i,j) = sum_lm (-i)^l ap(lm,i,j) yr_lm(g^) qrad(g,l,i,j)
  !
  !
  USE kinds, ONLY: DP
  USE us, ONLY: qrad=>qrad_d
  !USE uspp_param, ONLY: lmaxq=>lmaxq_d, nbetam=>nbetam_d
  USE uspp, ONLY: lpl=>lpl_d, lpx=>lpx_d, ap=>ap_d, indv=>indv_d, nhtolm=>nhtolm_d
  implicit none
  !
  ! Input variables
  !
  REAL(DP), PARAMETER:: dq = 0.01D0  
  integer,intent(IN) :: ngy, ih, jh, np, lmaxq, nbetam, nlx
  attributes(value):: ngy, ih, jh, np, lmaxq, nbetam, nlx
  ! ngy   :   number of G vectors to compute
  ! ih, jh:   first and second index of Q
  ! np    :   index of pseudopotentials
  !
  real(DP),intent(IN) :: ylmk0 (ngy, lmaxq * lmaxq), qmod (ngy)
  ! ylmk0 :  spherical harmonics
  ! qmod  :  moduli of the q+g vectors
  !
  ! output: the fourier transform of interest
  !
  real(DP),intent(OUT) :: qg (2,ngy)
  attributes(device):: ylmk0,qmod,qg
  !
  !     here the local variables
  !
  real (DP) :: sig
  ! the nonzero real or imaginary part of (-i)^L 
  real (DP), parameter :: sixth = 1.d0 / 6.d0
  !
  integer :: nb, mb, ijv, ivl, jvl, ig, lp, l, lm, i0, i1, i2, i3, ind
  real(DP) :: dqi, qm, px, ux, vx, wx, uvx, pwx, work, qm1
  real(DP) :: uvwx,pwvx,pwux,pxuvx
  ig= threadIdx%x+BlockDim%x*(BlockIdx%x-1)
  if (ig <= ngy) then
  !     compute the indices which correspond to ih,jh
  dqi = 1.0_DP / dq
  nb = indv (ih, np)
  mb = indv (jh, np)
  if (nb.ge.mb) then
     ijv = nb * (nb - 1) / 2 + mb
  else
     ijv = mb * (mb - 1) / 2 + nb
  endif
  ivl = nhtolm(ih, np)
  jvl = nhtolm(jh, np)
  if (nb > nbetam .OR. mb > nbetam) &
       print *,"errorei 1"
       !call errore (' qvan2 ', ' wrong dimensions (1)', MAX(nb,mb))
  if (ivl > nlx .OR. jvl > nlx) &
       print *,"errore 2"
       !call errore (' qvan2 ', ' wrong dimensions (2)', MAX(ivl,jvl))
  qg(1,ig) = 0.d0
  qg(2,ig) = 0.d0

           qm = qmod (ig) * dqi
           px = qm - int (qm)
           ux = 1.d0 - px
           vx = 2.d0 - px
           wx = 3.d0 - px
           i0 = INT( qm ) + 1
           i1 = i0 + 1
           i2 = i0 + 2
           i3 = i0 + 3
           uvx = ux * vx * sixth
           pwx = px * wx * 0.5d0

  do lm = 1, lpx (ivl, jvl)
     lp = lpl (ivl, jvl, lm)
      if (lp == 1) then
         l = 1
         sig = 1.0d0
         ind = 1
      elseif ( lp <= 4) then
         l = 2
         sig =-1.0d0
         ind = 2
      elseif ( lp <= 9 ) then
         l = 3
         sig =-1.0d0
         ind = 1
      elseif ( lp <= 16 ) then
         l = 4
          sig = 1.0d0
         ind = 2
      elseif ( lp <= 25 ) then
         l = 5
         sig = 1.0d0
         ind = 1
      elseif ( lp <= 36 ) then
         l = 6
         sig =-1.0d0
         ind = 2
      else
         l = 7
         sig =-1.0d0
         ind = 1
      endif
      sig = sig * ap (lp, ivl, jvl)
           work = qrad (i0, ijv, l, np) * uvx * wx + &
                  qrad (i1, ijv, l, np) * pwx * vx - &
                  qrad (i2, ijv, l, np) * pwx * ux + &
                  qrad (i3, ijv, l, np) * px * uvx
        qg (ind,ig) = qg (ind,ig) + sig * ylmk0 (ig, lp) * work
  end do

  end if

  return
end subroutine qvan2_kernel

subroutine qvan2_gpu(ngy, ih, jh, np, qmod, qg, ylmk0)
  use cudafor
  use kinds, ONLY: DP
  use uspp_param, ONLY: lmaxq, nbetam
  use uspp, ONLY: nlx
  implicit none
  integer,intent(IN) :: ngy, ih, jh, np
  !
  real(DP),intent(IN) :: ylmk0 (ngy, lmaxq * lmaxq), qmod (ngy)
  !
  !complex(DP),intent(OUT) :: qg (ngy)
  !pgi$ ignore_tkr (tkr) qg
  real(DP),device,intent(OUT) :: qg (2,ngy)
  attributes(device):: ylmk0,qmod

  type(dim3):: grid,tBlock

  tBlock = dim3(256,1,1)
  grid = dim3(ceiling(real(ngy)/tBlock%x),1,1)
  call qvan2_kernel<<<grid,tBlock>>>(ngy, ih, jh, np, qmod, qg, ylmk0,lmaxq,nbetam,nlx)
end subroutine qvan2_gpu

end module qvan2_gpu_m

!----------------------------------------------------------------------
SUBROUTINE addusdens_gpu(rho, rho_h)
  !----------------------------------------------------------------------
  !
  USE realus,               ONLY : addusdens_r
  USE control_flags,        ONLY : tqr
  USE noncollin_module,     ONLY : nspin_mag
  USE fft_base,             ONLY : dfftp
  USE kinds,                ONLY : DP
  !
  IMPLICIT NONE
  !
  !
  REAL(kind=dp), DEVICE, INTENT(inout) :: rho(dfftp%nnr,nspin_mag)
  REAL(kind=dp), INTENT(inout) :: rho_h(dfftp%nnr,nspin_mag)
  !
  IF ( tqr ) THEN
     ! Falling back to CPU path here. Can port addusdens_r to GPU if needed for performance.
     rho_h = rho
     CALL addusdens_r(rho_h,.true.)
     rho = rho_h
  ELSE
     CALL addusdens_g_gpu(rho)
     rho_h = rho
  ENDIF
  !
  RETURN
  !
END SUBROUTINE addusdens_gpu
!
!----------------------------------------------------------------------
SUBROUTINE addusdens_g_gpu(rho)
  !----------------------------------------------------------------------
  !
  !  This routine adds to the charge density the part which is due to
  !  the US augmentation.
  !
  USE kinds,                ONLY : DP
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft
  USE gvect,                ONLY : ngm, nl=>nl_d, nlm, gg=>gg_d, g=>g_d, &
                                   eigts1=>eigts1_d, eigts2=>eigts2_d, eigts3=>eigts3_d, mill=>mill_d
  USE noncollin_module,     ONLY : noncolin, nspin_mag
  USE uspp,                 ONLY : becsum=>becsum_d, okvan
  USE uspp_param,           ONLY : upf, lmaxq, nh
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic=>psic_d
  !
  USE ylmr2_gpu
  USE qvan2_gpu_m
  !
  USE cublas
  !
  IMPLICIT NONE
  !
  REAL(kind=dp), DEVICE, INTENT(inout) :: rho(dfftp%nnr,nspin_mag)
  !
  !     here the local variables
  !

  INTEGER :: ig, na, nt, ih, jh, ijh, is, nab, nb, nij, ii,jj,kk
  ! counters

  REAL(DP), DEVICE, ALLOCATABLE :: tbecsum(:,:,:)
  ! \sum_kv <\psi_kv|\beta_l><beta_m|\psi_kv> for each species of atoms
  REAL(DP), DEVICE, ALLOCATABLE :: qmod (:), ylmk0 (:,:)
  ! modulus of G, spherical harmonics
  COMPLEX(DP), DEVICE, ALLOCATABLE :: skk(:,:), aux2(:,:)
  ! structure factors, US contribution to rho
  COMPLEX(DP), DEVICE, ALLOCATABLE ::  aux (:,:), qgm(:)
  ! work space for rho(G,nspin), Fourier transform of q

  IF (.not.okvan) RETURN

  CALL start_clock ('addusdens')

  ALLOCATE (aux ( ngm, nspin_mag) )
  ALLOCATE (qmod( ngm), qgm( ngm) )
  ALLOCATE (ylmk0( ngm, lmaxq * lmaxq) )

  aux (:,:) = (0.d0, 0.d0)

  CALL ylmr2_d (lmaxq * lmaxq, ngm, g, gg, ylmk0)

!$cuf kernel do(1) <<<*,*>>>
  DO ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  ENDDO

  !
  DO nt = 1, ntyp
     IF ( upf(nt)%tvanp ) THEN
        !
        ! nij = max number of (ih,jh) pairs per atom type nt
        !
        nij = nh(nt)*(nh(nt)+1)/2
        !
        ! count max number of atoms of type nt
        !
        nab = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) nab = nab + 1
        ENDDO
        !
        ALLOCATE ( skk(ngm,nab), tbecsum(nij,nab,nspin_mag), aux2(ngm,nij) )
        !
        nb = 0
        DO na = 1, nat
           IF ( ityp(na) == nt ) THEN
              nb = nb + 1

!$cuf kernel do(2) <<<*,*>>>
              DO kk = 1, nspin_mag
                 DO ii = 1, nij
                   tbecsum(ii,nb,kk) = becsum(ii,na,kk)
                 ENDDO
              ENDDO

!$cuf kernel do(1) <<<*,*>>>
              DO ig = 1, ngm
                 skk(ig,nb) = eigts1 (mill (1,ig), na) * &
                                eigts2 (mill (2,ig), na) * &
                                eigts3 (mill (3,ig), na)
              ENDDO

           ENDIF
        ENDDO

        DO is = 1, nspin_mag
           ! sum over atoms
           CALL cublasDgemm( 'N', 'T', 2*ngm, nij, nab, 1.0_dp, skk, 2*ngm, &
                tbecsum(1,1,is), nij, 0.0_dp, aux2, 2*ngm )

           ! sum over lm indices of Q_{lm}
           ijh = 0
           DO ih = 1, nh (nt)
              DO jh = ih, nh (nt)
                 ijh = ijh + 1

                 CALL qvan2_gpu (ngm, ih, jh, nt, qmod, qgm, ylmk0)
                 
!$cuf kernel do(1) <<<*,*>>>
                 DO ig = 1, ngm
                    aux(ig,is) = aux(ig,is) + aux2(ig,ijh)*qgm(ig)
                 ENDDO

             ENDDO
           ENDDO
        ENDDO
        DEALLOCATE (aux2, tbecsum, skk )
     ENDIF
  ENDDO
  !
  DEALLOCATE (ylmk0)
  DEALLOCATE (qgm, qmod)
  !
  !     convert aux to real space and add to the charge density
  !
#ifdef DEBUG_ADDUSDENS
  CALL start_clock ('addus:fft')
#endif

  DO is = 1, nspin_mag

     psic(:) = (0.d0, 0.d0)

!$cuf kernel do(1) <<<*,*>>>
     DO ii = lbound(nl,1),ubound(nl,1)
       psic( nl(ii) ) = aux(ii,is)
     ENDDO

     IF (gamma_only) print *,"GAMMA ONLY NOT SUPPORTED YET!!!" !psic( nlm(:) ) = CONJG (aux(:,is))

     CALL invfft ('Dense', psic, dfftp)

!$cuf kernel do(1) <<<*,*>>>
     DO ii = lbound(rho,1),ubound(rho,1)
       rho(ii, is) = rho(ii, is) +  DBLE (psic (ii) )
     ENDDO

  ENDDO
#ifdef DEBUG_ADDUSDENS
  CALL stop_clock ('addus:fft')
#endif
  DEALLOCATE (aux)

  CALL stop_clock ('addusdens')
  RETURN
END SUBROUTINE addusdens_g_gpu

#endif
