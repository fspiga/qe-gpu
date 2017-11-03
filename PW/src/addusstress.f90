! Copyright (C) 2001-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE addusstres (sigmanlc)
  !----------------------------------------------------------------------
  !
  !   This routine computes the part of the atomic force which is due
  !   to the dependence of the Q function on the atomic position.
  !   Adds contribution to input sigmanlc, does not sum contributions
  !   from various processors (sum is performed by calling routine)
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE cell_base,  ONLY : omega, tpiba
  USE fft_base,   ONLY : dfftp
  USE gvect,      ONLY : ngm, nl, nlm, gg, g, eigts1, eigts2, eigts3, mill
  USE lsda_mod,   ONLY : nspin
  USE scf,        ONLY : v, vltot
  USE uspp,       ONLY : becsum, okvan
  USE uspp_param, ONLY : upf, lmaxq, nh, nhm
  USE control_flags, ONLY : gamma_only
  USE fft_interfaces,ONLY : fwfft
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(inout) :: sigmanlc (3, 3)
  ! the nonlocal stress
  INTEGER :: ig, nt, ih, jh, ijh, ipol, jpol, is, na, nij
  ! counters
  COMPLEX(DP), ALLOCATABLE :: aux(:), aux1(:,:), aux2(:,:), vg(:,:), qgm(:,:)
  ! work space (complex)
  COMPLEX(DP)              :: cfac
  REAL(dp)                 :: fac(3,nspin), sus(3,3)
  ! auxiliary variables
  REAL(DP) , ALLOCATABLE :: qmod(:), ylmk0(:,:), dylmk0(:,:), tbecsum(:,:)
  ! work space (real)
  !
  !
  sus(:,:) = 0.d0
  !
  ALLOCATE ( aux1(ngm,3), aux2(ngm,nspin), qmod(ngm) )
  ALLOCATE ( ylmk0(ngm,lmaxq*lmaxq), dylmk0(ngm,lmaxq*lmaxq) )
  !
  CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
!$omp parallel do default(shared) private(ig)
  DO ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  ENDDO
!$omp end parallel do
  !
  ! fourier transform of the total effective potential
  !
  ALLOCATE ( vg(ngm,nspin))
  ALLOCATE ( aux(dfftp%nnr) )
  DO is = 1, nspin
     IF ( nspin == 4 .and. is /= 1 ) THEN
        aux(:) = v%of_r(:,is)
     ELSE
        aux(:) = vltot(:) + v%of_r(:,is)
     ENDIF
     CALL fwfft ('Dense', aux, dfftp)
     DO ig = 1, ngm
        vg (ig, is) = aux (nl (ig) )
     ENDDO
  ENDDO
  DEALLOCATE ( aux )
  !
  ! here we compute the integral Q*V for each atom,
  !       I = sum_G i G_a exp(-iR.G) Q_nm v^*
  ! (no contribution from G=0)
  !
  DO ipol = 1, 3
     CALL dylmr2 (lmaxq * lmaxq, ngm, g, gg, dylmk0, ipol)
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp ) THEN
           nij = nh(nt)*(nh(nt)+1)/2
           ALLOCATE (qgm(ngm,nij), tbecsum(nij,nspin) )
           ijh = 0
           DO ih = 1, nh (nt)
              DO jh = ih, nh (nt)
                 ijh = ijh + 1
                 CALL dqvan2 (ngm, ih, jh, nt, qmod, qgm(1,ijh), ylmk0, &
                      dylmk0, ipol)
              ENDDO
           ENDDO
           !
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                 !
                 tbecsum(:,:) = becsum(1:nij,na,1:nspin)
                 !
                 CALL dgemm( 'N', 'N', 2*ngm, nspin, nij, 1.0_dp, &
                      qgm, 2*ngm, tbecsum, nij, 0.0_dp, aux2, 2*ngm )
                 !
!$omp parallel do default(shared) private(is, ig)
                 DO is = 1, nspin
                    DO ig = 1, ngm
                       aux2(ig,is) = aux2(ig,is) * CONJG(vg (ig, is))
                    END DO
                 END DO
!$omp end parallel do
!$omp parallel do default(shared) private(ig, cfac)
                 DO ig = 1, ngm
                    cfac = CONJG( eigts1 (mill (1,ig), na) * &
                                  eigts2 (mill (2,ig), na) * &
                                  eigts3 (mill (3,ig), na) )
                     aux1 (ig,1) = cfac * g (1, ig)
                     aux1 (ig,2) = cfac * g (2, ig)
                     aux1 (ig,3) = cfac * g (3, ig)
                ENDDO
!$omp end parallel do
                CALL DGEMM('T','N', 3, nspin, 2*ngm, 1.0_dp, aux1, 2*ngm, &
                           aux2, 2*ngm, 0.0_dp, fac, 3 )    
                DO is = 1, nspin
                   DO jpol = 1, ipol
                      sus (ipol, jpol) = sus (ipol, jpol) - omega * &
                          fac (jpol, is)
                   ENDDO
                ENDDO
              ENDIF
           ENDDO
           DEALLOCATE ( tbecsum, qgm )
        ENDIF
     ENDDO

  ENDDO

  IF (gamma_only) THEN
     sigmanlc(:,:) = sigmanlc(:,:) + 2.0_dp*sus(:,:)
  ELSE
     sigmanlc(:,:) = sigmanlc(:,:) + sus(:,:)
  ENDIF
  DEALLOCATE (ylmk0, dylmk0)
  DEALLOCATE (aux1, aux2, vg, qmod)

  RETURN

END SUBROUTINE addusstres

#ifdef USE_CUDA


module compute_sigmanl_m
contains
attributes(global) subroutine sigmanl_kernel (nt,nat,ngm, nij,nspin,ityp,qgm,vg,fac)
  USE kinds,      ONLY : DP
  USE gvect,      ONLY : g=>g_d, mill=>mill_d
  USE gvect,      ONLY : eigts1=>eigts1_d, eigts2=>eigts2_d, eigts3=>eigts3_d
  USE uspp,       ONLY : becsum=>becsum_d !copy in shared memory?
  USE  cudafor
  implicit none
  integer, value:: nt,nat,ngm, nij,nspin
  integer:: ityp(nat)
  complex(DP):: qgm(ngm,nij),vg(ngm,nspin)
  real(DP):: fac(3,nspin,nat)
  !
  integer:: ig,is,tx,j,na,laneID,istat
  complex(DP):: aux1(3),aux2(4),cfac !4=max nspin
  real(DP):: rv1,rv2,fac1(4),fac2(4),fac3(4) !4=max nspin

  ! Each block is computing a different atom
  na=blockIdx%y  
  ! If atom type is different from nt, return
  if (ityp(na) .ne. nt) return

  laneID = iand(threadIdx%x, 31)
  tx=threadIdx%x+(blockIdx%x-1)*blockDim%x

  

  do is=1,nspin
    fac1(is) = 0.d0
    fac2(is) = 0.d0
    fac3(is) = 0.d0
  end do

  !if(ig <= ngm) then
  do ig=tx,ngm,blockDim%x*gridDim%x
   ! 
    do is=1,nspin
     aux2(is)=cmplx(0.d0,0.d0)
     do j=1,nij
      aux2(is)=aux2(is)+qgm(ig,j)*becsum(j,na,is)  
     end do
    end do 
   ! 
   do is=1,nspin
    aux2(is)=aux2(is)* conjg(vg(ig,is))
   end do
   !
   cfac = CONJG( eigts1 (mill (1,ig), na) * &
                 eigts2 (mill (2,ig), na) * &
                 eigts3 (mill (3,ig), na) )
   aux1 (1) = cfac * g (1, ig)
   aux1 (2) = cfac * g (2, ig)
   aux1 (3) = cfac * g (3, ig)
   !
   do is=1,nspin
   fac1(is) = fac1(is) + aux1(1)%re*aux2(is)%re + aux1(1)%im*aux2(is)%im
   fac2(is) = fac2(is) + aux1(2)%re*aux2(is)%re + aux1(2)%im*aux2(is)%im
   fac3(is) = fac3(is) + aux1(3)%re*aux2(is)%re + aux1(3)%im*aux2(is)%im 
   end do

   end do
  !endif

 ! Partial sum within warps using shuffle
   do is=1,nspin
      rv1 = fac1(is)
      rv2 = __shfl_down(rv1,1)
      rv1 = rv1 + rv2
      rv2 = __shfl_down(rv1,2)
      rv1 = rv1 + rv2
      rv2 = __shfl_down(rv1,4)
      rv1 = rv1 + rv2
      rv2 = __shfl_down(rv1,8)
      rv1 = rv1 + rv2
      rv2 = __shfl_down(rv1,16)
      rv1 = rv1 + rv2
      if (laneID == 1) then
        istat = atomicadd(fac(1,is,na), rv1)
      endif
      rv1 = fac2(is)
      rv2 = __shfl_down(rv1,1)
      rv1 = rv1 + rv2
      rv2 = __shfl_down(rv1,2)
      rv1 = rv1 + rv2
      rv2 = __shfl_down(rv1,4)
      rv1 = rv1 + rv2
      rv2 = __shfl_down(rv1,8)
      rv1 = rv1 + rv2
      rv2 = __shfl_down(rv1,16)
      rv1 = rv1 + rv2
      if (laneID == 1) then
        istat = atomicadd(fac(2,is,na), rv1)
      endif
      rv1 = fac3(is)
      rv2 = __shfl_down(rv1,1)
      rv1 = rv1 + rv2
      rv2 = __shfl_down(rv1,2)
      rv1 = rv1 + rv2
      rv2 = __shfl_down(rv1,4)
      rv1 = rv1 + rv2
      rv2 = __shfl_down(rv1,8)
      rv1 = rv1 + rv2
      rv2 = __shfl_down(rv1,16)
      rv1 = rv1 + rv2
      if (laneID == 1) then
        istat = atomicadd(fac(3,is,na), rv1)
      endif
    end do

end subroutine sigmanl_kernel 

subroutine  sigmanl_gpu (nt,nat,ngm, nij,nspin,qgm,vg,fac)
     USE kinds,      ONLY : DP
     USE ions_base,  ONLY : ityp_d,ityp
     use cudafor
     implicit none
     integer::nt,nat,ngm,nij,nspin
     complex(DP),device::qgm(ngm,nij),vg(ngm,nspin)
     real(DP),device::fac(3,nspin,nat)
     type(dim3) :: blocks,grid
     type(cudadeviceprop):: prop(2) ! CUDA C keeps changing the size of the struct 
                                    ! (version 17.10 will work with CUDA 9.0, but
                                    ! previuous ones are going to access out-of-bound  
                                    ! memory ). We overallocate 
     integer:: i,ntot,istat,nsm

     
     ntot=0
     do i = 1, nat
       if (ityp (i) == nt) ntot=ntot+1
     end do
 
     ! Assuming all the GPUs are the same
     istat=cudaGetDeviceProperties(prop(1),0)
     nsm=prop(1)%multiprocessorCount

     ! If the number of atoms of the required type is  small compared to the
     ! number of SMs launch a tall grid, otherwise just 2 blocks per atom.
     ! A better scheme could probably be devised.
     blocks=dim3(256,1,1)
     grid=dim3( max(2, ceiling(real(nsm)/ntot)),nat,1)

     call sigmanl_kernel<<<grid,blocks>>> (nt,nat,ngm, nij,nspin,ityp_d,qgm,vg,fac)


end subroutine  

end module compute_sigmanl_m

SUBROUTINE addusstres_gpu (sigmanlc)
  !----------------------------------------------------------------------
  !
  !   This routine computes the part of the atomic force which is due
  !   to the dependence of the Q function on the atomic position.
  !   Adds contribution to input sigmanlc, does not sum contributions
  !   from various processors (sum is performed by calling routine)
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE cell_base,  ONLY : omega, tpiba
  USE fft_base,   ONLY : dfftp
  USE gvect,      ONLY : ngm, nl=>nl_d, nlm, gg=>gg_d, g=>g_d
  USE gvect,      ONLY : eigts1=>eigts1_d, eigts2=>eigts2_d,eigts3=>eigts3_d, mill=>mill_d
  USE lsda_mod,   ONLY : nspin
  USE scf,        ONLY : v, vltot=>vltot_d
  USE uspp,       ONLY : becsum=>becsum_d, okvan
  USE uspp_param, ONLY : upf, lmaxq, nh, nhm
  USE control_flags, ONLY : gamma_only
  USE fft_interfaces,ONLY : fwfft
  USE ylmr2_gpu
  USE compute_sigmanl_m
  use cublas
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(inout) :: sigmanlc (3, 3)
  ! the nonlocal stress
  INTEGER :: ig, nt, ih, jh, ijh, ipol, jpol, is, na, nij, i, j
  ! counters
  COMPLEX(DP), ALLOCATABLE, device :: aux(:), vg(:,:), qgm(:,:)
  ! work space (complex)
  COMPLEX(DP)              :: cfac
  REAL(dp)                 :: fac(3,nspin), sus(3,3), fac1,fac2,fac3
  REAL(dp),device          :: fac_d(3,nspin)
  COMPLEX(DP) :: temp
  ! auxiliary variables
  REAL(DP) , ALLOCATABLE, device :: qmod(:), ylmk0(:,:), dylmk0(:,:), fact_d(:,:,:)
  REAL(DP) , ALLOCATABLE :: fact_h(:,:,:)
  REAL(DP), POINTER, device :: v_of_r_d(:,:)
  ! work space (real)
  !
  !
  sus(:,:) = 0.d0
  !
  ALLOCATE ( qmod(ngm) )
  ALLOCATE ( ylmk0(ngm,lmaxq*lmaxq), dylmk0(ngm,lmaxq*lmaxq) )
  !
  CALL ylmr2_d(lmaxq * lmaxq, ngm, g, gg, ylmk0)
!$cuf kernel do (1) <<<*,*>>>
  DO ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  ENDDO
  !
  ! fourier transform of the total effective potential
  !
  ALLOCATE ( vg(ngm,nspin))
  ALLOCATE ( aux(dfftp%nnr) )

  v_of_r_d => v%of_r_d
  DO is = 1, nspin
     IF ( nspin == 4 .and. is /= 1 ) THEN
      !$cuf kernel do (1) <<<*,*>>>
        do i = lbound(aux,1), ubound(aux,1) 
          aux(i) = v_of_r_d(i,is)
        end do
     ELSE
      !$cuf kernel do (1) <<<*,*>>>
        do i = lbound(aux,1), ubound(aux,1)
          aux(i) = vltot(i) + v_of_r_d(i,is)
        end do
     ENDIF
     CALL fwfft ('Dense', aux, dfftp)
    !$cuf kernel do (1) <<<*,*>>>
     DO ig = 1, ngm
        vg (ig, is) = aux (nl (ig) )
     ENDDO
  ENDDO
  DEALLOCATE ( aux )
  !
  ! here we compute the integral Q*V for each atom,
  !       I = sum_G i G_a exp(-iR.G) Q_nm v^*
  ! (no contribution from G=0)
  !
  ALLOCATE ( fact_d(3,nspin,nat),fact_h(3,nspin,nat) )
  DO ipol = 1, 3
     CALL dylmr2_gpu (lmaxq * lmaxq, ngm, g, gg, dylmk0, ipol)
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp ) THEN
           nij = nh(nt)*(nh(nt)+1)/2
           ALLOCATE (qgm(ngm,nij))
           ijh = 0
           DO ih = 1, nh (nt)
              DO jh = ih, nh (nt)
                 ijh = ijh + 1
                 CALL dqvan2_gpu (ngm, ih, jh, nt, qmod, qgm(1,ijh), ylmk0, &
                      dylmk0, ipol)
              ENDDO
           ENDDO
           !
 
           fact_d=0.d0
           call sigmanl_gpu (nt,nat,ngm, nij,nspin,qgm,vg,fact_d)
           fact_h=fact_d

            
            DO na = 1, nat
                DO is = 1, nspin
                   DO jpol = 1, ipol
                      sus (ipol, jpol) = sus (ipol, jpol) - omega * &
                          fact_h (jpol, is,na)
                   ENDDO
                ENDDO
           ENDDO
          ENDIF
           DEALLOCATE (  qgm )
     ENDDO

  ENDDO
  DEALLOCATE ( fact_d, fact_h )

  IF (gamma_only) THEN
     sigmanlc(:,:) = sigmanlc(:,:) + 2.0_dp*sus(:,:)
  ELSE
     sigmanlc(:,:) = sigmanlc(:,:) + sus(:,:)
  ENDIF
  DEALLOCATE (ylmk0, dylmk0)
  DEALLOCATE (vg, qmod)

  RETURN

END SUBROUTINE addusstres_gpu



! Previuos version done with CUF kernels.
!----------------------------------------------------------------------
SUBROUTINE addusstres_gpu_old(sigmanlc)
  !----------------------------------------------------------------------
  !
  !   This routine computes the part of the atomic force which is due
  !   to the dependence of the Q function on the atomic position.
  !   Adds contribution to input sigmanlc, does not sum contributions
  !   from various processors (sum is performed by calling routine)
  !
  USE kinds,      ONLY : DP
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp
  USE cell_base,  ONLY : omega, tpiba
  USE fft_base,   ONLY : dfftp
  USE gvect,      ONLY : ngm, nl=>nl_d, nlm, gg=>gg_d, g=>g_d
  USE gvect,      ONLY : eigts1=>eigts1_d, eigts2=>eigts2_d,eigts3=>eigts3_d, mill=>mill_d
  USE lsda_mod,   ONLY : nspin
  USE scf,        ONLY : v, vltot=>vltot_d
  USE uspp,       ONLY : becsum=>becsum_d, okvan
  USE uspp_param, ONLY : upf, lmaxq, nh, nhm
  USE control_flags, ONLY : gamma_only
  USE fft_interfaces,ONLY : fwfft
  USE ylmr2_gpu
  use cublas
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(inout) :: sigmanlc (3, 3)
  ! the nonlocal stress
  INTEGER :: ig, nt, ih, jh, ijh, ipol, jpol, is, na, nij, i, j
  ! counters
  COMPLEX(DP), ALLOCATABLE, device :: aux(:), aux1(:,:), aux2(:,:), vg(:,:), qgm(:,:)
  ! work space (complex)
  COMPLEX(DP)              :: cfac
  REAL(dp)                 :: fac(3,nspin), sus(3,3), fac1,fac2,fac3
  REAL(dp),device          :: fac_d(3,nspin)
  COMPLEX(DP) :: temp
  ! auxiliary variables
  REAL(DP) , ALLOCATABLE, device :: qmod(:), ylmk0(:,:), dylmk0(:,:), tbecsum(:,:)
  REAL(DP), POINTER, device :: v_of_r_d(:,:)
  ! work space (real)
  !
  !
  sus(:,:) = 0.d0
  !
  ALLOCATE ( aux1(ngm,3), aux2(ngm,nspin), qmod(ngm) )
  ALLOCATE ( ylmk0(ngm,lmaxq*lmaxq), dylmk0(ngm,lmaxq*lmaxq) )
  !
  CALL ylmr2_d(lmaxq * lmaxq, ngm, g, gg, ylmk0)
!$cuf kernel do (1) <<<*,*>>>
  DO ig = 1, ngm
     qmod (ig) = sqrt (gg (ig) )
  ENDDO
  !
  ! fourier transform of the total effective potential
  !
  ALLOCATE ( vg(ngm,nspin))
  ALLOCATE ( aux(dfftp%nnr) )

  v_of_r_d => v%of_r_d
  DO is = 1, nspin
     IF ( nspin == 4 .and. is /= 1 ) THEN
      !$cuf kernel do (1) <<<*,*>>>
        do i = lbound(aux,1), ubound(aux,1) 
          aux(i) = v_of_r_d(i,is)
        end do
     ELSE
      !$cuf kernel do (1) <<<*,*>>>
        do i = lbound(aux,1), ubound(aux,1)
          aux(i) = vltot(i) + v_of_r_d(i,is)
        end do
     ENDIF
     CALL fwfft ('Dense', aux, dfftp)
    !$cuf kernel do (1) <<<*,*>>>
     DO ig = 1, ngm
        vg (ig, is) = aux (nl (ig) )
     ENDDO
  ENDDO
  DEALLOCATE ( aux )
  !
  ! here we compute the integral Q*V for each atom,
  !       I = sum_G i G_a exp(-iR.G) Q_nm v^*
  ! (no contribution from G=0)
  !
  DO ipol = 1, 3
     CALL dylmr2_gpu (lmaxq * lmaxq, ngm, g, gg, dylmk0, ipol)
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp ) THEN
           nij = nh(nt)*(nh(nt)+1)/2
           ALLOCATE (qgm(ngm,nij), tbecsum(nij,nspin) )
           ijh = 0
           DO ih = 1, nh (nt)
              DO jh = ih, nh (nt)
                 ijh = ijh + 1
                 CALL dqvan2_gpu (ngm, ih, jh, nt, qmod, qgm(1,ijh), ylmk0, &
                      dylmk0, ipol)
              ENDDO
           ENDDO
           !
           DO na = 1, nat
              IF (ityp (na) == nt) THEN
                 !
                ! tbecsum(:,:) = becsum(1:nij,na,1:nspin)
                 !$cuf kernel do(2) <<<*, *>>>
                 do j=1,nspin
                  do i=1,nij
                   tbecsum(i,j)=becsum(i,na,j)
                  end do
                 end do

              if( nspin .gt. 1 ) then
                 !
                 CALL dgemm( 'N', 'N', 2*ngm, nspin, nij, 1.0_dp, &
                      qgm, 2*ngm, tbecsum, nij, 0.0_dp, aux2, 2*ngm )

              else
#if 0
                 CALL my_dgemv( 2*ngm, nij, qgm, tbecsum, aux2 )
                                !( 'N', 2*ngm, nij, 1.0_dp, qgm, 2*ngm, tbecsum, 1, 0.0_dp, aux2, 1 )
#else
                 !$cuf kernel do(1) <<<*, *>>>
                 do j=1,ngm
                   temp = (0.d0, 0.d0)        
                   do i=1,nij
                     temp = temp + qgm(j,i)*tbecsum(i,1)
                   end do
                   aux2(j,1) = temp
                 end do
#endif
               end if
                 !
                 DO is = 1, nspin
                 !$cuf kernel do(1) <<<*, *>>>
                    DO ig = 1, ngm
                       aux2(ig,is) = aux2(ig,is) * CONJG(vg (ig, is))
                    END DO
                 END DO
                 !$cuf kernel do(1) <<<*, *>>>
                 DO ig = 1, ngm
                    cfac = CONJG( eigts1 (mill (1,ig), na) * &
                                  eigts2 (mill (2,ig), na) * &
                                  eigts3 (mill (3,ig), na) )
                     aux1 (ig,1) = cfac * g (1, ig)
                     aux1 (ig,2) = cfac * g (2, ig)
                     aux1 (ig,3) = cfac * g (3, ig)
                ENDDO

              if( nspin .gt. 1 ) then

                CALL DGEMM('T','N', 3, nspin, 2*ngm, 1.0_dp, aux1, 2*ngm, &
                           aux2, 2*ngm, 0.0_dp, fac_d, 3 )    
                fac=fac_d

              else
                
                !CALL myDdotAB( 2*ngm, aux1(:,1), aux2, fac(1,1) )
                !CALL myDdotAB( 2*ngm, aux1(:,2), aux2, fac(2,1) )
                !CALL myDdotAB( 2*ngm, aux1(:,3), aux2, fac(3,1) )
                fac1=0.d0
                fac2=0.d0
                fac3=0.d0

                 !$cuf kernel do(1) <<<*, *>>>
                 do j=1,ngm
                  fac1 = fac1 + aux1(j,1)%re*aux2(j,1)%re + aux1(j,1)%im*aux2(j,1)%im
                  fac2 = fac2 + aux1(j,2)%re*aux2(j,1)%re + aux1(j,2)%im*aux2(j,1)%im
                  fac3 = fac3 + aux1(j,3)%re*aux2(j,1)%re + aux1(j,3)%im*aux2(j,1)%im
                 end do
                
                fac(1,1) = fac1
                fac(2,1) = fac2
                fac(3,1) = fac3

               end if

                DO is = 1, nspin
                   DO jpol = 1, ipol
                      sus (ipol, jpol) = sus (ipol, jpol) - omega * &
                          fac (jpol, is)
                   ENDDO
                ENDDO
              ENDIF
           ENDDO
           DEALLOCATE ( tbecsum, qgm )
        ENDIF
     ENDDO

  ENDDO

  IF (gamma_only) THEN
     sigmanlc(:,:) = sigmanlc(:,:) + 2.0_dp*sus(:,:)
  ELSE
     sigmanlc(:,:) = sigmanlc(:,:) + sus(:,:)
  ENDIF
  DEALLOCATE (ylmk0, dylmk0)
  DEALLOCATE (aux1, aux2, vg, qmod)

  RETURN

END SUBROUTINE addusstres_gpu_old
#endif
