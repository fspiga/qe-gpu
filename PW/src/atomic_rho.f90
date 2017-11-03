!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine atomic_rho (rhoa, nspina)
  !-----------------------------------------------------------------------
  ! This routine calculates rhoa as the superposition of atomic charges.
  !
  ! nspina is the number of spin components to be calculated
  !
  ! if nspina = 1 the total atomic charge density is calculated
  ! if nspina = 2 the spin up and spin down atomic charge densities are
  !               calculated assuming an uniform atomic spin-polarization
  !               equal to starting_magnetization(nt)
  ! if nspina = 4 noncollinear case. The total density is calculated
  !               in the first component and the magnetization vector 
  !               in the other three.
  !
  ! NB: nspina may not be equal to nspin because in some cases (as in update)
  ! the total charge only could be needed, even in a LSDA calculation.
  !
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE atom,                 ONLY : rgrid, msh
  USE ions_base,            ONLY : ntyp => nsp
  USE cell_base,            ONLY : tpiba, omega
  USE gvect,                ONLY : ngm, ngl, gstart, nl, nlm, gl, igtongl
  USE lsda_mod,             ONLY : starting_magnetization, lsda
  USE vlocal,               ONLY : strf
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic
  USE noncollin_module,     ONLY : angle1, angle2
  USE uspp_param,           ONLY : upf
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft

  !
  implicit none
  !
  integer :: nspina
  ! the number of spin polarizations
  real(DP) :: rhoa (dfftp%nnr, nspina)
  ! the output atomic charge
  !
  ! local variables
  !
  real(DP) :: rhoneg, rhoima, gx
  real(DP), allocatable :: rhocgnt (:), aux (:)
  complex(DP), allocatable :: rhocg (:,:)
  integer :: ir, is, ig, igl, nt, ndm
  !
  ! superposition of atomic charges contained in the array rho_at
  ! (read from pseudopotential files)
  !
  ! allocate work space (psic must already be allocated)
  !
  CALL start_clock( 'atomic_rho' )
  allocate (rhocg(  ngm, nspina))    
  ndm = MAXVAL ( msh(1:ntyp) )
  allocate (aux(ndm))    
  allocate (rhocgnt( ngl))    
  rhoa(:,:) = 0.d0
  rhocg(:,:) = (0.d0,0.d0)

  do nt = 1, ntyp
     !
     ! Here we compute the G=0 term
     !
     if (gstart == 2) then
        do ir = 1, msh (nt)
           aux (ir) = upf(nt)%rho_at (ir)
        enddo
        call simpson (msh (nt), aux, rgrid(nt)%rab, rhocgnt (1) )
     endif
     !
     ! Here we compute the G<>0 term
     !
     do igl = gstart, ngl
        gx = sqrt (gl (igl) ) * tpiba
        do ir = 1, msh (nt)
           if (rgrid(nt)%r(ir) < 1.0d-8) then
              aux(ir) = upf(nt)%rho_at(ir)
           else
              aux(ir) = upf(nt)%rho_at(ir) * &
                        sin(gx*rgrid(nt)%r(ir)) / (rgrid(nt)%r(ir)*gx)
           endif
        enddo
        call simpson (msh (nt), aux, rgrid(nt)%rab, rhocgnt (igl) )
     enddo
     !
     ! we compute the 3D atomic charge in reciprocal space
     !
     if (nspina == 1) then
        do ig = 1, ngm
           rhocg(ig,1) = rhocg(ig,1) + &
                         strf(ig,nt) * rhocgnt(igtongl(ig)) / omega
        enddo
     else if (nspina == 2) then
        do ig = 1, ngm
           rhocg(ig,1) = rhocg(ig,1) + &
                         0.5d0 * ( 1.d0 + starting_magnetization(nt) ) * &
                         strf(ig,nt) * rhocgnt(igtongl(ig)) / omega
           rhocg(ig,2) = rhocg(ig,2) + &
                         0.5d0 * ( 1.d0 - starting_magnetization(nt) ) * &
                         strf(ig,nt) * rhocgnt(igtongl(ig)) / omega
        enddo
     else
!
!    Noncolinear case
!
        do ig = 1,ngm
           rhocg(ig,1) = rhocg(ig,1) + &
                strf(ig,nt)*rhocgnt(igtongl(ig))/omega

           ! Now, the rotated value for the magnetization

           rhocg(ig,2) = rhocg(ig,2) + &
                starting_magnetization(nt)* &
                sin(angle1(nt))*cos(angle2(nt))* &
                strf(ig,nt)*rhocgnt(igtongl(ig))/omega
           rhocg(ig,3) = rhocg(ig,3) + &
                starting_magnetization(nt)* &
                sin(angle1(nt))*sin(angle2(nt))* &
                strf(ig,nt)*rhocgnt(igtongl(ig))/omega
           rhocg(ig,4) = rhocg(ig,4) + &
                starting_magnetization(nt)* &
                cos(angle1(nt))* &
                strf(ig,nt)*rhocgnt(igtongl(ig))/omega
        end do
     endif
  enddo

  deallocate (rhocgnt)
  deallocate (aux)

  do is = 1, nspina
     !
     ! and we return to real space
     !
     psic(:) = (0.d0,0.d0)
     psic (nl (:) ) = rhocg (:, is)
     if (gamma_only) psic ( nlm(:) ) = CONJG( rhocg (:, is) )
     CALL invfft ('Dense', psic, dfftp)
     !
     ! we check that everything is correct
     !
     rhoneg = 0.d0
     rhoima = 0.d0
     do ir = 1, dfftp%nnr
        rhoneg = rhoneg + MIN (0.d0,  DBLE (psic (ir)) )
        rhoima = rhoima + abs (AIMAG (psic (ir) ) )
     enddo
     rhoneg = omega * rhoneg / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
     rhoima = omega * rhoima / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
     !
     call mp_sum(  rhoneg, intra_bgrp_comm )
     call mp_sum(  rhoima, intra_bgrp_comm )
     !
     IF ( rhoima > 1.0d-4 ) THEN
        WRITE( stdout,'(5x,"Check: imaginary charge or magnetization=",&
          & f12.6," (component ",i1,") set to zero")') rhoima, is
     END IF
     IF ( (is == 1) .OR. lsda ) THEN
        !
        IF ( (rhoneg < -1.0d-4) ) THEN
           IF ( lsda ) THEN 
              WRITE( stdout,'(5x,"Check: negative starting charge=", &
                   &"(component",i1,"):",f12.6)') is, rhoneg
           ELSE
              WRITE( stdout,'(5x,"Check: negative starting charge=", &
          &          f12.6)') rhoneg
           END IF
        END IF
     END IF
     !
     ! set imaginary terms to zero - negative terms are not set to zero
     ! because it is basically useless to do it in real space: negative
     ! charge will re-appear when Fourier-transformed back and forth
     !
     DO ir = 1, dfftp%nnr
        rhoa (ir, is) =  DBLE (psic (ir))
     END DO
     !
  enddo

  deallocate (rhocg)
  CALL stop_clock( 'atomic_rho' )
  return
end subroutine atomic_rho

#ifdef USE_CUDA
module compute_rhocgnt_gpu_m
contains
  attributes(global) subroutine  compute_rhocgnt_gpu(n, mesh, tpiba, gl, r, rho_at, rab, rhocgnt)
  use cudafor
  use kinds, ONLY: DP
  implicit none
 
  integer, value :: n, mesh
  real(DP), value :: tpiba
  real(DP), device, intent(in) :: gl(n), r(mesh), rho_at(mesh), rab(mesh)
  real(DP), device, intent(out) :: rhocgnt(n)
 
  integer :: tx, ty, igl, ir
  real(DP):: mysum, val, gx, x
 
  tx = threadIdx%x
  ty = threadIdx%y
 
  igl = (blockIdx%x - 1) * blockDim%y + ty
 
  if (igl > n) return
 
  gx = sqrt(gl(igl)) * tpiba
  mysum = 0.d0
 
  do ir = tx, mesh, blockDim%x 
    val = rho_at(ir) * rab(ir)

    if (r(ir)  .ge. 1.0d-8) then
      val = val * sin (gx*r(ir)) / (gx *r(ir)) 
    end if

    if (ir == 1 .or. ir == mesh) then
      mysum = mysum + val
    else if (mod(ir,2)) then
      mysum = mysum + 2.d0*val
    else
      mysum = mysum + 4.d0*val
    endif
  end do 
 
 ! Reduce by warp
       val = __shfl_down(mysum,1)
       mysum = mysum + val
       val = __shfl_down(mysum,2)
       mysum = mysum + val
       val = __shfl_down(mysum,4)
       mysum = mysum + val
       val = __shfl_down(mysum,8)
       mysum = mysum + val
       val = __shfl_down(mysum,16)
       mysum = mysum + val
 
       if (tx == 1) then
         rhocgnt(igl) =  mysum / 3.d0 
       endif 
 
  end subroutine compute_rhocgnt_gpu
end module compute_rhocgnt_gpu_m

subroutine atomic_rho_gpu (rhoa, nspina)
  !-----------------------------------------------------------------------
  ! This routine calculates rhoa as the superposition of atomic charges.
  !
  ! nspina is the number of spin components to be calculated
  !
  ! if nspina = 1 the total atomic charge density is calculated
  ! if nspina = 2 the spin up and spin down atomic charge densities are
  !               calculated assuming an uniform atomic spin-polarization
  !               equal to starting_magnetization(nt)
  ! if nspina = 4 noncollinear case. The total density is calculated
  !               in the first component and the magnetization vector 
  !               in the other three.
  !
  ! NB: nspina may not be equal to nspin because in some cases (as in update)
  ! the total charge only could be needed, even in a LSDA calculation.
  !
  !
  USE kinds,                ONLY : DP
  USE io_global,            ONLY : stdout
  USE atom,                 ONLY : rgrid, msh
  USE ions_base,            ONLY : ntyp => nsp
  USE cell_base,            ONLY : tpiba, omega
  USE gvect,                ONLY : ngm, ngl, gstart, nl=>nl_d, nlm, gl=>gl_d, igtongl=>igtongl_d
  USE lsda_mod,             ONLY : starting_magnetization, lsda
  USE vlocal,               ONLY : strf=>strf_d
  USE control_flags,        ONLY : gamma_only
  USE wavefunctions_module, ONLY : psic=>psic_d
  USE noncollin_module,     ONLY : angle1, angle2
  USE uspp_param,           ONLY : upf
  USE mp_bands,             ONLY : intra_bgrp_comm
  USE mp,                   ONLY : mp_sum
  USE fft_base,             ONLY : dfftp
  USE fft_interfaces,       ONLY : invfft
  USE compute_rhocgnt_gpu_m
  USE cudafor

  !
  implicit none
  !
  integer :: nspina
  ! the number of spin polarizations
  real(DP), device :: rhoa (dfftp%nnr, nspina)
  ! the output atomic charge
  !
  ! local variables
  !
  real(DP) :: rhoneg, rhoima, gx
  real(DP) :: rhocgnt1, sm
  real(DP), allocatable, device :: rhocgnt (:), aux (:)
  real(DP), pointer, device :: rho_at_d(:)
  complex(DP), allocatable, device :: rhocg (:,:)
  integer :: i, ir, is, ig, igl, nt, ndm, blocks
  type(dim3) :: threads
  !
  ! superposition of atomic charges contained in the array rho_at
  ! (read from pseudopotential files)
  !
  ! allocate work space (psic must already be allocated)
  !
  CALL start_clock( 'atomic_rho' )
  allocate (rhocg(  ngm, nspina))    
  ndm = MAXVAL ( msh(1:ntyp) )
  allocate (aux(ndm))    
  allocate (rhocgnt( ngl))    
  rhoa(:,:) = 0.d0
  rhocg(:,:) = (0.d0,0.d0)

  do nt = 1, ntyp
     rho_at_d => upf(nt)%rho_at_d
     !
     ! Here we compute the G=0 term
     !
     if (gstart == 2) then
        !$cuf kernel do (1) <<<*, *>>>
        do ir = 1, msh (nt)
           aux (ir) = rho_at_d (ir)
        enddo
        call simpson_gpu (msh (nt), aux, rgrid(nt)%rab_d, rhocgnt1 )

        rhocgnt(1) = rhocgnt1
     endif
     !
     ! Here we compute the G<>0 term
     !
     threads = dim3(32, 8, 1)
     blocks = ceiling(real(ngl - gstart + 1)/8)
     call compute_rhocgnt_gpu<<<blocks, threads>>>(ngl - gstart + 1, msh(nt), tpiba, gl(gstart), rgrid(nt)%r_d, & 
       rho_at_d, rgrid(nt)%rab_d, rhocgnt(gstart))
     !
     ! we compute the 3D atomic charge in reciprocal space
     !
     if (nspina == 1) then
        !$cuf kernel do (1) <<<*, *>>>
        do ig = 1, ngm
           rhocg(ig,1) = rhocg(ig,1) + &
                         strf(ig,nt) * rhocgnt(igtongl(ig)) / omega
        enddo
     else if (nspina == 2) then
        sm = starting_magnetization(nt)
        !$cuf kernel do (1) <<<*, *>>>
        do ig = 1, ngm
           rhocg(ig,1) = rhocg(ig,1) + &
                         0.5d0 * ( 1.d0 + sm ) * &
                         strf(ig,nt) * rhocgnt(igtongl(ig)) / omega
           rhocg(ig,2) = rhocg(ig,2) + &
                         0.5d0 * ( 1.d0 - sm ) * &
                         strf(ig,nt) * rhocgnt(igtongl(ig)) / omega
        enddo
     else
       print*, "atomic_rho: NONCOLINEAR NOT IMPLEMENTED!"
       flush(6); stop
!
!    Noncolinear case
!
!        do ig = 1,ngm
!           rhocg(ig,1) = rhocg(ig,1) + &
!                strf(ig,nt)*rhocgnt(igtongl(ig))/omega
!
!           ! Now, the rotated value for the magnetization
!
!           rhocg(ig,2) = rhocg(ig,2) + &
!                starting_magnetization(nt)* &
!                sin(angle1(nt))*cos(angle2(nt))* &
!                strf(ig,nt)*rhocgnt(igtongl(ig))/omega
!           rhocg(ig,3) = rhocg(ig,3) + &
!                starting_magnetization(nt)* &
!                sin(angle1(nt))*sin(angle2(nt))* &
!                strf(ig,nt)*rhocgnt(igtongl(ig))/omega
!           rhocg(ig,4) = rhocg(ig,4) + &
!                starting_magnetization(nt)* &
!                cos(angle1(nt))* &
!                strf(ig,nt)*rhocgnt(igtongl(ig))/omega
!        end do
     endif
  enddo

  deallocate (rhocgnt)
  deallocate (aux)

  do is = 1, nspina
     !
     ! and we return to real space
     !
     !$cuf kernel do (1) <<<*, *>>>
     do i = lbound(psic, 1), ubound(psic, 1)
       psic(i) = (0.d0,0.d0)
     end do

     !$cuf kernel do (1) <<<*, *>>>
     do i = lbound(nl, 1), ubound(nl, 1)
       psic (nl (i) ) = rhocg (i, is)
     end do

     if (gamma_only) then
       !psic ( nlm(:) ) = CONJG( rhocg (:, is) )
       print*, "atomic_rho: GAMMA ONLY NOT SUPPORTED!"
       flush(6); stop
     endif

     CALL invfft ('Dense', psic, dfftp)
     !
     ! we check that everything is correct
     !
     rhoneg = 0.d0
     rhoima = 0.d0
     !$cuf kernel do (1) <<<*, *>>>
     do ir = 1, dfftp%nnr
        rhoneg = rhoneg + MIN (0.d0,  DBLE (psic (ir)) )
        rhoima = rhoima + abs (AIMAG (psic (ir) ) )
     enddo
     rhoneg = omega * rhoneg / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
     rhoima = omega * rhoima / (dfftp%nr1 * dfftp%nr2 * dfftp%nr3)
     !
     call mp_sum(  rhoneg, intra_bgrp_comm )
     call mp_sum(  rhoima, intra_bgrp_comm )
     !
     IF ( rhoima > 1.0d-4 ) THEN
        WRITE( stdout,'(5x,"Check: imaginary charge or magnetization=",&
          & f12.6," (component ",i1,") set to zero")') rhoima, is
     END IF
     IF ( (is == 1) .OR. lsda ) THEN
        !
        IF ( (rhoneg < -1.0d-4) ) THEN
           IF ( lsda ) THEN 
              WRITE( stdout,'(5x,"Check: negative starting charge=", &
                   &"(component",i1,"):",f12.6)') is, rhoneg
           ELSE
              WRITE( stdout,'(5x,"Check: negative starting charge=", &
          &          f12.6)') rhoneg
           END IF
        END IF
     END IF
     !
     ! set imaginary terms to zero - negative terms are not set to zero
     ! because it is basically useless to do it in real space: negative
     ! charge will re-appear when Fourier-transformed back and forth
     !
     !$cuf kernel do (1) <<<*, *>>>
     DO ir = 1, dfftp%nnr
        rhoa (ir, is) =  DBLE (psic (ir))
     END DO
     !
  enddo

  deallocate (rhocg)
  CALL stop_clock( 'atomic_rho' )
  return
end subroutine atomic_rho_gpu
#endif
