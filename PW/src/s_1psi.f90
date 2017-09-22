!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#ifdef USE_GPU
#define MY_ROUTINE(x)  x##_gpu
#else
#define MY_ROUTINE(x)  x##_cpu
#endif
!----------------------------------------------------------------------------
SUBROUTINE MY_ROUTINE(s_1psi)( npwx, n, psi, spsi )
  !----------------------------------------------------------------------------
  !
  ! ... spsi = S*psi for one wavefunction
  ! ... Wrapper routine - calls calbec and s_psi
  !
  USE kinds,  ONLY : DP
  ! USE uspp,   ONLY : vkb, nkb
  #ifdef USE_GPU
    USE uspp,     ONLY : vkb=>vkb_d
  #else
    USE uspp,     ONLY : vkb, nkb
  #endif
  USE becmod, ONLY : bec_type, becp, calbec
  USE control_flags,    ONLY : gamma_only 
  USE noncollin_module, ONLY : noncolin, npol 
  USE realus, ONLY : real_space, &
              invfft_orbital_gamma, fwfft_orbital_gamma, calbec_rs_gamma, s_psir_gamma, &
              invfft_orbital_k, fwfft_orbital_k, calbec_rs_k, s_psir_k
  USE wvfct,  ONLY: nbnd
  USE cpu_gpu_interface, ONLY : s_psi
  !
  IMPLICIT NONE
  !
  INTEGER     :: npwx, n, ibnd
  COMPLEX(DP) :: psi(npwx*npol,1), spsi(npwx*npol,1)
  #ifdef USE_GPU
    ATTRIBUTES( DEVICE ) :: psi, spsi
  #endif
  !
  !
  CALL start_clock( 's_1psi' )
  !
  IF ( real_space) then
#ifndef USE_GPU
     IF ( gamma_only ) then
        do ibnd=1,nbnd,2
           ! transform the orbital to real space
           call invfft_orbital_gamma(psi,ibnd,nbnd) 
           ! global becp%r is updated
           call calbec_rs_gamma(ibnd,nbnd,becp%r) 
        enddo
        call s_psir_gamma(1,1)
        call fwfft_orbital_gamma(spsi,1,1)
        !
     ELSE
        do ibnd=1,nbnd
           ! transform the orbital to real space
           call invfft_orbital_k(psi,ibnd,nbnd) 
           ! global becp%r is updated
           call calbec_rs_k(ibnd,nbnd) 
        enddo
        call s_psir_k(1,1)
        call fwfft_orbital_k(spsi,1,1)
     END IF
#else
    print *,"REAL_SPACE not implemented!!"; call flush(6); STOP
#endif
  ELSE
     !
     CALL calbec( n, vkb, psi, becp )
     CALL s_psi( npwx, n, 1, psi, spsi )
     !
  END IF
  !
  CALL stop_clock( 's_1psi' )
  !
  RETURN
  !
END SUBROUTINE MY_ROUTINE(s_1psi)
