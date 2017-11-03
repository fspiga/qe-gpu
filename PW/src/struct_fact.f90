!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine struc_fact (nat, tau, ntyp, ityp, ngm, g, bg, nr1, nr2, &
     nr3, strf, eigts1, eigts2, eigts3)
  !----------------------------------------------------------------------
  !
  !   calculate the structure factors for each type of atoms in the unit
  !   cell
  !
  USE kinds
  USE constants, ONLY : tpi
#ifdef USE_CUDA
  USE vlocal, ONLY : strf_d      
  USE gvect,        ONLY : eigts1_d, eigts2_d, eigts3_d
#endif

  implicit none
  !
  !   Here the dummy variables
  !

  integer :: nat, ntyp, ityp (nat), ngm, nr1, nr2, nr3
  ! input: the number of atom in the unit cel
  ! input: the number of atom types
  ! input: for each atom gives the type
  ! input: the number of G vectors
  ! input: fft dimension along x
  ! input: fft dimension along y
  ! input: fft dimension along z

  real(DP) :: bg (3, 3), tau (3, nat), g (3, ngm)
  ! input: reciprocal crystal basis vectors
  ! input: the positions of the atoms in the c
  ! input: the coordinates of the g vectors

  complex(DP) :: strf (ngm, ntyp),        &
                      eigts1 ( -nr1:nr1, nat), &
                      eigts2 ( -nr2:nr2, nat), &
                      eigts3 ( -nr3:nr3, nat)
  ! output: the structure factor
  !
  ! output: the phases e^{-iG\tau_s}
  !
  !
  !    here the local variables
  !
  integer :: nt, na, ng, n1, n2, n3, ipol
  ! counter over atom type
  ! counter over atoms
  ! counter over G vectors
  ! counter over fft dimension along x
  ! counter over fft dimension along y
  ! counter over fft dimension along z
  ! counter over polarizations

  real(DP) :: arg, bgtau (3)
  ! the argument of the exponent
  ! scalar product of bg and tau

  CALL start_clock( 'struct_fact' )

  strf(:,:) = (0.d0,0.d0)
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) .eq.nt) then
           do ng = 1, ngm
              arg = (g (1, ng) * tau (1, na) + g (2, ng) * tau (2, na) &
                   + g (3, ng) * tau (3, na) ) * tpi
              strf (ng, nt) = strf (ng, nt) + CMPLX(cos (arg), -sin (arg),kind=DP)
           enddo
        endif
     enddo
  enddo

  do na = 1, nat
     do ipol = 1, 3
        bgtau (ipol) = bg (1, ipol) * tau (1, na) + &
                       bg (2, ipol) * tau (2, na) + &
                       bg (3, ipol) * tau (3, na)
     enddo
     do n1 = - nr1, nr1
        arg = tpi * n1 * bgtau (1)
        eigts1 (n1, na) = CMPLX(cos (arg), - sin (arg) ,kind=DP)
     enddo
     do n2 = - nr2, nr2
        arg = tpi * n2 * bgtau (2)
        eigts2 (n2, na) = CMPLX(cos (arg), - sin (arg) ,kind=DP)
     enddo
     do n3 = - nr3, nr3
        arg = tpi * n3 * bgtau (3)
        eigts3 (n3, na) = CMPLX(cos (arg), - sin (arg) ,kind=DP)
     enddo
  enddo

#ifdef USE_CUDA
  strf_d   = strf
  eigts1_d = eigts1
  eigts2_d = eigts2
  eigts3_d = eigts3
#endif

  CALL stop_clock( 'struct_fact' )
  return
end subroutine struc_fact

#ifdef USE_CUDA
subroutine struc_fact_gpu (nat, tau, ntyp, ityp, ngm, g, bg, nr1, nr2, &
     nr3, strf, eigts1, eigts2, eigts3)
  !----------------------------------------------------------------------
  !
  !   calculate the structure factors for each type of atoms in the unit
  !   cell
  !
  USE kinds
  USE constants, ONLY : tpi

  implicit none
  !
  !   Here the dummy variables
  !

  integer :: nat, ntyp, ityp (nat), ngm, nr1, nr2, nr3
  ! input: the number of atom in the unit cel
  ! input: the number of atom types
  ! input: for each atom gives the type
  ! input: the number of G vectors
  ! input: fft dimension along x
  ! input: fft dimension along y
  ! input: fft dimension along z

  real(DP) :: bg (3, 3), tau (3, nat), g (3, ngm)
  ! input: reciprocal crystal basis vectors
  ! input: the positions of the atoms in the c
  ! input: the coordinates of the g vectors

  complex(DP) :: strf (ngm, ntyp),        &
                      eigts1 ( -nr1:nr1, nat), &
                      eigts2 ( -nr2:nr2, nat), &
                      eigts3 ( -nr3:nr3, nat)
  ! output: the structure factor
  !
  ! output: the phases e^{-iG\tau_s}
  !
  attributes(device) :: tau, g, strf, eigts1, eigts2, eigts3
  !
  !    here the local variables
  !
  integer :: nt, na, ng, n1, n2, n3, ipol
  ! counter over atom type
  ! counter over atoms
  ! counter over G vectors
  ! counter over fft dimension along x
  ! counter over fft dimension along y
  ! counter over fft dimension along z
  ! counter over polarizations

  real(DP) :: arg
  ! the argument of the exponent
  ! scalar product of bg and tau
  real(DP) :: bgtau1, bgtau2, bgtau3
  real(DP) :: bg11, bg12, bg13
  real(DP) :: bg21, bg22, bg23
  real(DP) :: bg31, bg32, bg33

  CALL start_clock( 'struct_fact' )

  bg11 = bg(1, 1); bg12 = bg(1, 2); bg13 = bg(1, 3)
  bg21 = bg(2, 1); bg22 = bg(2, 2); bg23 = bg(2, 3)
  bg31 = bg(3, 1); bg32 = bg(3, 2); bg33 = bg(3, 3)

  strf(:,:) = (0.d0,0.d0)
  do nt = 1, ntyp
     do na = 1, nat
        if (ityp (na) .eq.nt) then
           !$cuf kernel do (1) <<<*, *>>>
          do ng = 1, ngm
              arg = (g (1, ng) * tau (1, na) + g (2, ng) * tau (2, na) &
                   + g (3, ng) * tau (3, na) ) * tpi
              strf (ng, nt) = strf (ng, nt) + CMPLX(cos (arg), -sin (arg),kind=DP)
           enddo
        endif
     enddo
  enddo

  !$cuf kernel do (1) <<<*, *>>>
  do na = 1, nat
     bgtau1 = bg11 * tau(1, na) + bg21 * tau(2, na) + bg31 * tau(3, na)
     bgtau2 = bg12 * tau(1, na) + bg22 * tau(2, na) + bg32 * tau(3, na)
     bgtau3 = bg13 * tau(1, na) + bg23 * tau(2, na) + bg33 * tau(3, na)

     do n1 = - nr1, nr1
        arg = tpi * n1 * bgtau1
        eigts1 (n1, na) = CMPLX(cos (arg), - sin (arg) ,kind=DP)
     enddo
     do n2 = - nr2, nr2
        arg = tpi * n2 * bgtau2
        eigts2 (n2, na) = CMPLX(cos (arg), - sin (arg) ,kind=DP)
     enddo
     do n3 = - nr3, nr3
        arg = tpi * n3 * bgtau3
        eigts3 (n3, na) = CMPLX(cos (arg), - sin (arg) ,kind=DP)
     enddo
  enddo

  CALL stop_clock( 'struct_fact' )
  return
end subroutine struc_fact_gpu
#endif

