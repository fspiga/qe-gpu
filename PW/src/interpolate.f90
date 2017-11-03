!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
subroutine interpolate (v, vs, iflag)
  !
  !     This subroutine interpolates :
  !     from the smooth mesh (vs) to a  thicker mesh (v)  (iflag>0)
  !        vs is unchanged on output
  !     from the  thick mesh (v ) to a smoother mesh (vs) (iflag<=0)
  !        v  is unchanged on output
  !     V and Vs are real and in real space . V and Vs may coincide
  !
  USE kinds, ONLY: DP
  USE gvect,  ONLY: nl, nlm
  USE gvecs,ONLY: ngms, nls, nlsm, doublegrid
  USE control_flags, ONLY: gamma_only
  USE fft_base,      ONLY : dfftp, dffts
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  implicit none
  real(DP) :: v (dfftp%nnr), vs (dffts%nnr)
  ! function on thick mesh
  ! function on smooth mesh

  complex(DP), allocatable :: aux (:), auxs (:)
  ! work array on thick mesh
  ! work array on smooth mesh

  integer :: iflag
  ! gives the direction of the interpolation

  integer :: ig, ir

  call start_clock ('interpolate')
  if (iflag <= 0) then
     !
     !    from thick to smooth
     !
     if (doublegrid) then
        allocate (aux( dfftp%nnr))    
        allocate (auxs(dffts%nnr))    
        aux (:) = v (:)
        CALL fwfft ('Dense', aux, dfftp)
        auxs (:) = (0.d0, 0.d0)
        do ig = 1, ngms
           auxs (nls (ig) ) = aux (nl (ig) )
        enddo
        if (gamma_only) then
           do ig = 1, ngms
              auxs (nlsm(ig) ) = aux (nlm(ig) )
           enddo
        end if
        CALL invfft ('Smooth', auxs, dffts)
        vs (:) = auxs (:)
        deallocate (auxs)
        deallocate (aux)
     else
        do ir = 1, dfftp%nnr
           vs (ir) = v (ir)
        enddo
     endif
  else
     !
     !   from smooth to thick
     !
     if (doublegrid) then
        allocate (aux( dfftp%nnr))    
        allocate (auxs(dffts%nnr))    
        auxs (:) = vs (:)
        CALL fwfft ('Smooth', auxs, dffts)
        aux (:) = (0.d0, 0.d0)
        do ig = 1, ngms
           aux (nl (ig) ) = auxs (nls (ig) )
        enddo
        if (gamma_only) then
           do ig = 1, ngms
              aux (nlm(ig) ) = auxs (nlsm(ig) )
           enddo
        end if
        CALL invfft ('Dense', aux, dfftp)
        v (:) = aux (:)
        deallocate (auxs)
        deallocate (aux)
     else
        do ir = 1, dfftp%nnr
           v (ir) = vs (ir)
        enddo
     endif
  endif
  call stop_clock ('interpolate')
  return
end subroutine interpolate
!
subroutine cinterpolate (v, vs, iflag)
  !
  !     This subroutine interpolates :
  !     from the smooth mesh (vs) to a  thicker mesh (v)  (iflag>0)
  !        vs is unchanged on output
  !     from the  thick mesh (v ) to a smoother mesh (vs) (iflag<=0)
  !        v  is unchanged on output
  !     V and Vs are complex and in real space . V and Vs may coincide
  !
  USE kinds, ONLY: DP
  USE gvect,  ONLY: nl, nlm
  USE gvecs,ONLY: ngms, nls, nlsm, doublegrid
  USE control_flags, ONLY: gamma_only
  USE fft_base,      ONLY : dfftp, dffts
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  IMPLICIT NONE
  complex(DP) :: v (dfftp%nnr), vs (dffts%nnr)
  ! function on thick mesh
  ! function on smooth mesh

  integer :: iflag
  ! gives the direction of the interpolation

  complex(DP), allocatable :: aux (:), auxs (:)
  ! work array on thick mesh
  ! work array on smooth mesh

  integer :: ig

  if (gamma_only) call errore ('cinterpolate','not allowed', 1)
  call start_clock ('interpolate')
  if (iflag <= 0) then
     !
     !    from thick to smooth
     !
     if (doublegrid) then
        allocate (aux ( dfftp%nnr))    
        aux (:) = v(:) 
        CALL fwfft ('Dense', aux, dfftp)
        vs (:) = (0.d0, 0.d0)
        do ig = 1, ngms
           vs (nls (ig) ) = aux (nl (ig) )
        enddo
        CALL invfft ('Smooth', vs, dffts)
        deallocate (aux)
     else
        call zcopy (dfftp%nnr, v, 1, vs, 1)
     endif
  else
     !
     !   from smooth to thick
     !
     if (doublegrid) then
        allocate (auxs (dffts%nnr))    
        auxs (:) = vs(:)
        CALL fwfft ('Smooth', auxs, dffts)
        v (:) = (0.d0, 0.d0)
        do ig = 1, ngms
           v (nl (ig) ) = auxs (nls (ig) )
        enddo
        CALL invfft ('Dense', v, dfftp)
        deallocate (auxs)
     else
        call zcopy (dfftp%nnr, vs, 1, v, 1)
     endif
  endif
  call stop_clock ('interpolate')
  return
end subroutine cinterpolate

#ifdef USE_CUDA

subroutine interpolate_gpu (v, vs, iflag)
  !
  !     This subroutine interpolates :
  !     from the smooth mesh (vs) to a  thicker mesh (v)  (iflag>0)
  !        vs is unchanged on output
  !     from the  thick mesh (v ) to a smoother mesh (vs) (iflag<=0)
  !        v  is unchanged on output
  !     V and Vs are real and in real space . V and Vs may coincide
  !
  USE kinds, ONLY: DP
  USE gvect,  ONLY: nl, nl_d, nlm
  USE gvecs,ONLY: ngms, nls, nls_d, nlsm, doublegrid
  USE control_flags, ONLY: gamma_only
  USE fft_base,      ONLY : dfftp, dffts
  USE fft_interfaces,ONLY : fwfft, invfft
  !
  USE cudafor
!  USE ep_debug, ONLY: MPI_Wtime
!  USE ep_debug, ONLY: compare
  !
  implicit none
  real(DP), DEVICE :: v (dfftp%nnr), vs (dffts%nnr)
  ! function on thick mesh
  ! function on smooth mesh

  complex(DP), DEVICE, allocatable :: aux (:), auxs (:)
  ! work array on thick mesh
  ! work array on smooth mesh
!  complex(DP), allocatable :: aux_h (:), auxs_h (:)

  integer :: iflag
  ! gives the direction of the interpolation

  integer :: ig, ir

  real(DP) :: timer

  call start_clock ('interpolate')
!  timer = MPI_Wtime()
  if (iflag <= 0) then
     !
     !    from thick to smooth
     !
     if (doublegrid) then
        allocate (aux( dfftp%nnr))    
        allocate (auxs(dffts%nnr))    

        !allocate (aux_d( dfftp%nnr))
        !allocate (auxs_d(dffts%nnr))

!$cuf kernel do(1) <<<*,*>>>
        do ir = 1, dfftp%nnr
           aux (ir) = v (ir)
        enddo

        !aux_d = aux
        CALL fwfft ('Dense', aux, dfftp)

        !CALL fwfft ('Dense', aux_d, dfftp)

        !call compare(aux, aux_d, "intrp aux fw")

        auxs (:) = (0.d0, 0.d0)
!$cuf kernel do(1) <<<*,*>>>
        do ig = 1, ngms
           auxs (nls_d (ig) ) = aux (nl_d (ig) )
        enddo
#if 0
        if (gamma_only) then
           do ig = 1, ngms
              auxs (nlsm(ig) ) = aux (nlm(ig) )
           enddo
        end if
#endif
        CALL invfft ('Smooth', auxs, dffts)
        !vs (:) = auxs (:)
!$cuf kernel do(1) <<<*,*>>>
        do ir = 1, dffts%nnr
           vs (ir) = auxs (ir)
        enddo

        deallocate (auxs)
        deallocate (aux)
     else
!$cuf kernel do(1) <<<*,*>>>
        do ir = 1, dfftp%nnr
           vs (ir) = v (ir)
        enddo
     endif
  else
     !
     !   from smooth to thick
     !
     if (doublegrid) then
        allocate (aux( dfftp%nnr))    
        allocate (auxs(dffts%nnr))  

!        allocate (aux_h( dfftp%nnr))
!        allocate (auxs_h(dffts%nnr))
 
!$cuf kernel do(1) <<<*,*>>>
        do ir = 1, dffts%nnr
           auxs (ir) = vs (ir)
        enddo
 
!        auxs_h = auxs
! print *,"calling fwfft in interpolate smooth to thick"
        CALL fwfft ('Smooth', auxs, dffts)

!        CALL fwfft ('Smooth', auxs_h, dffts)

!        call compare(auxs_h, auxs, "intrp auxs fw")

        aux (:) = (0.d0, 0.d0)
!        aux_h = aux

!$cuf kernel do(1) <<<*,*>>>
        do ig = 1, ngms
           aux (nl_d (ig) ) = auxs (nls_d (ig) )
        enddo

!        do ig = 1, ngms
!           aux_h (nl (ig) ) = auxs_h (nls (ig) )
!        enddo

!        call compare(aux_h, aux, "intrp aux b4")


#if 0
        if (gamma_only) then
print *,"gamma ONLY!!!"
           do ig = 1, ngms
              aux (nlm(ig) ) = auxs (nlsm(ig) )
           enddo
        end if
#endif
!aux_d = aux

!print *,"input aux(1) = ",aux(1)
! print *,"calling invfft in interpolate smooth to thick"

!        CALL invfft ('Dense', aux_h, dfftp)

        CALL invfft ('Dense', aux, dfftp)


!print *,"cpu aux(1) = ",aux(1)

!        CALL invfft ('Dense', aux_d, dfftp)

!aux = aux_d

!print *,"gpu aux(1) = ",aux(1)


!        call compare(aux_h, aux, "intrp aux fft")

!  deallocate(auxs_h,aux_h)
!        v (:) = aux (:)
!$cuf kernel do(1) <<<*,*>>>
        do ir = 1, dfftp%nnr
           v (ir) = aux (ir)
        enddo

        deallocate (auxs)
        deallocate (aux)
     else
!$cuf kernel do(1) <<<*,*>>>
        do ir = 1, dfftp%nnr
           v (ir) = vs (ir)
        enddo
     endif
  endif
  call stop_clock ('interpolate')
!  timer = MPI_Wtime() - timer
!  print *,"interpolate in ",timer," secs"
  return
end subroutine interpolate_gpu
!

#endif
