!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine deriv_drhoc (ngl, gl, omega, tpiba2, mesh, r, rab, rhoc, drhocg)
  !-----------------------------------------------------------------------
  USE kinds
  USE constants, ONLY : pi, fpi
  implicit none
  !
  !    first the dummy variables
  !

  integer :: ngl, mesh
  ! input: the number of g shell
  ! input: the number of radial mesh points

  real(DP), intent(in) :: gl (ngl), r (mesh), rab (mesh), rhoc (mesh), &
                          omega, tpiba2
  real(DP), intent(out) :: drhocg (ngl)
  ! input: the number of G shells
  ! input: the radial mesh
  ! input: the derivative of the radial mesh
  ! input: the radial core charge
  ! input: the volume of the unit cell
  ! input: 2 times pi / alat
  ! output: fourier transform of d Rho_c/dG
  !
  !     here the local variables
  !
  real(DP) :: gx, rhocg1
  ! the modulus of g for a given shell
  ! the fourier transform
  real(DP), allocatable :: aux (:)
  ! auxiliary memory for integration

  integer :: ir, igl, igl0
  ! counter on radial mesh points
  ! counter on g shells
  ! lower limit for loop on ngl

  !
  ! G=0 term
  !
  if (gl (1) < 1.0d-8) then
     drhocg (1) = 0.0d0
     igl0 = 2
  else
     igl0 = 1
  endif
  !
  ! G <> 0 term
  !
  allocate (aux( mesh))    
  do igl = igl0, ngl
     gx = sqrt (gl (igl) * tpiba2)
     do ir = 1, mesh
        aux (ir) = r (ir) * rhoc (ir) * (r (ir) * cos (gx * r (ir) ) &
             / gx - sin (gx * r (ir) ) / gx**2)
     enddo
     call simpson (mesh, aux, rab, rhocg1)
     drhocg (igl) = fpi / omega * rhocg1
  enddo
  deallocate (aux)

  return
end subroutine deriv_drhoc

#ifdef USE_CUDA
module compute_deriv_drhocg_gpu_m
 contains
  attributes(global) subroutine compute_deriv_rhocg_gpu(ngl, gl, omega, tpiba2, mesh, r, rab, rhoc, drhocg)
      use cudafor
      use kinds
      use constants, ONLY : pi, fpi
      implicit none
      integer, value :: ngl, mesh
  ! input: the number of g shell
  ! input: the number of radial mesh points
      real(DP), intent(in), device :: gl (ngl), r (mesh), rab (mesh), rhoc (mesh)
      real(DP), intent(in), value :: omega, tpiba2
      real(DP), intent(out), device :: drhocg (ngl)

      integer :: tx, ty, igl, ir
      real(DP) :: mysum, val, gx, x

      tx = threadIdx%x
      ty = threadIdx%y

      igl = (blockIdx%x - 1) * blockDim%y + ty

      if (igl > ngl) return
      gx = sqrt(gl(igl) * tpiba2)
      mysum = 0.d0

      do ir = tx, mesh, blockDim%x 
       val = r (ir) * rhoc (ir) * (r (ir) * cos (gx * r (ir) ) &
             / gx - sin (gx * r (ir) ) / gx**2) * &
              rab(ir) 

         if (ir == 1 .or. ir == mesh) then
           mysum = mysum + val
         else if (mod(ir,2)) then
           mysum = mysum + 4.d0*val
         else
           mysum = mysum + 2.d0*val
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
         drhocg(igl) =  (mysum /3.0d0) *fpi /omega
         ! the  G=0 component is not computed
          if (igl ==1 .and. gl(1) < 1.0d-8)   drhocg(1) =  0.d0
       endif  



  end subroutine compute_deriv_rhocg_gpu

end module compute_deriv_drhocg_gpu_m

subroutine deriv_drhoc_gpu (ngl, gl, omega, tpiba2, mesh, r, rab, rhoc, drhocg)
  !-----------------------------------------------------------------------
  !
  USE kinds
  USE constants, ONLY : pi, fpi
  use cudafor
  use compute_rhocg_gpu_m
  use compute_deriv_drhocg_gpu_m
  implicit none
  integer :: ngl, mesh
  ! input: the number of g shell
  ! input: the number of radial mesh points

  real(DP), intent(in), device :: gl(ngl), r(mesh), rab(mesh), rhoc(mesh)
  ! input: the number of G shells
  ! input: the radial mesh
  ! input: the derivative of the radial mesh
  ! input: the radial core charge
  real(DP), intent(in) :: omega, tpiba2
  ! input: the volume of the unit cell
  ! input: 2 times pi / alat
  real(DP), intent(out), device :: drhocg (ngl)
  ! output: fourier transform of d Rho_c/dG
  !
  !     here the local variables
  !
   integer:: blocks
   type(dim3) :: threads
   threads = dim3(32, 8, 1)
   blocks = ceiling(real(ngl)/8)
   call compute_deriv_rhocg_gpu<<<blocks, threads>>>(ngl, gl, omega, tpiba2, mesh, r, rab, rhoc, drhocg)

 end subroutine

#endif
