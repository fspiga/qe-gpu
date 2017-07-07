!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
SUBROUTINE simpson(mesh, func, rab, asum)
  !-----------------------------------------------------------------------
  !
  !     simpson's rule integration. On input:
  !       mesh = the number of grid points (should be odd)
  !       func(i)= function to be integrated
  !       rab(i) = r(i) * dr(i)/di * di
  !     For the logarithmic grid not including r=0 :
  !       r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
  !     For the logarithmic grid including r=0 :
  !       r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
  !     Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr
  !     where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
  !
  USE kinds, ONLY: DP
  IMPLICIT NONE
  INTEGER, INTENT(in) :: mesh
  real(DP), INTENT(in) :: rab (mesh), func (mesh)
  real(DP), INTENT(out):: asum
  !
  real(DP) :: f1, f2, f3, r12
  INTEGER :: i
  !
  asum = 0.0d0
  r12 = 1.0d0 / 3.0d0
  f3 = func (1) * rab (1) * r12

  DO i = 2, mesh - 1, 2
     f1 = f3
     f2 = func (i) * rab (i) * r12
     f3 = func (i + 1) * rab (i + 1) * r12
     asum = asum + f1 + 4.0d0 * f2 + f3
  ENDDO
  !
  ! if mesh is not odd, use open formula instead:
  ! ... 2/3*f(n-5) + 4/3*f(n-4) + 13/12*f(n-3) + 0*f(n-2) + 27/12*f(n-1)
  !!! Under testing
  !
  !IF ( MOD(mesh,2) == 0 ) THEN
  !   print *, 'mesh even: correction:', f1*5.d0/4.d0-4.d0*f2+23.d0*f3/4.d0, &
  !                                      func(mesh)*rab(mesh), asum
  !   asum = asum + f1*5.d0/4.d0 - 4.d0*f2 + 23.d0*f3/4.d0
  !END IF

  RETURN
END SUBROUTINE simpson
#ifdef USE_CUDA
SUBROUTINE simpson_gpu(mesh, func_d, rab_d, asum)
  !-----------------------------------------------------------------------
  !
  !     simpson's rule integration. On input:
  !       mesh = the number of grid points (should be odd)
  !       func(i)= function to be integrated
  !       rab(i) = r(i) * dr(i)/di * di
  !     For the logarithmic grid not including r=0 :
  !       r(i) = r_0*exp((i-1)*dx) ==> rab(i)=r(i)*dx
  !     For the logarithmic grid including r=0 :
  !       r(i) = a(exp((i-1)*dx)-1) ==> rab(i)=(r(i)+a)*dx
  !     Output in asum = \sum_i c_i f(i)*rab(i) = \int_0^\infty f(r) dr
  !     where c_i are alternativaly 2/3, 4/3 except c_1 = c_mesh = 1/3
  !
  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  INTEGER, INTENT(in) :: mesh
  real(DP), INTENT(in), device :: rab_d (mesh), func_d (mesh)
  real(DP), INTENT(out):: asum
  !
  real(DP) :: f1, f2, f3, r12
  INTEGER :: i
  !
  asum = 0.0d0
  r12 = 1.0d0 / 3.0d0

  !$cuf kernel do (1) <<<*, *>>>
  do i = 2, mesh-1, 2
    asum = asum + func_d(i-1)*rab_d(i-1) + 4.0d0*func_d(i)*rab_d(i) + func_d(i+1)*rab_d(i+1)
  end do

  asum = asum*r12
  !
  ! if mesh is not odd, use open formula instead:
  ! ... 2/3*f(n-5) + 4/3*f(n-4) + 13/12*f(n-3) + 0*f(n-2) + 27/12*f(n-1)
  !!! Under testing
  !
  !IF ( MOD(mesh,2) == 0 ) THEN
  !   print *, 'mesh even: correction:', f1*5.d0/4.d0-4.d0*f2+23.d0*f3/4.d0, &
  !                                      func(mesh)*rab(mesh), asum
  !   asum = asum + f1*5.d0/4.d0 - 4.d0*f2 + 23.d0*f3/4.d0
  !END IF

  RETURN
END SUBROUTINE simpson_gpu
#endif

!=-----------------------------------------------------------------------
SUBROUTINE simpson_cp90( mesh, func, rab, asum )
  !-----------------------------------------------------------------------
  !
  !    This routine computes the integral of a function defined on a
  !    logaritmic mesh, by using the open simpson formula given on
  !    pag. 109 of Numerical Recipes. In principle it is used to
  !    perform integrals from zero to infinity. The first point of
  !    the function should be the closest to zero but not the value
  !    in zero. The formula used here automatically includes the
  !    contribution from the zero point and no correction is required.
  !
  !    Input as "simpson". At least 8 integrating points are required.
  !
  !    last revised 12 May 1995 by Andrea Dal Corso
  !
  USE kinds, ONLY: DP
  IMPLICIT NONE
  INTEGER, INTENT(in) :: mesh
  real(DP), INTENT(in) :: rab (mesh), func (mesh)
  real(DP), INTENT(out):: asum
  !
  real(DP) :: c(4)
  INTEGER ::i
  !
  IF ( mesh < 8 ) CALL errore ('simpson_cp90','few mesh points',8)

  c(1) = 109.0d0 / 48.d0
  c(2) = -5.d0 / 48.d0
  c(3) = 63.d0 / 48.d0
  c(4) = 49.d0 / 48.d0

  asum = ( func(1)*rab(1) + func(mesh  )*rab(mesh  ) )*c(1) &
       + ( func(2)*rab(2) + func(mesh-1)*rab(mesh-1) )*c(2) &
       + ( func(3)*rab(3) + func(mesh-2)*rab(mesh-2) )*c(3) &
       + ( func(4)*rab(4) + func(mesh-3)*rab(mesh-3) )*c(4)
  DO i=5,mesh-4
     asum = asum + func(i)*rab(i)
  ENDDO

  RETURN
END SUBROUTINE simpson_cp90
!
!-----------------------------------------------------------------------
SUBROUTINE herman_skillman_int(mesh,func,rab,asum)
!-----------------------------------------------------------------------
  !     simpson rule integration for herman skillman mesh (obsolescent)
  !     Input as in "simpson". BEWARE: "func" is overwritten!!!
  !
  USE kinds, ONLY: DP
  IMPLICIT NONE
  INTEGER, INTENT(in) :: mesh
  real(DP), INTENT(in) :: rab (mesh)
  real(DP), INTENT(inout) :: func (mesh)
  real(DP), INTENT(out):: asum
  !
  INTEGER :: i, j, k, i1, nblock
  REAL(DP) :: a1, a2e, a2o, a2es
  !
  a1=0.0d0
  a2e=0.0d0
  asum=0.0d0
  nblock=mesh/40
  i=1
  func(1)=0.0d0
  DO j=1,nblock
     DO k=1,20
        i=i+2
        i1=i-1
        a2es=a2e
        a2o=func(i1)/12.0d0
        a2e=func(i)/12.0d0
        a1=a1+5.0d0*a2es+8.0d0*a2o-a2e
        func(i1)=asum+a1*rab(i1)
        a1=a1-a2es+8.0d0*a2o+5.0d0*a2e
        func(i)=asum+a1*rab(i)
     ENDDO
     asum=func(i)
     a1=0.0d0
  ENDDO
  !
  RETURN
END SUBROUTINE herman_skillman_int

#ifdef USE_CUDA
module simpson_gpu_m
  contains
    attributes(global) subroutine simpson_gpu_1_kernel(n, tpiba2, omega, gl, r, rhoc, rab, mesh, rhocg)
      use cudafor
      use kinds
      use constants, ONLY : pi, fpi, eps14
      implicit none

      integer, value :: n, mesh
      real(DP), value :: tpiba2, omega
      real(DP), device, intent(in) :: gl(n), r(mesh), rhoc(mesh), rab(mesh)
      real(DP), device, intent(out) :: rhocg(n)

      integer :: tx, ty, igl, ir
      real(DP) :: mysum, val, gx, x

      tx = threadIdx%x
      ty = threadIdx%y

      igl = (blockIdx%x - 1) * blockDim%y + ty

      if (igl > n) return

      gx = sqrt(gl(igl) * tpiba2)
      mysum = 0.d0


      if (abs(gx) < eps14) then
        do ir = tx, mesh, blockDim%x 
          val = r(ir) * r(ir) * rhoc(ir) * rab(ir)

          if (ir == 1 .or. ir == mesh) then
            mysum = mysum + val
          else if (mod(ir,2)) then
            mysum = mysum + 4.d0*val
          else
            mysum = mysum + 2.d0*val
          endif
        end do
      else
        do ir = tx, mesh, blockDim%x 
          x = gx * r(ir)
          if (abs(x) > 0.5_DP) then
            val = sin (x) / (x) * r(ir) * r(ir) * rhoc(ir)
          else
            val = ( 1.0_dp - x*x/6.0_dp * &
                  ( 1.0_dp - x*x/20.0_dp * &
                  ( 1.0_dp - x*x/42.0_dp * &
                  ( 1.0_dp - x*x/72.0_dp ) ) ) ) * r(ir) * r(ir) * rhoc(ir)
          endif

          val = val * rab(ir)

          if (ir == 1 .or. ir == mesh) then
            mysum = mysum + val
          else if (mod(ir,2)) then
            mysum = mysum + 4.d0*val
          else
            mysum = mysum + 2.d0*val
          endif
        end do
      endif

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
        rhocg(igl) = fpi * mysum / (3.d0 * omega)
      endif

    end subroutine simpson_gpu_1_kernel
end module
#endif
