!
! Copyright (C) 2002 CP90 group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!=======================================================================
!
   subroutine runcg_uspp( nfi, tfirst, tlast, eigr, bec, irb, eigrb, &
      rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac, fion, ema0bg, becdr, &
      lambdap, lambda  )

      use kinds, only: dp
      use control_flags, only: iprint, thdyn, tpre, tbuff, iprsta, trhor, &
            tfor, tvlocw, trhow, taurdr, tprnfor
      use control_flags, only: ndr, ndw, nbeg, nomore, tsde, tortho, tnosee, &
            tnosep, trane, tranp, tsdp, tcp, tcap, ampre, amprp, tnoseh

      use atom, only: nlcc
      use core, only: nlcc_any
!---ensemble-DFT
      use energies, only: eht, epseu, exc, etot, eself, enl, ekin,          &
     &                    atot, entropy, egrand
      use electrons_base, only: f, nspin, nel, iupdwn, nupdwn, nudx, nelt, &
                                nx => nbspx, n => nbsp, ispin

      use ensemble_dft, only: tens, tgrand, ninner, ismear, etemp, ef,       &
     &                tdynz, tdynf, zmass, fmass, fricz, fricf, z0, c0diag,  &
                      becdiag, fmat0, becdrdiag, becm, bec0, fion2, atot0,   &
                      etot0, h0c0, c0hc0, epsi0, e0, dval, z1, f1, dfmat, fmat1, &
                      ef1, enocc, f0, fmatx, fx, zaux, zx, ex, zxt, atot1, etot1, &
                      dedx1, dentdx1, dx, dadx1, faux, eqc, eqa, atotmin, xmin, &
                      entropy2, f2, etot2, eqb, compute_entropy2, compute_entropy_der, &
                      compute_entropy
!---
      use gvecp, only: ngm
      use gvecs, only: ngs
      use gvecb, only: ngb
      use gvecw, only: ngw
      use reciprocal_vectors, only: ng0 => gstart
      use cvan, only: nvb, ish
      use ions_base, only: na, nat, pmass, nax, nsp, rcmax
      use grid_dimensions, only: nnr => nnrx, nr1, nr2, nr3
      use cell_base, only: ainv, a1, a2, a3
      use cell_base, only: omega, alat
      use cell_base, only: h, hold, deth, wmass, tpiba2
      use smooth_grid_dimensions, only: nnrsx, nr1s, nr2s, nr3s
      use smallbox_grid_dimensions, only: nnrb => nnrbx, nr1b, nr2b, nr3b
      use local_pseudo, only: vps, rhops
      use io_global, ONLY: io_global_start, stdout, ionode, ionode_id
      use mp_global, ONLY: mp_global_start
      use para_mod
      use dener
      use derho
      use cdvan
      use stre
      use parameters, only: nacx, natx, nsx, nbndxx
      use constants, only: pi, factem
      use io_files, only: psfile, pseudo_dir
      use io_files, only: outdir

      use uspp, only : nhsa=> nkb, betae => vkb, rhovan => becsum, deeq,qq
      use uspp_param, only: nh
      use cg_module, only: ltresh, itercg, etotnew, etotold, tcutoff, &
          restartcg, passof, passov, passop, ene_ok, numok, maxiter, &
          enever, etresh, ene0, hpsi, gi, hi, esse, essenew, dene0, spasso, &
          ene1, passo, iter3, enesti, ninner_ef, emme
      use ions_positions, only: tau0
      use wavefunctions_module, only: c0, cm, phi => cp
      use efield_module, only: tefield, evalue, ctable, qmat, detq, ipolp, &
            berry_energy, ctabin, gqq, gqqm, df, pberryel
      use mp, only: mp_sum,mp_bcast
      use cp_electronic_mass,       ONLY : emass_cutoff
      use orthogonalize_base,       ONLY : calphi

!
      implicit none
!
      integer :: nfi
      logical :: tfirst , tlast
      complex(dp) :: eigr(ngw,nat)
      real(dp) :: bec(nhsa,n)
      real(dp) :: becdr(nhsa,n,3)
      integer irb(3,nat)
      complex(dp) :: eigrb(ngb,nat)
      real(dp) :: rhor(nnr,nspin)
      real(dp) :: rhog(ngm,nspin)
      real(dp) :: rhos(nnrsx,nspin)
      real(dp) :: rhoc(nnr)
      complex(dp) :: ei1(-nr1:nr1,nat)
      complex(dp) :: ei2(-nr2:nr2,nat)
      complex(dp) :: ei3(-nr3:nr3,nat)
      complex(dp) :: sfac( ngs, nsp )
      real(dp) :: fion(3,natx)
      real(dp) :: ema0bg(ngw)
      real(dp) :: lambdap(nx,nx)
      real(dp) :: lambda(nx,nx)
!
!
      integer :: i, j, ig, k, is, iss,ia, iv, jv, il
      integer :: inl, jnl, niter, istart, nss
      real(dp) :: enb, enbi, x
      complex(dp) :: c2( ngw )
      complex(dp) :: c3( ngw )
      real(dp) :: gamma, entmp, sta
      complex(dp),allocatable :: hpsi0(:,:)
      real(dp) :: sca
      logical  :: newscheme
!
!
      newscheme=.true.

      allocate(hpsi0(ngw,n))
      fion2=0.d0

      if(ionode) open(37,file='convergence.dat',status='unknown')!for debug and tuning purposes
      if(tfirst.and.ionode) write(stdout,*) 'PERFORMING CONJUGATE GRADIENT MINIMIZATION OF EL. STATES'

!set tpa preconditioning

      call  emass_precond_tpa( ema0bg, tpiba2, emass_cutoff )
     
      call prefor(eigr,betae) 

      ltresh    = .false.
      itercg    = 1
      etotold   = 1.d8
      tcutoff   = .false.
      restartcg = .true.
      passof = passop
      ene_ok = .false.


      !orthonormalize c0

      call calbec(1,nsp,eigr,c0,bec)

      call gram(betae,bec,nhsa,c0,ngw,n)

      call calbec(1,nsp,eigr,c0,bec)

      !calculates phi for pcdaga

      ! call calphiid(c0,bec,betae,phi)
      CALL calphi( c0, SIZE(c0,1), bec, nhsa, betae, phi, n )

      !calculates the factors for S and K inversion in US case
      if(nvb.gt.0) then
        call  set_s_minus1(betae)
        call  set_k_minus1(betae, ema0bg)
      endif  


!-------------verifica
!      do i=1,n
!      do ig=1,ngw
!      	phi(ig,i,1,1)=phi(ig,i,1,1)+c0(ig,i,1,1)*(1.d0/ema0bg(ig)-1.d0)
!      enddo
!      enddo
!      !call calbec(1,nsp,eigr,phi,becm)
!      !call sminus1(phi,becm,betae)
!      call kminus1(phi,betae,ema0bg)
!      call calbec(1,nsp,eigr,phi,becm)
!      do i=1,n
!         do j=1,n
!             sca=0.d0
!                do ig=1,ngw
!                 sca=sca+2*DBLE(CONJG(phi(ig,i,1,1))*phi(ig,j,1,1))
!              enddo
!              if (ng0.eq.2) then
!                 sca=sca-DBLE(CONJG(phi(1,i,1,1))*phi(1,j,1,1))
!              endif


!           if (nvb.gt.0) then
!              do is=1,nvb
!                 do iv=1,nh(is)
!                    do jv=1,nh(is)
!                       do ia=1,na(is)
!                          inl=ish(is)+(iv-1)*na(is)+ia
!                          jnl=ish(is)+(jv-1)*na(is)+ia
!                          sca=sca+ qq(iv,jv,is)*becm(inl,i)*becm(jnl,j)
!                       end do
!                    end do
!                 end do
!              end do
!           endif
!            write(6,*) 'VERIFCA S :',i,j,sca
!         enddo
!       enddo






 
      !set index on number of converged iterations

      numok = 0

      do while ( itercg .lt. maxiter .and. (.not.ltresh) )


        ENERGY_CHECK: if(.not. ene_ok ) then
          call calbec(1,nsp,eigr,c0,bec)
          if(.not.tens) then
            call rhoofr(nfi,c0,irb,eigrb,bec,rhovan,rhor,rhog,rhos,enl,ekin)
          else

           if(newscheme) then 
               call  inner_loop( nfi, tfirst, tlast, eigr,  irb, eigrb, &
                 rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,c0,bec  )
           endif
            !     calculation of the rotated quantities
            call rotate(z0,c0(:,:,1,1),bec,c0diag,becdiag)
            !     calculation of rho corresponding to the rotated wavefunctions
            call rhoofr(nfi,c0diag,irb,eigrb,becdiag                         &
                     &                    ,rhovan,rhor,rhog,rhos,enl,ekin)
          endif

          !calculates the potential
          !
          !     put core charge (if present) in rhoc(r)
          !
          if (nlcc_any) call set_cc(irb,eigrb,rhoc)

          !
          !---ensemble-DFT

          if (.not.tens) then

            call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                  &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
            etotnew=etot

          else

            call compute_entropy2( entropy, f, n, nspin )
           
            call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,                &
                     &            ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
            etotnew=etot+entropy

          end if

          if(tefield  ) then!just in this case calculates elfield stuff at zeo field-->to be bettered
            
             call berry_energy( enb, enbi, bec, c0(:,:,1,1), fion )
             etot=etot+enb+enbi
          endif
        else

          etot=enever
          if(.not.tens) then 
             etotnew=etot
          else
             etotnew=etot+entropy
          endif
          ene_ok=.false.

        end if ENERGY_CHECK
        if(ionode) write(37,*)itercg, etotnew,pberryel!for debug and tuning purposes


        

        if(abs(etotnew-etotold).lt.etresh) then
           numok=numok+1
        else 
           numok=0
        endif

        if(numok.ge.4) then
           ltresh=.true.
        endif



        etotold=etotnew
        ene0=etot
        if(tens .and. newscheme) then
          ene0=ene0+entropy
        endif


        !update d

        call newd(rhor,irb,eigrb,rhovan,fion)


        call prefor(eigr,betae)!ATTENZIONE

        do i=1,n,2
          call dforce(bec,betae,i,c0(1,i,1,1),c0(1,i+1,1,1),c2,c3,rhos)
          if(tefield .and. (evalue.ne.0.d0)) then
            call dforceb(c0, i, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            c2(1:ngw)=c2(1:ngw)+evalue*df(1:ngw)
            call dforceb(c0, i+1, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            c3(1:ngw)=c3(1:ngw)+evalue*df(1:ngw)
          endif
          hpsi(1:ngw,  i)=c2(1:ngw)
          hpsi(1:ngw,i+1)=c3(1:ngw)
          if (ng0.eq.2) then
            hpsi(1,  i)=CMPLX(DBLE(hpsi(1,  i)), 0.d0)
            hpsi(1,i+1)=CMPLX(DBLE(hpsi(1,i+1)), 0.d0)
          end if
        enddo

               
        call pcdaga2(c0,phi,hpsi)
               
        hpsi0(1:ngw,1:n)=hpsi(1:ngw,1:n)
        gi(1:ngw,1:n) = hpsi(1:ngw,1:n)
        
        call calbec(1,nsp,eigr,hpsi,becm)
        call sminus1(hpsi,becm,betae)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!look if the following two lines are really needed
        call calbec(1,nsp,eigr,hpsi,becm)
        call pc2(c0,bec,hpsi,becm)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        call kminus1(gi,betae,ema0bg)
        call calbec(1,nsp,eigr,gi,becm)
        call pc2(c0,bec,gi,becm)

        
        if(tens) call calcmt(f,z0,fmat0)

        call calbec(1,nsp,eigr,gi,becm)
        call calbec(1,nsp,eigr,hpsi,bec0) 

!  calculates gamma
        gamma=0.d0
        
        if(.not.tens) then
           
           do i=1,n
              do ig=1,ngw
                 gamma=gamma+2*DBLE(CONJG(gi(ig,i))*hpsi(ig,i))
              enddo
              if (ng0.eq.2) then
                 gamma=gamma-DBLE(CONJG(gi(1,i))*hpsi(1,i))
              endif
           enddo
           
           call mp_sum(gamma)
           
           if (nvb.gt.0) then
              do is=1,nvb
                 do iv=1,nh(is)
                    do jv=1,nh(is)
                       do ia=1,na(is)
                          inl=ish(is)+(iv-1)*na(is)+ia
                          jnl=ish(is)+(jv-1)*na(is)+ia
                          gamma=gamma+ qq(iv,jv,is)*becm(inl,i)*bec0(jnl,i)
                       end do
                    end do
                 end do
              end do
           endif

        else

           do iss=1,nspin
              nss=nupdwn(iss)
              istart=iupdwn(iss)
              do i=1,nss
                 do j=1,nss
                    do ig=1,ngw
                       gamma=gamma+2*DBLE(CONJG(gi(ig,i+istart-1))*hpsi(ig,j+istart-1))*fmat0(j,i,iss)
                    enddo
                    if (ng0.eq.2) then
                       gamma=gamma-DBLE(CONJG(gi(1,i+istart-1))*hpsi(1,j+istart-1))*fmat0(j,i,iss)
                    endif
                    
                 enddo
              enddo
           enddo

           call mp_sum(gamma)
           
           if(nvb.gt.0) then
              do iss=1,nspin
                 nss=nupdwn(iss)
                 istart=iupdwn(iss)
                 do i=1,nss
                    do j=1,nss
                       do is=1,nvb
                          do iv=1,nh(is)
                             do jv=1,nh(is)
                                do ia=1,na(is)
                                   inl=ish(is)+(iv-1)*na(is)+ia
                                   jnl=ish(is)+(jv-1)*na(is)+ia
                                   gamma=gamma+ qq(iv,jv,is)*becm(inl,i+istart-1)*bec0(jnl,j+istart-1)*fmat0(j,i,iss)
                                end do
                             end do
                          end do
                       enddo
                    enddo
                 enddo
              enddo
          endif
        endif



        !case of first iteration

        if(itercg==1.or.(mod(itercg,20).eq.1).or.restartcg) then

          restartcg=.false.
          passof=passop
          hi(1:ngw,1:n)=gi(1:ngw,1:n)!hi is the search direction
          esse=gamma


        else

          !find direction hi for general case 
          !calculates gamma for general case, not using Polak Ribiere
          
          essenew=gamma
          gamma=gamma/esse
          esse=essenew

          hi(1:ngw,1:n)=gi(1:ngw,1:n)+gamma*hi(1:ngw,1:n)

        endif
!note that hi, is saved  on gi, because we need it before projection on conduction states

        !find minimum along direction hi:

        !project hi on conduction sub-space

        call calbec(1,nsp,eigr,hi,bec0)
        call pc2(c0,bec,hi,bec0)
        
        call calbec(1,nsp,eigr,hi,bec0)

        !do quadratic minimization
        !             
        !calculate derivative with respect to  lambda along direction hi

        dene0=0.
        if(.not.tens) then
          do i=1,n               
            do ig=1,ngw
              dene0=dene0-4.d0*DBLE(CONJG(hi(ig,i))*hpsi0(ig,i))
            enddo
            if (ng0.eq.2) then
              dene0=dene0+2.d0*DBLE(CONJG(hi(1,i))*hpsi0(1,i))
            endif
          end do
        else
          !in the ensamble case the derivative is Sum_ij (<hi|H|Psi_j>+ <Psi_i|H|hj>)*f_ji
          !     calculation of the kinetic energy x=xmin      
          call calcmt(f,z0,fmat0)
         do is=1,nspin
             nss=nupdwn(is)
             istart=iupdwn(is)
             do i=1,nss
                do j=1,nss
                   do ig=1,ngw
                    dene0=dene0-2.d0*DBLE(CONJG(hi(ig,i+istart-1))*hpsi0(ig,j+istart-1))*fmat0(j,i,is)
                    dene0=dene0-2.d0*DBLE(CONJG(hpsi0(ig,i+istart-1))*hi(ig,j+istart-1))*fmat0(j,i,is)
                   enddo
                   if (ng0.eq.2) then
                    dene0=dene0+DBLE(CONJG(hi(1,i+istart-1))*hpsi0(1,j+istart-1))*fmat0(j,i,is)
                    dene0=dene0+DBLE(CONJG(hpsi0(1,i+istart-1))*hi(1,j+istart-1))*fmat0(j,i,is)
                   end if
                  enddo
            enddo
          enddo
        endif

        call mp_sum(dene0)

        !if the derivative is positive, search along opposite direction
        if(dene0.gt.0.d0) then
          spasso=-1.D0
        else
          spasso=1.d0
        endif

        !calculates wave-functions on a point on direction hi

        cm(1:ngw,1:n,1,1)=c0(1:ngw,1:n,1,1)+spasso*passof*hi(1:ngw,1:n)


        !orthonormalize

        call calbec(1,nsp,eigr,cm,becm)
        call gram(betae,becm,nhsa,cm,ngw,n)
        call calbec(1,nsp,eigr,cm,becm)
               
        !calculate energy
        if(.not.tens) then
          call rhoofr(nfi,cm,irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,ekin)
        else
          if(newscheme)  call  inner_loop( nfi, tfirst, tlast, eigr,  irb, eigrb, &
           rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm  )

          !     calculation of the rotated quantities
          call rotate(z0,cm(:,:,1,1),becm,c0diag,becdiag)
          !     calculation of rho corresponding to the rotated wavefunctions
          call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,ekin)
        endif

        !calculate potential
        !
        !     put core charge (if present) in rhoc(r)
        !
        if (nlcc_any) call set_cc(irb,eigrb,rhoc)
!
        if (.not.tens) then
          call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                      &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
        else
          call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                      &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
        end if

        if( tefield  ) then!to be bettered
          call berry_energy( enb, enbi, becm, cm(:,:,1,1), fion )
          etot=etot+enb+enbi
        endif
        ene1=etot
        if(tens.and.newscheme) then
          call compute_entropy2( entropy, f, n, nspin )
          ene1=ene1+entropy
        endif
              
            
        !find the minimum

        call minparabola(ene0,spasso*dene0,ene1,passof,passo,enesti)

        if(iprsta.gt.1) write(6,*) ene0,dene0,ene1,passo, gamma, esse

        !set new step

        passov=passof
        passof=2.d0*passo
              
        !calculates wave-functions at minimum

        cm(1:ngw,1:n,1,1)=c0(1:ngw,1:n,1,1)+spasso*passo*hi(1:ngw,1:n)
        if(ng0.eq.2) then
          cm(1,:,1,1)=0.5d0*(cm(1,:,1,1)+CONJG(cm(1,:,1,1)))
        endif

        call calbec(1,nsp,eigr,cm,becm)
        call gram(betae,becm,nhsa,cm,ngw,n)

        !test on energy: check the energy has really diminished

        call calbec(1,nsp,eigr,cm,becm)
        if(.not.tens) then
          call rhoofr(nfi,cm,irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,ekin)
        else
          if(newscheme)  call  inner_loop( nfi, tfirst, tlast, eigr,  irb, eigrb, &
              rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm  )
          !     calculation of the rotated quantities
          call rotate(z0,cm(:,:,1,1),becm,c0diag,becdiag)
          !     calculation of rho corresponding to the rotated wavefunctions
          call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,ekin)
        endif

        !calculates the potential
        !
        !     put core charge (if present) in rhoc(r)
        !
        if (nlcc_any) call set_cc(irb,eigrb,rhoc)
!
        if (.not.tens) then
          call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                       &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
        else
          call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                       &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
        end if
        if( tefield )  then!to be bettered
          call berry_energy( enb, enbi, becm, cm(:,:,1,1), fion )
          etot=etot+enb+enbi
        endif
        enever=etot
        if(tens.and.newscheme) then
          call compute_entropy2( entropy, f, n, nspin )
          enever=enever+entropy
        endif
        if(tens.and.newscheme) then
          if(ionode) write(37,'(a3,4f20.10)') 'CG1',ene0,ene1,enesti,enever
          if(ionode) write(37,'(a3,4f10.7)')  'CG2',spasso,passov,passo,(enever-ene0)/passo/dene0
        else
          if(ionode) write(37,'(a3,4f20.10)') 'CG1',ene0+entropy,ene1+entropy,enesti+entropy,enever+entropy
          if(ionode) write(37,'(a3,4f10.7)')  'CG2',spasso,passov,passo,(enever-ene0)/passo/dene0
        endif
        !check with  what supposed

        if(ionode) then
            if(iprsta.gt.1) then
                 write(stdout,*) 'cg_sub: estimate :'  , (enesti-enever)/(ene0-enever)
                 write(stdout,*) 'cg_sub: minmum   :'  , enever,passo,passov
             endif
        endif

        !if the energy has diminished with respect to  ene0 and ene1 , everything ok
        if( (enever.lt.ene0) .and. (enever.lt.ene1)) then
          c0(:,:,1,1)=cm(:,:,1,1)
          bec(:,:)=becm(:,:)
          ene_ok=.true.
        elseif( (enever.ge.ene1) .and. (enever.lt.ene0)) then
          if(ionode) then
             write(stdout,*) 'cg_sub: missed minimum, case 1, iteration',itercg
          endif
          c0(1:ngw,1:n,1,1)=c0(1:ngw,1:n,1,1)+spasso*passov*hi(1:ngw,1:n)
          restartcg=.true.
          call calbec(1,nsp,eigr,c0,bec)
          call gram(betae,bec,nhsa,c0,ngw,n)
          ene_ok=.false.
          !if  ene1 << energy <  ene0; go to  ene1
        else if( (enever.ge.ene0).and.(ene0.gt.ene1)) then
          if(ionode) then
             write(stdout,*) 'cg_sub: missed minimum, case 2, iteration',itercg
          endif  
          c0(1:ngw,1:n,1,1)=c0(1:ngw,1:n,1,1)+spasso*passov*hi(1:ngw,1:n)
          restartcg=.true.!ATTENZIONE
          call calbec(1,nsp,eigr,c0,bec)
          call gram(betae,bec,nhsa,c0,ngw,n)
          !if ene > ene0,en1 do a steepest descent step
          ene_ok=.false.
        else if((enever.ge.ene0).and.(ene0.le.ene1)) then
        if(ionode) then
             write(stdout,*) 'cg_sub: missed minimum, case 3, iteration',itercg
         endif

          iter3=0
          do while(enever.gt.ene0 .and. iter3.lt.4)
            iter3=iter3+1
            passov=passov*0.5d0
            cm(1:ngw,1:n,1,1)=c0(1:ngw,1:n,1,1)+spasso*passov*hi(1:ngw,1:n)
            ! chenge the searching direction
            spasso=spasso*(-1.d0)
            call calbec(1,nsp,eigr,cm,becm)
            call gram(betae,bec,nhsa,cm,ngw,n)
            call calbec(1,nsp,eigr,cm,becm)
            if(.not.tens) then
              call rhoofr(nfi,cm,irb,eigrb,becm,rhovan,rhor,rhog,rhos,enl,ekin)
            else
              if(newscheme)  call  inner_loop( nfi, tfirst, tlast, eigr,  irb, eigrb, &
              rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,cm,becm  )
              !     calculation of the rotated quantities
              call rotate(z0,cm(:,:,1,1),becm,c0diag,becdiag)
              !     calculation of rho corresponding to the rotated wavefunctions
              call rhoofr(nfi,c0diag,irb,eigrb,becdiag,rhovan,rhor,rhog,rhos,enl,ekin)
            endif
  
            !calculates the potential
            !
            !     put core charge (if present) in rhoc(r)
            !
            if (nlcc_any) call set_cc(irb,eigrb,rhoc)
  !
            if (.not.tens) then
              call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                        &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
                    
            else
              call vofrho(nfi,rhor,rhog,rhos,rhoc,tfirst,tlast,             &
                        &        ei1,ei2,ei3,irb,eigrb,sfac,tau0,fion)
            end if
            if( tefield)  then !to be bettered
              call berry_energy( enb, enbi, becm, cm(:,:,1,1), fion )
              etot=etot+enb+enbi
            endif
            enever=etot
           if(tens.and.newscheme) then
             call compute_entropy2( entropy, f, n, nspin )
             enever=enever+entropy
           endif

          end do
  
          c0(:,:,1,1)=cm(:,:,1,1)
          restartcg=.true.
          ene_ok=.false.
        end if
        
        if(tens.and.newscheme) enever=enever-entropy
 
        call calbec (1,nsp,eigr,c0,bec)
        !calculates phi for pc_daga
        !call calphiid(c0,bec,betae,phi)
        CALL calphi( c0, SIZE(c0,1), bec, nhsa, betae, phi, n )
  
        !=======================================================================
        !
        !                 start of the inner loop
        !                 (Uij degrees of freedom)
        !
        !=======================================================================
        if(tens.and. .not.newscheme) call  inner_loop( nfi, tfirst, tlast, eigr,  irb, eigrb, &
           rhor, rhog, rhos, rhoc, ei1, ei2, ei3, sfac,c0,bec  )

 
            
      
          !=======================================================================
          !                 end of the inner loop
          !=======================================================================
  
  
        itercg=itercg+1

!   restore hi
!        hi(:,:)=gi(:,:) 
 
      end do!on conjugate gradient iterations
      !calculates atomic forces and lambda
      call newd(rhor,irb,eigrb,rhovan,fion)
      if (.not.tens) then
        if (tfor .or. tprnfor) call nlfq(c0,eigr,bec,becdr,fion)
      else
        if (tfor .or. tprnfor) call nlfq(c0diag,eigr,becdiag,becdrdiag,fion)
      endif
  
        call prefor(eigr,betae)
        do i=1,n,2
          call dforce(bec,betae,i,c0(1,i,1,1),c0(1,i+1,1,1),c2,c3,rhos)
          if(tefield.and.(evalue .ne. 0.d0)) then
            call dforceb &
               (c0, i, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            do ig=1,ngw
              c2(ig)=c2(ig)+evalue*df(ig)
            enddo
            call dforceb &
               (c0, i+1, betae, ipolp, bec ,ctabin(1,1,ipolp), gqq, gqqm, qmat, deeq, df)
            do ig=1,ngw
              c3(ig)=c3(ig)+evalue*df(ig)
            enddo
          endif
          do ig=1,ngw
            gi(ig,  i)=c2(ig)
            gi(ig,i+1)=c3(ig)
          end do
          if (ng0.eq.2) then
            gi(1,  i)=CMPLX(DBLE(gi(1,  i)),0.d0)
            gi(1,i+1)=CMPLX(DBLE(gi(1,i+1)),0.d0)
          end if
        enddo
        do i=1,n
          do j=i,n
            lambda(i,j)=0.d0
            do ig=1,ngw
              lambda(i,j)=lambda(i,j)-2.d0*DBLE(CONJG(c0(ig,i,1,1))*gi(ig,j))
            enddo
            if(ng0.eq.2) then
              lambda(i,j)=lambda(i,j)+DBLE(CONJG(c0(1,i,1,1))*gi(1,j))
            endif
            lambda(j,i)=lambda(i,j)
          enddo
        enddo
  
        call mp_sum(lambda)
  
        if(tens) then!in the ensemble case matrix labda must be multiplied with f
          do i=1,n
            do j=1,n
              lambdap(i,j)=0.d0
              do k=1,n
                lambdap(i,j)=lambdap(i,j)+lambda(i,k)*fmat0(k,j,1)
              end do
            end do
          enddo
          do i=1,n
            do j=1,n
              sta=lambda(i,j)
              lambda(i,j)=lambdap(i,j)
              lambdap(i,j)=sta
            enddo
          enddo
        call nlsm2(ngw,nhsa,n,eigr,c0(:,:,1,1),becdr,.true.)
        endif
        call nlfl(bec,becdr,lambda,fion)
          
        ! bforceion adds the force term due to electronic berry phase
        ! only in US-case
          
        if( tefield.and.(evalue .ne. 0.d0) ) then
           call bforceion(fion,tfor.or.tprnfor,ipolp, qmat,bec,becdr,gqq,evalue)

        endif

        deallocate( hpsi0) 
       if(ionode) close(37)!for debug and tuning purposes
END SUBROUTINE
