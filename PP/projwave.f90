!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine do_projwfc (nodenumber)
  !-----------------------------------------------------------------------
  !
  ! projects wavefunctions onto orthogonalized atomic wavefunctions
  ! calculates Lowdin charges, spilling parameter
  ! input: namelist "&inputpp", with variables
  ! prefix      prefix of input files saved by program pwscf
  ! outdir      temporary directory where files resides
  !
  use pwcom
  use io
#ifdef __PARA
  use para, only: me
  use mp
#endif
  implicit none
  character (len=3)  :: nodenumber
  character (len=8)  :: io_choice
  character(len=256) :: outdir
  real (kind=DP)     :: Emin, Emax, DeltaE, smoothing
  integer :: ios, ionode_id = 0 
  namelist / inputpp / outdir, prefix, io_choice, &
             Emin, Emax, DeltaE, smoothing
  !
  nd_nmbr = nodenumber
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  outdir = './'
  io_choice ='both'
  Emin   =-1000000.
  Emax   =+1000000.
  DeltaE = 0.01
  smoothing  = 0.d0
  !
#ifdef __PARA
  if (me == 1)  then
#endif
  read (5, inputpp, err = 200, iostat = ios)
200 call errore ('projwave', 'reading inputpp namelist', abs (ios) )
  !
  tmp_dir = trim(outdir)
  !
#ifdef __PARA
  end if
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id )
  CALL mp_bcast( prefix, ionode_id )
  CALL mp_bcast( io_choice, ionode_id )
  CALL mp_bcast( smoothing, ionode_id )
  CALL mp_bcast( DeltaE, ionode_id )
  CALL mp_bcast( Emin, ionode_id )
  CALL mp_bcast( Emax, ionode_id )
#endif
  if ( smoothing .lt. DeltaE ) smoothing= DeltaE
  if (io_choice.ne.'standard' .and. io_choice.ne.'files' .and.  &
      io_choice.ne.'both') &
      call errore ('projwave','io_choice definition is invalid',1)

  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  call read_file
  call openfil
  !
  call projwave (io_choice,Emin, Emax, DeltaE, smoothing)
  !
  return
end subroutine do_projwfc

!-----------------------------------------------------------------------
subroutine projwave (io_choice,Emin, Emax, DeltaE, smoothing)
  !-----------------------------------------------------------------------
  !
#include "machine.h"
  use pwcom
  use becmod
  use io
#ifdef __PARA
  use para
#endif
  implicit none
  character (len=80) :: filproj
  character (len=8)  :: io_choice
  character (len=33) :: filextension
  character (len=1)  :: l_label(0:3)=(/'s','p','d','f'/)
  !
  type wfc_label
     integer na, n, l, m
  end type wfc_label
  type(wfc_label), allocatable :: nlmchi(:)
  !
  integer :: ik, ibnd, i, j, k, na, nb, nt, isym, n,  m, m1, l, lm, nwfc,&
       nwfc1, lmax_wfc, c_tab, ne, ie_mid, ie_delta, ie, is
  logical :: exst
  real(kind=DP) :: psum, totcharge, Emin, Emax, DeltaE, smoothing, etev, &
             delta, w0gauss, Elw, Eup
  real(kind=DP), allocatable :: e (:), proj (:,:,:), charges(:,:), pdos(:,:,:)
  complex(kind=DP), allocatable :: wfcatom (:,:),  overlap (:,:), &
       work (:,:), work1(:), proj0(:,:)
  integer, allocatable :: index(:)
  external w0gauss
  !
  !
  write (6, '(/5x,"Calling projwave .... ")')
  if (io_choice.eq.'standard' ) &
     write (6, '(5x,"Projections are written on standard output")')
  if (io_choice.eq.'files' ) &
     write (6, '(5x,"Projections are written on files")')
  if (io_choice.eq.'both' ) &
     write (6, '(5x,"Projections are written on both standard output and file")')

  !
  allocate(swfcatom (npwx , natomwfc ) )
  allocate(wfcatom (npwx, natomwfc) )
  allocate(proj (natomwfc, nbnd, nkstot) )
  allocate(overlap (natomwfc, natomwfc) )
  allocate(e (natomwfc) )

  proj   = 0.d0
  overlap= (0.d0,0.d0)
  !
  ! initialize D_Sl for l=1, l=2 and l=3, for l=0 D_S0 is 1
  !
  call d_matrix (d1, d2, d3)  
  !
  ! fill structure nlmchi
  !
  allocate (nlmchi(natomwfc))
  nwfc=0
  lmax_wfc = 0
  do na = 1, nat
     nt = ityp (na)
     do n = 1, nchi (nt)
        if (oc (n, nt) .gt.0.d0.or..not.newpseudo (nt) ) then
           l = lchi (n, nt)
           lmax_wfc = max (lmax_wfc, l )
           do m = 1, 2 * l + 1
              nwfc=nwfc+1
              nlmchi(nwfc)%na = na
              nlmchi(nwfc)%n  =  n
              nlmchi(nwfc)%l  =  l
              nlmchi(nwfc)%m  =  m
           enddo
        endif
     enddo
  enddo

  if (lmax_wfc.gt.3) call errore ('projwave', 'l > 3 not yet implemented', 1)
  if (nwfc.ne.natomwfc) call errore ('projwave', 'wrong # of atomic wfcs?', 1)
  !
  !    loop on k points
  !
  call init_us_1
  !
  do ik = 1, nks
     npw = npwx
     call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     call davcio (evc, nwordwfc, iunwfc, ik, - 1)

     call atomic_wfc (ik, wfcatom)
     call init_us_2 (npw, igk, xk (1, ik), vkb)

     call ccalbec (nkb, npwx, npw, natomwfc, becp, vkb, wfcatom)

     call s_psi (npwx, npw, natomwfc, wfcatom, swfcatom)
     !
     ! wfcatom = |phi_i> , swfcatom = \hat S |phi_i>
     ! calculate overlap matrix O_ij = <phi_i|\hat S|\phi_j>
     !
     call ZGEMM ('c', 'n', natomwfc, natomwfc, npw, (1.d0, 0.d0) , &
          wfcatom, npwx, swfcatom, npwx, (0.d0, 0.d0) , overlap, natomwfc)
#ifdef __PARA
     call reduce (2 * natomwfc * natomwfc, overlap)
#endif
     !
     ! calculate O^{-1/2}
     !
     allocate(work (natomwfc, natomwfc) )
     call cdiagh (natomwfc, overlap, natomwfc, e, work)
     do i = 1, natomwfc
        e (i) = 1.d0 / dsqrt (e (i) )
     enddo
     do i = 1, natomwfc
        do j = i, natomwfc
           overlap (i, j) = (0.d0, 0.d0)
           do k = 1, natomwfc
              overlap (i, j) = overlap (i, j) + e (k) * work (j, k) * conjg (work (i, k) )
           enddo
           if (j.ne.i) overlap (j, i) = conjg (overlap (i, j))
        enddo
     enddo
     deallocate (work)
     !
     ! calculate wfcatom = O^{-1/2} \hat S | phi>
     !
     call ZGEMM ('n', 't', npw, natomwfc, natomwfc, (1.d0, 0.d0) , &
          swfcatom, npwx,  overlap, natomwfc, (0.d0, 0.d0), wfcatom, npwx)
     !
     ! make the projection <psi_i| O^{-1/2} \hat S | phi_j>
     !
     allocate(proj0(natomwfc,nbnd) )
     call ZGEMM ('c', 'n', natomwfc, nbnd, npw, (1.d0, 0.d0) , &
          wfcatom, npwx, evc, npwx, (0.d0, 0.d0) , proj0, natomwfc)
#ifdef __PARA
     call reduce (2 * natomwfc * nbnd, proj0)
#endif
     !
     ! symmetrize the projections
     !
     allocate(work1 (nbnd) )
     do nwfc = 1, natomwfc
        !
        !  atomic wavefunction nwfc is on atom na
        !
        na= nlmchi(nwfc)%na
        n = nlmchi(nwfc)%n
        l = nlmchi(nwfc)%l
        m = nlmchi(nwfc)%m
        !
        do isym = 1, nsym
           nb = irt (isym, na)
           do nwfc1 =1, natomwfc
              if (nlmchi(nwfc1)%na.eq. nb             .and. &
                   nlmchi(nwfc1)%n .eq. nlmchi(nwfc)%n .and. &
                   nlmchi(nwfc1)%l .eq. nlmchi(nwfc)%l .and. &
                   nlmchi(nwfc1)%m .eq. 1 ) go to 10
           end do
           call errore('projwave','cannot symmetrize',1)
10         nwfc1=nwfc1-1
           !
           !  nwfc1 is the first rotated atomic wfc corresponding to nwfc
           !
           if (l.eq.0) then
              work1(:) = proj0 (nwfc1 + 1,:)
           else if (l.eq.1) then 
              work1(:) = 0.d0  
              do m1 = 1, 3  
                 work1(:) = work1(:) + d1 (m1, m, isym) * proj0 (nwfc1 + m1,:)
              enddo
           else if (l.eq.2) then 
              work1(:) = 0.d0  
              do m1 = 1, 5  
                 work1(:) = work1(:) + d2 (m1, m, isym) * proj0 (nwfc1 + m1,:)
              enddo
           else if (l.eq.3) then 
              work1(:) = 0.d0  
              do m1 = 1, 7  
                 work1(:) = work1(:) + d3 (m1, m, isym) * proj0 (nwfc1 + m1,:)
              enddo
           endif
           do ibnd = 1, nbnd
              proj (nwfc, ibnd, ik) = proj (nwfc, ibnd, ik) + &
                   work1(ibnd) * conjg (work1(ibnd)) / nsym
           enddo
        enddo
     enddo
     deallocate (work1)
     deallocate (proj0 )
     ! on k-points
  enddo
#ifdef __PARA
  !
  !   recover the vector proj over the pools
  !
  call poolrecover (et, nbnd, nkstot, nks)
  call poolrecover (proj, nbnd * natomwfc, nkstot, nks)
  !
  if (me.eq.1.and.mypool.eq.1) then
#endif
     !
     ! write on the standard output file
     !
     if (io_choice.eq.'standard' .or.io_choice.eq.'both' ) then 
        write(6,'(/"Projection on atomic states:"/)')
        do nwfc = 1, natomwfc
           write(6,'(5x,"state #",i3,": atom ",i3," (",a3,"), wfc ",i2, &
                &       " (l=",i1," m=",i2,")")') &
                nwfc, nlmchi(nwfc)%na, atm(ityp(nlmchi(nwfc)%na)), &
                nlmchi(nwfc)%n, nlmchi(nwfc)%l, nlmchi(nwfc)%m
        end do
        !
        allocate(index (natomwfc) )
        do ik = 1, nkstot
           write (6, '(/" k = ",3f14.10)') (xk (i, ik) , i = 1, 3)
           do ibnd = 1, nbnd
              write (6, '(5x,"e = ",f14.10," eV")') et (ibnd, ik) * rytoev
              !
              ! sort projections by magnitude, in decreasing order
              !
              do nwfc = 1, natomwfc
                 index (nwfc) = 0
                 e (nwfc) = - proj (nwfc, ibnd, ik)
              end do
              call hpsort (natomwfc, e, index)
              !
              !  only projections that are larger than 0.001 are written
              !
              do nwfc = 1, natomwfc
                 e (nwfc) = - e(nwfc)
                 if ( abs (e(nwfc)).lt.0.001 ) go to 20
              end do
              nwfc = natomwfc + 1
20            nwfc = nwfc -1
              !
              ! fancy (?!?) formatting
              !
              write (6, '(5x,"psi = ",5(f5.3,"*[#",i3,"]+"))') &
                  (e (i), index(i), i = 1, min(5,nwfc))
              do j = 1, (nwfc-1)/5
                 write (6, '(10x,"+",5(f5.3,"*[#",i3,"]+"))') &
                   (e (i), index(i), i = 5*j+1, min(5*(j+1),nwfc))
              end do
              psum = 0.d0
              do nwfc = 1, natomwfc
                 psum = psum + proj (nwfc, ibnd, ik)
              end do
              write (6, '(4x,"|psi|^2 = ",f5.3)') psum
              !
           enddo
        enddo
        deallocate (index)
     end if
     !
     ! estimate partial charges (Loewdin) on each atom
     !
     allocate ( charges (nat, 0:lmax_wfc ) )
     charges = 0.0
     do ik = 1, nkstot
        do ibnd = 1, nbnd
           do nwfc = 1, natomwfc
              na= nlmchi(nwfc)%na
              l = nlmchi(nwfc)%l
              charges(na,l) = charges(na,l) + wg (ibnd,ik) * &
                              proj (nwfc, ibnd, ik)
           enddo
        end do
     end do
     !
     write (6, '(/"Lowdin Charges: "/)')
     !
     psum = 0.0
     do na = 1, nat
        totcharge = 0.d0
        do l = 0, lmax_wfc
           totcharge = totcharge + charges(na,l)
        end do
        psum = psum + totcharge
        write (6, '(5x,"Atom # ",i3,": total charge = ",f8.4, &
      &                ", s, p, d, f = ",4f8.4    )') &
             na, totcharge, ( charges(na,l), l= 0,lmax_wfc)
     end do
     psum = psum / nelec
     write (6, '(5x,"Spilling Parameter: ",f8.4)') 1.0 - psum
     !
     ! Sanchez-Portal et al., Sol. State Commun.  95, 685 (1995).
     ! The spilling parameter measures the ability of the basis provided by
     ! the pseudo-atomic wfc to represent the PW eigenstates,
     ! by measuring how much of the subspace of the Hamiltonian
     ! eigenstates falls outside the subspace spanned by the atomic basis
     !
     deallocate (charges)

     if (io_choice.eq.'files' .or. io_choice.eq.'both') then
        !
        ! find band extrema
        !
        Elw = et (1, 1)
        Eup = et (nbnd, 1)
        do ik = 2, nkstot
           Elw = min (Elw, et (1, ik) )
           Eup = max (Eup, et (nbnd, ik) )
        enddo
        Emin = max (Emin, Elw*rytoev - 5*smoothing )
        Emax = min (Emax, Eup*rytoev + 5*smoothing )

        ne = nint( (Emax-Emin)/DeltaE )
        
        allocate (pdos(0:ne,0:natomwfc+1,nspin))
        pdos(:,:,:) = 0.d0
        current_spin = 1
        ie_delta = 5 * smoothing / DeltaE+1
        do ik = 1,nkstot
           if ( nspin.eq.2 ) current_spin = isk ( ik )
           do ibnd = 1, nbnd
              etev = et(ibnd,ik) * rytoev
              ie_mid = nint( (etev-Emin)/DeltaE )
              do ie = max(ie_mid-ie_delta, 0), min(ie_mid+ie_delta, ne)
                 delta = w0gauss((Emin+DeltaE*ie-etev)/smoothing,0)/smoothing
                 do nwfc = 1, natomwfc
                    pdos(ie,nwfc,current_spin) = pdos(ie,nwfc,current_spin) + &
                                 wk(ik) * delta * proj (nwfc, ibnd, ik) 
                 end do
                 pdos(ie,0,current_spin) = pdos(ie,0,current_spin) + &
                                 wk(ik) * delta
              end do
           end do
        end do

        do is=1,nspin
           do ie=0,ne
              pdos(ie,natomwfc+1,is) = sum(pdos(ie,1:natomwfc,is))
           end do
        end do

        do nwfc = 1, natomwfc
           if (nlmchi(nwfc)%m .eq. 1) then
              filextension='.pdos_atm#'
        !                   12345678901

              c_tab = 11
              if (nlmchi(nwfc)%na.lt.10) then
                 write (filextension( c_tab : c_tab ),'(i1)') nlmchi(nwfc)%na
                 c_tab = c_tab + 1
              else if (nlmchi(nwfc)%na.lt.100) then
                 write (filextension( c_tab : c_tab+1 ),'(i2)') nlmchi(nwfc)%na
                 c_tab = c_tab + 2
              else if (nlmchi(nwfc)%na.lt.1000) then
                 write (filextension( c_tab : c_tab+2 ),'(i3)') nlmchi(nwfc)%na
                 c_tab = c_tab + 3
              else
                 call errore('projwave',&
                             'file extension not supporting so many atoms', &
                              nwfc)
              endif
              write (filextension(c_tab:c_tab+4),'(a1,a)') &
                    '(',trim(atm(ityp(nlmchi(nwfc)%na)))
              c_tab = c_tab + len_trim(atm(ityp(nlmchi(nwfc)%na))) + 1
              if (nlmchi(nwfc)%n.ge.10) &
                 call errore('projwave',&
                             'file extension not supporting so many atmic wfc',&
                              nwfc)
              if (nlmchi(nwfc)%l.gt.3) &
                 call errore('projwave',&
                             'file extension not supporting so many l', &
                              nwfc)
              write (filextension(c_tab:),'(")_wfc#",i1,"(",a1,")")')  &
                    nlmchi(nwfc)%n, l_label(nlmchi(nwfc)%l)
              open (4,file=trim(prefix)//filextension,form='formatted', &
                      status='unknown')

              write (4,'("# E (eV) ",$)')
              do m=1,2 * nlmchi(nwfc)%l + 1
                 if (nspin.eq.1) then
                    write(4,'(" dos(E)    ",$)') 
                 else
                    write(4,'(" dosup(E)  ",$)') 
                    write(4,'(" dosdw(E)  ",$)') 
                 end if
              end do
              write(4,*)

              do ie= 0, ne
                 etev = Emin + ie * DeltaE
                 write (4,'(f7.3,14e11.3)') etev,  &
                       ((pdos(ie,nwfc+m-1,is), is=1,nspin), &
                                               m=1,2*nlmchi(nwfc)%l+1)
              end do
              close (4)
           end if
        end do
        open (4,file=trim(prefix)//".pdos_tot",form='formatted', &
                status='unknown')
        if (nspin.eq.1) then
           write (4,'("# E (eV)  dos(E)    pdos(E)")') 
        else
           write (4,'("# E (eV)  dosup(E)   dosdw(E)  pdosup(E)  pdosdw(E)")') 
        end if
        do ie= 0, ne
           etev = Emin + ie * DeltaE
           write (4,'(f7.3,4e11.3)') etev, (pdos(ie,0,is), is=1,nspin), &
                 (pdos(ie,natomwfc+1,is), is=1,nspin)
        end do
        close (4)
        deallocate (pdos)
     end if

#ifdef __PARA
  endif
#endif
  deallocate (nlmchi)
  deallocate (e)
  deallocate (overlap)
  deallocate (proj)
  deallocate (wfcatom)
  deallocate (swfcatom)

  return

end subroutine projwave
