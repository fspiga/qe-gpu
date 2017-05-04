#ifdef TRACK_FLOPS
module flops_tracker
  implicit none
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  REAL(DP) :: fft_ops=0.d0

end module flops_tracker
#endif

#ifdef USE_CUDA
module mpiDeviceUtil
  implicit none
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif

  interface
     subroutine quicksort(base, nmemb, elemsize, compar) &
          bind(C,name='qsort')
       use iso_c_binding
       implicit none
       !pgi$ ignore_tkr base,nmemb,elemsize,compar
       type(C_PTR), value :: base
       integer(C_SIZE_T), value :: nmemb, elemsize
       type(C_FUNPTR), value :: compar
     end subroutine quicksort

     integer function strcmp(a,b) bind(C,name='strcmp')
       use iso_c_binding
       implicit none
       !pgi$ ignore_tkr a,b
       type(C_PTR), value :: a, b
     end function strcmp
  end interface
contains
  subroutine assignDevice(dev)
    use cudafor
    implicit none
    integer :: dev
    character (len=MPI_MAX_PROCESSOR_NAME), allocatable :: hosts(:)
    character (len=MPI_MAX_PROCESSOR_NAME) :: hostname
    integer :: namelength, color, i, j
    integer :: nProcs, myrank, newComm, newRank, ierr

    call MPI_COMM_SIZE(MPI_COMM_WORLD, nProcs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, myrank, ierr)

    ! allocate array of hostnames
    allocate(hosts(0:nProcs-1))
  
    ! Every process collects the hostname of all the nodes
    call MPI_GET_PROCESSOR_NAME(hostname, namelength, ierr)
    hosts(myrank)=hostname(1:namelength)

    do i=0,nProcs-1
       call MPI_BCAST(hosts(i),MPI_MAX_PROCESSOR_NAME,MPI_CHARACTER,i, &
            MPI_COMM_WORLD,ierr)
    end do
  
    ! sort the list of names
    call quicksort(hosts,nProcs,MPI_MAX_PROCESSOR_NAME,strcmp)

    ! assign the same color to the same node
    color=0
    do i=0,nProcs-1
       if (i > 0) then
          if ( lne(hosts(i-1),hosts(i)) ) color=color+1
       end if
       if ( leq(hostname,hosts(i)) ) exit
    end do
  
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,0,newComm,ierr)
    call MPI_COMM_RANK(newComm, newRank, ierr)

    dev = newRank
    ierr = cudaSetDevice(dev)
    
#if 0
    do i=0,nProcs-1
      if(myrank == i) then
          write(6,"(A8,I4,A8,A12,A12,I2)") "Rank: ",myrank,"Host: ",hostname(1:namelength),"Using GPU: ",dev
      endif
      do j=0,1000
        call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      end do
    end do
#endif
    deallocate(hosts)
  end subroutine assignDevice

  ! lexical .eq.
  function leq(s1, s2) result(res)
    implicit none
    character (len=*) :: s1, s2
    logical :: res    
    res = .false.
    if (lle(s1,s2) .and. lge(s1,s2)) res = .true.
  end function leq

  ! lexical .ne.
  function lne(s1, s2) result(res)
    implicit none
    character (len=*) :: s1, s2
    logical :: res    
    res = .not. leq(s1, s2)
  end function lne
end module mpiDeviceUtil
#endif

! ----
! nvtx
! ----

module nvtx
  use iso_c_binding
#ifdef USE_CUDA
  use cudafor
#endif
  implicit none
#ifdef USE_NVTX
  integer,private :: col(7) = [ Z'0000ff00', Z'000000ff', Z'00ffff00',Z'00ff00ff',Z'0000ffff', &
                                Z'00ff0000', Z'00ffffff']
  character(len=256),private :: tempName
!  logical, save :: use_nvtx=.false.
  type, bind(C):: nvtxEventAttributes
     integer(C_INT16_T):: version=1
     integer(C_INT16_T):: size=48 !
     integer(C_INT):: category=0
     integer(C_INT):: colorType=1 ! NVTX_COLOR_ARGB = 1
     integer(C_INT):: color
     integer(C_INT):: payloadType=0 ! NVTX_PAYLOAD_UNKNOWN = 0
     integer(C_INT):: reserved0
     integer(C_INT64_T):: payload   ! union uint,int,double
     integer(C_INT):: messageType=1  ! NVTX_MESSAGE_TYPE_ASCII     = 1 
     type(C_PTR):: message  ! ascii char
  end type nvtxEventAttributes

  interface nvtxRangePush
     ! push range with custom label and standard color
     subroutine nvtxRangePushA(name) bind(C, name='nvtxRangePushA')
       use iso_c_binding
       character(kind=C_CHAR,len=*) :: name
     end subroutine nvtxRangePushA

     ! push range with custom label and custom color
     subroutine nvtxRangePushEx(event) bind(C, name='nvtxRangePushEx')
       use iso_c_binding
       import:: nvtxEventAttributes
       type(nvtxEventAttributes):: event
     end subroutine nvtxRangePushEx
  end interface nvtxRangePush

  interface nvtxRangePop
     subroutine nvtxRangePop() bind(C, name='nvtxRangePop')
     end subroutine nvtxRangePop
  end interface nvtxRangePop
#endif

contains

  subroutine nvtxStartRange(name,id)
    character(kind=c_char,len=*) :: name
    integer, optional:: id
#ifdef USE_NVTX
    type(nvtxEventAttributes):: event
#ifdef USE_CUDA
    integer :: istat
    istat = cudaDeviceSynchronize()
#endif

    tempName=trim(name)//c_null_char

    if ( .not. present(id)) then
       call nvtxRangePush(tempName)
    else
       event%color=col(mod(id,7)+1)
       event%message=c_loc(tempName)
       call nvtxRangePushEx(event)
    end if
#endif
  end subroutine nvtxStartRange

  subroutine nvtxStartRangeAsync(name,id)
    character(kind=c_char,len=*) :: name
    integer, optional:: id
#ifdef USE_NVTX
    type(nvtxEventAttributes):: event

    tempName=trim(name)//c_null_char

    if ( .not. present(id)) then
       call nvtxRangePush(tempName)
    else
       event%color=col(mod(id,7)+1)
       event%message=c_loc(tempName)
       call nvtxRangePushEx(event)
    end if
#endif
  end subroutine nvtxStartRangeAsync


  subroutine nvtxEndRange
#ifdef USE_NVTX
#ifdef USE_CUDA
    integer :: istat
    istat = cudaDeviceSynchronize()
#endif
    call nvtxRangePop
#endif
  end subroutine nvtxEndRange

  subroutine nvtxEndRangeAsync
#ifdef USE_NVTX
    call nvtxRangePop
#endif
  end subroutine nvtxEndRangeAsync

end module nvtx

#ifdef USE_CUDA

module ep_debug
#ifdef USE_CUDA
  use cudafor
#endif
  implicit none
#if defined(__MPI)
  INCLUDE 'mpif.h'
#endif
  INTEGER, PARAMETER :: DP = selected_real_kind(14,200)
  integer, save :: counter=0,counter2=0
  interface compare
     module procedure compare_cpu_gpu_complex_1d
     module procedure compare_cpu_gpu_complex_2d
     module procedure compare_cpu_gpu_complex_3d
     module procedure compare_cpu_gpu_real_1d
     module procedure compare_cpu_gpu_real_2d
!     module procedure compare_cpu_1d
!     module procedure compare_cpu_3d
!     module procedure compare_cpu_complex_3d
!#ifdef USE_CUDA
!     module procedure compare_gpu_1d
!     module procedure compare_gpu_3d
!     module procedure compare_gpu_complex_1d
!     module procedure compare_gpu_complex_3d
!#endif
  end interface compare

  interface write_compare
    module procedure write_complex_1d
  end interface write_compare

  interface read_compare
    module procedure read_complex_1d
    module procedure read_complex_2d
    module procedure read_complex_3d
    module procedure read_real_1d
    module procedure read_real_2d
    module procedure read_real_3d
  end interface read_compare

contains


  subroutine write_complex_1d(A,filename)
#ifdef USE_NVTX
  USE nvtx
#endif
    implicit none
    complex(DP), dimension(:) :: A,A_d
    character (len=*), OPTIONAL, intent(IN) :: filename
    character (len=4) :: itcount
    character (len=1) :: proc
    integer :: nrank, ierror, nproc
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
    !print *,"rank,nproc,ierror: ",nrank,nproc,ierror
    counter = counter + 1
    print *,"!!!!!  counter = ",counter," !!!!"
    write(itcount,'(i4)') counter
    write(proc,'(i1)') nrank
    open(unit=10, status='replace',file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)), form='unformatted')
    write(10) A
    close(10)
  end subroutine write_complex_1d

#if 0
  subroutine read_complex_1d(A_d,filename)
#ifdef USE_NVTX
  USE nvtx
#endif
    implicit none
    complex(DP), dimension(:) :: A_d
    attributes( device ) :: A_d
    character (len=*), OPTIONAL, intent(IN) :: filename
    !real(fp_kind), dimension(:), allocatable :: A_d
    complex(DP), pinned, dimension(:), allocatable :: A_h, A
    real(DP) :: maxerr,perr,l2normerr,norm,buf, rmaxerr,l2normerr_loc
    complex(DP) :: rAimax,rA_himax
    integer :: i,imax,rimax,myid,cntr
    integer :: mpistat(MPI_STATUS_SIZE)
    character (len=4) :: itcount
    character (len=1) :: proc
    integer :: nrank, ierror, nproc
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
    !print *,"rank,nproc,ierror: ",nrank,nproc,ierror
    counter = counter + 1
    write(itcount,'(i4)') counter-7
    write(proc,'(i1)') nrank
#ifdef USE_NVTX
  CALL nvtxStartRange("COMP",0)
#endif
    allocate(A, source=A_d)
    allocate(A_h, source=A_d)
    A=0.0d0
    open(unit=10, status='old', file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)),form='unformatted')
    read(10) A
    close(10)
    l2normerr = 0.d0
    norm = 0.d0
    maxerr = 0.d0
    imax=1
    cntr = 0
        do i=lbound(A,1),ubound(A,1)
          if(abs(A(i)) >= 1e-20) then
            perr = abs(A(i) - A_h(i))/abs(A(i))*100.d0
            norm = norm + abs(A(i)*A(i));
            l2normerr = l2normerr + abs((A(i) - A_h(i))*(A(i) - A_h(i)))
            if((perr > 0.001d0 .or. cntr<20)  .and. nrank==0) then
cntr = cntr + 1
        write(*,"(A10,I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,A6,2X,E20.14,1X,E20.14,2X,A6,2X,E20.14,1X,E20.14,1X,E20.14)") &
        filename,nrank,"l2norm error",l2normerr/norm,"error",perr,"% at ",i,"cpu=",REAL(A(i)),AIMAG(A(i)),"gpu=",REAL(A_h(i)),AIMAG(A_h(i))

              call flush(6)
if(cntr > 100) stop
            endif
          else
            perr = 0.d0
          endif
          if(perr>maxerr .and. abs(a(i))/=0.0d0 .and. abs(a_h(i))/=0.0d0) then
            maxerr = perr
            imax = i
          endif
        enddo

 call MPI_ALLREDUCE(norm,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
 norm = buf
 l2normerr_loc = l2normerr
 call MPI_ALLREDUCE(l2normerr,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
 l2normerr = buf

  norm = sqrt(norm)
  l2normerr = sqrt(l2normerr)

  if(l2normerr /= 0.d0) then
    l2normerr = l2normerr/norm
  do myid = 0,nproc-1

    if(nrank == 0) then

      if(myid == 0) then
        rmaxerr = maxerr
        rimax = imax
        rAimax = A(imax)
        rA_himax = A_h(imax)
      else
        !print *,"rank: ",nrank," recving",ierror
        call MPI_RECV(rmaxerr  ,1,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rimax    ,1,MPI_INTEGER         ,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rAimax   ,2,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rA_himax ,2,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
      endif
      IF( PRESENT( filename ) ) THEN
        write(*,"(A10,I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,A6,2X,E20.14,1X,E20.14,2X,A6,2X,E20.14,1X,E20.14,1X,E20.14)") &
        filename,myid,"l2norm error",l2normerr,"max error",rmaxerr,"% at ",rimax,"cpu=",REAL(rAimax),AIMAG(rAimax),"gpu=",REAL(rA_himax),AIMAG(rA_himax),l2normerr_loc
      ELSE
        write(*,"(I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,A6,2X,E20.14,1X,E20.14,2X,A6,2X,E20.14,1X,E20.14)") &
        myid,"l2norm error",l2normerr,"max error",rmaxerr,"% at ",rimax,"cpu=",REAL(rAimax),AIMAG(rAimax),"gpu=",REAL(rA_himax),AIMAG(rA_himax)
      END IF

    elseif(nrank == myid .and. myid /= 0) then

      !print *,"rank: ",nrank," sending",ierror
      call MPI_SEND(maxerr   ,1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(imax     ,1,MPI_INTEGER         ,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(A(imax)  ,2,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(A_h(imax),2,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)

     end if

  end do

if(nrank==0) then
        imax = 0
        do i=lbound(A,1),ubound(A,1)
            perr = abs(A(i) - A_h(i))/abs(A(i))*100.d0
            if(perr>0.001d0 .and. abs(a(i))/=0.0d0 .and. abs(a_h(i))/=0.0d0) then
              if(imax<100) then
                write(*,"(A12,ES10.3,A6,I8,A6,2X,E20.14,1X,E20.14,2X,A6,2X,E20.14,1X,E20.14,A6,E20.14)") &
                "error",perr,"% at ",i,"cpu=",REAL(A(i)),AIMAG(A(i)),"gpu=",REAL(A_h(i)),AIMAG(A_h(i)),"ratio ",REAL(A_h(i))/REAL(A(i))
              end if
              imax = imax + 1
            endif
        enddo
        print *,"NUMBER ERRORS = ",imax,"of",ubound(A,1)-lbound(A,1)+1
endif

  call flush(6)
  call MPI_BARRIER( MPI_COMM_WORLD, ierror)

   else

      IF( PRESENT( filename ) ) THEN
       write(*,"(A10,I4,A16)") filename, nrank, "EXACT MATCH"
      ELSE
       write(*,"(I4,A16)") nrank, "EXACT MATCH"
      ENDIF
   endif

  if(l2normerr >= 0.00001d0 ) then
    STOP
    !A_d = A
  endif
  deallocate(A_h)
#ifdef USE_NVTX
  CALL nvtxEndRange
#endif
  end subroutine read_complex_1d
#endif

  subroutine compare_cpu_gpu_real_1d(A,A_d,filename)
#ifdef USE_NVTX
  USE nvtx
#endif
    implicit none
    real(DP), dimension(:) :: A,A_d
    attributes( device ) :: A_d
    character (len=*), OPTIONAL, intent(IN) :: filename
    !real(fp_kind), dimension(:), allocatable :: A_d
    real(DP), pinned, dimension(:), allocatable :: A_h
    real(DP) :: maxerr,perr,l2normerr,norm,buf, rmaxerr,l2normerr_loc
    real(DP) :: rAimax,rA_himax
    integer :: i,imax,rimax,myid,cntr
    integer :: mpistat(MPI_STATUS_SIZE)
    !character (len=4) :: itcount
    !character (len=1) :: proc
    integer :: nrank, ierror, nproc
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
    !print *,"rank,nproc,ierror: ",nrank,nproc,ierror
    !counter = counter + 1
    !write(itcount,'(i4)') counter
    !write(proc,'(i1)') nrank
#ifdef USE_NVTX
  CALL nvtxStartRange("COMP",0)
#endif
    allocate(A_h, source=A_d)
    !A_h=0.0d0
    !open(unit=10, status='old', file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)),form='unformatted')
    !read(10) A_h
    !close(10)
    l2normerr = 0.d0
    norm = 0.d0
    maxerr = 0.d0
    imax=1
    cntr = 0
        do i=lbound(A,1),ubound(A,1)
          if(abs(A(i)) >= 1e-10) then
            perr = abs(A(i) - A_h(i))/abs(A(i))*100.d0
            norm = norm + abs(A(i)*A(i));
            l2normerr = l2normerr + abs((A(i) - A_h(i))*(A(i) - A_h(i)))
            if(perr>0.01 .and. nrank==0) then !if(l2normerr/norm >= 0.1 .and. nrank==0) then
cntr = cntr + 1
        write(*,"(A10,I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,A6,2X,E20.14,2X,A6,2X,E20.14,1X,E20.14)") &
        filename,nrank,"l2norm error",l2normerr/norm,"error",perr,"% at ",i,"cpu=",A(i),"gpu=",A_h(i)

              call flush(6)
if(cntr > 800) stop
            endif
          else
            perr = 0.d0
          endif
          if(perr>maxerr .and. abs(a(i))/=0.0d0 .and. abs(a_h(i))/=0.0d0) then
            maxerr = perr
            imax = i
          endif
        enddo

 call MPI_ALLREDUCE(norm,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
 norm = buf
 l2normerr_loc = l2normerr
 call MPI_ALLREDUCE(l2normerr,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
 l2normerr = buf

  norm = sqrt(norm)
  l2normerr = sqrt(l2normerr)

  if(l2normerr /= 0.d0) then
    l2normerr = l2normerr/norm
  do myid = 0,nproc-1

    if(nrank == 0) then

      if(myid == 0) then
        rmaxerr = maxerr
        rimax = imax
        rAimax = A(imax)
        rA_himax = A_h(imax)
      else
        !print *,"rank: ",nrank," recving",ierror
        call MPI_RECV(rmaxerr  ,1,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rimax    ,1,MPI_INTEGER         ,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rAimax   ,1,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rA_himax ,1,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
      endif
      IF( PRESENT( filename ) ) THEN
        write(*,"(A10,I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,A6,2X,E20.14,2X,A6,2X,E20.14,1X,E20.14)") &
        filename,myid,"l2norm error",l2normerr,"max error",rmaxerr,"% at ",rimax,"cpu=",rAimax,"gpu=",rA_himax,l2normerr_loc
      ELSE
        write(*,"(I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,A6,2X,E20.14,2X,A6,2X,E20.14)") &
        myid,"l2norm error",l2normerr,"max error",rmaxerr,"% at ",rimax,"cpu=",rAimax,"gpu=",rA_himax
      END IF

    elseif(nrank == myid .and. myid /= 0) then

      !print *,"rank: ",nrank," sending",ierror
      call MPI_SEND(maxerr   ,1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(imax     ,1,MPI_INTEGER         ,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(A(imax)  ,1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(A_h(imax),1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)

     end if

  end do

  call flush(6)
  call MPI_BARRIER( MPI_COMM_WORLD, ierror)

   else

      IF( PRESENT( filename ) ) THEN
       write(*,"(A10,I4,A16)") filename, nrank, "EXACT MATCH"
      ELSE
       write(*,"(I4,A16)") nrank, "EXACT MATCH"
      ENDIF
   endif

  if(l2normerr >= 0.00001d0 ) then
    STOP
    !A_d = A
  endif
  deallocate(A_h)
#ifdef USE_NVTX
  CALL nvtxEndRange
#endif
  end subroutine compare_cpu_gpu_real_1d

  subroutine compare_cpu_gpu_complex_1d(A,A_d,filename)
#ifdef USE_NVTX
  USE nvtx
#endif
    implicit none
    complex(DP), dimension(:) :: A,A_d
    attributes( device ) :: A_d
    character (len=*), OPTIONAL, intent(IN) :: filename
    !real(fp_kind), dimension(:), allocatable :: A_d
    complex(DP), pinned, dimension(:), allocatable :: A_h
    real(DP) :: maxerr,perr,l2normerr,norm,buf, rmaxerr,l2normerr_loc
    complex(DP) :: rAimax,rA_himax
    integer :: i,imax,rimax,myid,cntr
    integer :: mpistat(MPI_STATUS_SIZE)
    !character (len=4) :: itcount
    !character (len=1) :: proc
    integer :: nrank, ierror, nproc
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
    !print *,"rank,nproc,ierror: ",nrank,nproc,ierror
    !counter = counter + 1
    !write(itcount,'(i4)') counter
    !write(proc,'(i1)') nrank
#ifdef USE_NVTX
  CALL nvtxStartRange("COMP",0)
#endif
    allocate(A_h, source=A_d)
    !A_h=0.0d0
    !open(unit=10, status='old', file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)),form='unformatted')
    !read(10) A_h
    !close(10)
    l2normerr = 0.d0
    norm = 0.d0
    maxerr = 0.d0
    imax=1
    cntr = 0
        do i=lbound(A,1),ubound(A,1)
          if(abs(A(i)) >= 1e-10) then
            perr = abs(A(i) - A_h(i))/abs(A(i))*100.d0
            norm = norm + abs(A(i)*A(i));
            l2normerr = l2normerr + abs((A(i) - A_h(i))*(A(i) - A_h(i)))
            if((perr > 0.001d0 .or. cntr>20)  .and. nrank==0) then
cntr = cntr + 1
        write(*,"(A10,I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,A6,2X,E20.14,1X,E20.14,2X,A6,2X,E20.14,1X,E20.14,1X,E20.14)") &
        filename,nrank,"l2norm error",l2normerr/norm,"error",perr,"% at ",i,"cpu=",REAL(A(i)),AIMAG(A(i)),"gpu=",REAL(A_h(i)),AIMAG(A_h(i))

              call flush(6)
if(cntr > 100) stop
            endif
          else
            perr = 0.d0
          endif
          if(perr>maxerr .and. abs(a(i))/=0.0d0 .and. abs(a_h(i))/=0.0d0) then
            maxerr = perr
            imax = i
          endif
        enddo

 call MPI_ALLREDUCE(norm,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
 norm = buf
 l2normerr_loc = l2normerr
 call MPI_ALLREDUCE(l2normerr,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
 l2normerr = buf

  norm = sqrt(norm)
  l2normerr = sqrt(l2normerr)

  if(l2normerr /= 0.d0) then
    l2normerr = l2normerr/norm
  do myid = 0,nproc-1

    if(nrank == 0) then

      if(myid == 0) then
        rmaxerr = maxerr
        rimax = imax
        rAimax = A(imax)
        rA_himax = A_h(imax)
      else
        !print *,"rank: ",nrank," recving",ierror
        call MPI_RECV(rmaxerr  ,1,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rimax    ,1,MPI_INTEGER         ,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rAimax   ,2,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rA_himax ,2,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
      endif
      IF( PRESENT( filename ) ) THEN
        write(*,"(A10,I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,A6,2X,E20.14,1X,E20.14,2X,A6,2X,E20.14,1X,E20.14,1X,E20.14)") &
        filename,myid,"l2norm error",l2normerr,"max error",rmaxerr,"% at ",rimax,"cpu=",REAL(rAimax),AIMAG(rAimax),"gpu=",REAL(rA_himax),AIMAG(rA_himax),l2normerr_loc
      ELSE
        write(*,"(I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,A6,2X,E20.14,1X,E20.14,2X,A6,2X,E20.14,1X,E20.14)") &
        myid,"l2norm error",l2normerr,"max error",rmaxerr,"% at ",rimax,"cpu=",REAL(rAimax),AIMAG(rAimax),"gpu=",REAL(rA_himax),AIMAG(rA_himax)
      END IF

    elseif(nrank == myid .and. myid /= 0) then

      !print *,"rank: ",nrank," sending",ierror
      call MPI_SEND(maxerr   ,1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(imax     ,1,MPI_INTEGER         ,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(A(imax)  ,2,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(A_h(imax),2,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)

     end if

  end do
#if 0
if(nrank==0) then
        imax = 0
        do i=lbound(A,1),ubound(A,1)
            perr = abs(A(i) - A_h(i))/abs(A(i))*100.d0
            if(perr>0.001d0 .and. abs(a(i))/=0.0d0 .and. abs(a_h(i))/=0.0d0) then
              if(imax<100) then
                write(*,"(A12,ES10.3,A6,I8,A6,2X,E20.14,1X,E20.14,2X,A6,2X,E20.14,1X,E20.14,A6,E20.14)") &
                "error",perr,"% at ",i,"cpu=",REAL(A(i)),AIMAG(A(i)),"gpu=",REAL(A_h(i)),AIMAG(A_h(i)),"ratio ",REAL(A_h(i))/REAL(A(i))
              end if
              imax = imax + 1
            endif
        enddo
        print *,"NUMBER ERRORS = ",imax,"of",ubound(A,1)-lbound(A,1)+1
endif
#endif
  call flush(6)
  call MPI_BARRIER( MPI_COMM_WORLD, ierror)

   else

      IF( PRESENT( filename ) ) THEN
       write(*,"(A10,I4,A16)") filename, nrank, "EXACT MATCH"
      ELSE
       write(*,"(I4,A16)") nrank, "EXACT MATCH"
      ENDIF
   endif

  if(l2normerr >= 0.00001d0 ) then
    STOP
    !A_d = A
  endif
  deallocate(A_h)
#ifdef USE_NVTX
  CALL nvtxEndRange
#endif
  end subroutine compare_cpu_gpu_complex_1d


  subroutine compare_cpu_gpu_real_2d(A,A_d,filename)
#ifdef USE_NVTX
  USE nvtx
#endif
    implicit none
    REAL(DP), dimension(:,:) :: A,A_d
    attributes( device ) :: A_d
    character (len=*), OPTIONAL, intent(IN) :: filename
    !real(fp_kind), dimension(:), allocatable :: A_d
    REAL(DP), pinned, dimension(:,:), allocatable :: A_h
    real(DP) :: maxerr,perr,l2normerr,norm,buf, rmaxerr
    REAL(DP) :: rAimax,rA_himax
    integer :: i,imax,rimax,myid,cntr
    integer :: j,jmax,rjmax
    integer :: mpistat(MPI_STATUS_SIZE)
    !character (len=4) :: itcount
    !character (len=1) :: proc
    integer :: nrank, ierror, nproc
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
    !print *,"rank,nproc,ierror: ",nrank,nproc,ierror
    !counter = counter + 1
    !write(itcount,'(i4)') counter
    !write(proc,'(i1)') nrank
#ifdef USE_NVTX
  CALL nvtxStartRange("COMP",0)
#endif
    allocate(A_h, source=A_d)
    !A_h=0.0d0
    !open(unit=10, status='old', file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)),form='unformatted')
    !read(10) A_h
    !close(10)
    l2normerr = 0.d0
    norm = 0.d0
    maxerr = 0.d0
    imax=1; jmax=1
    cntr = 0
      do j=lbound(A,2),ubound(A,2)
        do i=lbound(A,1),ubound(A,1)
          if(abs(A(i,j)) >= 1e-10) then
            perr = abs(A(i,j) - A_h(i,j))/abs(A(i,j))*100.d0
            norm = norm + abs(A(i,j)*A(i,j));
            l2normerr = l2normerr + abs((A(i,j) - A_h(i,j))*(A(i,j) - A_h(i,j)))
            if(perr >= 0.01 .and. nrank==0) then
cntr = cntr + 1
        write(*,"(A10,I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,I8,A6,2X,E20.14,2X,A6,2X,E20.14)") &
        filename,nrank,"l2norm error",l2normerr/norm,"error",perr,"% at ",i,j,"cpu=",A(i,j),"gpu=",A_h(i,j)
              call flush(6)
if(cntr > 100) stop
             if(i>2) exit
            endif

          else
            perr = 0.d0
          endif
          if(perr>maxerr .and. abs(a(i,j))/=0.0d0 .and. abs(a_h(i,j))/=0.0d0) then
            maxerr = perr
            imax = i
            jmax = j
          endif
        enddo
      enddo
 call MPI_ALLREDUCE(norm,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
 norm = buf
 call MPI_ALLREDUCE(l2normerr,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
 l2normerr = buf

  norm = sqrt(norm)
  l2normerr = sqrt(l2normerr)

  if(l2normerr /= 0.d0) then
    l2normerr = l2normerr/norm

  do myid = 0,nproc-1

    if(nrank == 0) then

      if(myid == 0) then
        rmaxerr = maxerr
        rimax = imax
        rjmax = jmax
        rAimax = A(imax,jmax)
        rA_himax = A_h(imax,jmax)
      else
        !print *,"rank: ",nrank," recving",ierror
        call MPI_RECV(rmaxerr  ,1,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rimax    ,1,MPI_INTEGER         ,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rjmax    ,1,MPI_INTEGER         ,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rAimax   ,1,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rA_himax ,1,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
      endif
      IF( PRESENT( filename ) ) THEN
        write(*,"(A10,I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,I8,A6,2X,E20.14,2X,A6,2X,E20.14)") &
        filename,myid,"l2norm error",l2normerr,"max error",rmaxerr,"% at ",rimax,rjmax,"cpu=",rAimax,"gpu=",rA_himax
      ELSE
        write(*,"(    I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,I8,A6,2X,E20.14,2X,A6,2X,E20.14)") &
                 myid,"l2norm error",l2normerr,"max error",rmaxerr,"% at ",rimax,rjmax,"cpu=",rAimax,"gpu=",rA_himax
      END IF

    elseif(nrank == myid .and. myid /= 0) then

      !print *,"rank: ",nrank," sending"
      call MPI_SEND(maxerr        ,1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(imax          ,1,MPI_INTEGER         ,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(jmax          ,1,MPI_INTEGER         ,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(A(imax,jmax)  ,1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(A_h(imax,jmax),1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)

    end if

  end do

   else

      IF( PRESENT( filename ) ) THEN
       write(*,"(A10,I4,A16)") filename, nrank, "EXACT MATCH"
      ELSE
       write(*,"(I4,A16)") nrank, "EXACT MATCH"
      ENDIF

   endif

  deallocate(A_h)

  if(l2normerr >= 0.00001d0 ) STOP

#ifdef USE_NVTX
  CALL nvtxEndRange
#endif
  end subroutine compare_cpu_gpu_real_2d

  subroutine compare_cpu_gpu_complex_2d(A,A_d,filename)
#ifdef USE_NVTX
  USE nvtx
#endif
    implicit none
    complex(DP), dimension(:,:) :: A,A_d
    attributes( device ) :: A_d
    character (len=*), OPTIONAL, intent(IN) :: filename
    !real(fp_kind), dimension(:), allocatable :: A_d
    complex(DP), pinned, dimension(:,:), allocatable :: A_h
    real(DP) :: maxerr,perr,nperr,l2normerr,norm,buf, rmaxerr
    complex(DP) :: rAimax,rA_himax
    integer :: i,imax,rimax,myid,cntr
    integer :: j,jmax,rjmax
    integer :: mpistat(MPI_STATUS_SIZE)
    !character (len=4) :: itcount
    !character (len=1) :: proc
    integer :: nrank, ierror, nproc
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
    !print *,"rank,nproc,ierror: ",nrank,nproc,ierror
    !counter = counter + 1
    !write(itcount,'(i4)') counter
    !write(proc,'(i1)') nrank
#ifdef USE_NVTX
  CALL nvtxStartRange("COMP",0)
#endif
    allocate(A_h, source=A_d)
    !A_h=0.0d0
    !open(unit=10, status='old', file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)),form='unformatted')
    !read(10) A_h
    !close(10)
    l2normerr = 0.d0
    norm = 0.d0
    maxerr = 0.d0
    imax=1; jmax=1
    cntr = 0
!    print *,"check 2d bounds :"
!    print *,lbound(A,2),ubound(A,2)
!    print *,lbound(A,1),ubound(A,1)

      do j=lbound(A,2),ubound(A,2)
        do i=lbound(A,1),ubound(A,1)
          if(1) then !abs(A(i,j)) >= 1e-32) then
            if(abs(A(i,j)) >= 1e-10) then
              perr = abs(A(i,j) - A_h(i,j))/abs(A(i,j))*100.d0
             nperr = abs(A(i,j) + A_h(i,j))/abs(A(i,j))*100.d0
              if(nperr < perr) then 
                perr = nperr
                A_h(i,j) = -A_h(i,j)
              endif
            else
              perr = 0.0
            endif
            norm = norm + abs(A(i,j)*A(i,j));
            l2normerr = l2normerr + abs((A(i,j) - A_h(i,j))*(A(i,j) - A_h(i,j)))
            if(perr >= 0.01 .and. nrank==0) then
cntr = cntr + 1
        write(*,"(A10,I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,I8,A6,2X,E20.14,1X,E20.14,2X,A6,2X,E20.14,1X,E20.14,1X,E20.14)") &
        filename,nrank,"l2norm error",l2normerr/norm,"error",perr,"% at ",i,j,"cpu=",REAL(A(i,j)),AIMAG(A(i,j)),"gpu=",REAL(A_h(i,j)),AIMAG(A_h(i,j))
              call flush(6)
if(cntr > 100) stop
             if(i>16) exit
            endif

          else
            perr = 0.d0
          endif
          if(perr>maxerr .and. abs(a(i,j))/=0.0d0 .and. abs(a_h(i,j))/=0.0d0) then
            maxerr = perr
            imax = i
            jmax = j
          endif
        enddo
      enddo
 call MPI_ALLREDUCE(norm,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
 norm = buf
 call MPI_ALLREDUCE(l2normerr,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
 l2normerr = buf

  norm = sqrt(norm)
  l2normerr = sqrt(l2normerr)

  if(l2normerr /= 0.d0) then
    l2normerr = l2normerr/norm

  do myid = 0,nproc-1

    if(nrank == 0) then

      if(myid == 0) then
        rmaxerr = maxerr
        rimax = imax
        rjmax = jmax
        rAimax = A(imax,jmax)
        rA_himax = A_h(imax,jmax)
      else
        !print *,"rank: ",nrank," recving",ierror
        call MPI_RECV(rmaxerr  ,1,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rimax    ,1,MPI_INTEGER         ,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rjmax    ,1,MPI_INTEGER         ,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rAimax   ,2,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rA_himax ,2,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
      endif
      IF( PRESENT( filename ) ) THEN
        write(*,"(A10,I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,I8,A6,2X,E20.14,1X,E20.14,2X,A6,2X,E20.14,1X,E20.14)") &
        filename,myid,"l2norm error",l2normerr,"max error",rmaxerr,"% at ",rimax,rjmax,"cpu=",REAL(rAimax),AIMAG(rAimax),"gpu=",REAL(rA_himax),AIMAG(rA_himax)
      ELSE
        write(*,"(    I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,I8,A6,2X,E20.14,1X,E20.14,2X,A6,2X,E20.14,1X,E20.14)") &
                 myid,"l2norm error",l2normerr,"max error",rmaxerr,"% at ",rimax,rjmax,"cpu=",REAL(rAimax),AIMAG(rAimax),"gpu=",REAL(rA_himax),AIMAG(rA_himax)
      END IF

    elseif(nrank == myid .and. myid /= 0) then

      !print *,"rank: ",nrank," sending"
      call MPI_SEND(maxerr        ,1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(imax          ,1,MPI_INTEGER         ,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(jmax          ,1,MPI_INTEGER         ,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(A(imax,jmax)  ,2,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(A_h(imax,jmax),2,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)

    end if

  end do

   else

      IF( PRESENT( filename ) ) THEN
       write(*,"(A10,I4,A16)") filename, nrank, "EXACT MATCH"
      ELSE
       write(*,"(I4,A16)") nrank, "EXACT MATCH"
      ENDIF

   endif

  deallocate(A_h)

  if(l2normerr >= 0.00001d0 ) STOP

#ifdef USE_NVTX
  CALL nvtxEndRange
#endif
  end subroutine compare_cpu_gpu_complex_2d


  subroutine compare_cpu_gpu_complex_3d(A,A_d,filename)
#ifdef USE_NVTX
  USE nvtx
#endif
    implicit none
    complex(DP), dimension(:,:,:) :: A,A_d
    attributes( device ) :: A_d
    character (len=*), OPTIONAL, intent(IN) :: filename
    !real(fp_kind), dimension(:), allocatable :: A_d
    complex(DP), pinned, dimension(:,:,:), allocatable :: A_h
    real(DP) :: maxerr,perr,nperr,l2normerr,norm,buf, rmaxerr
    complex(DP) :: rAimax,rA_himax
    integer :: i,imax,rimax,myid,cntr
    integer :: j,jmax,rjmax
    integer :: k,kmax,rkmax
    integer :: mpistat(MPI_STATUS_SIZE)
    !character (len=4) :: itcount
    !character (len=1) :: proc
    integer :: nrank, ierror, nproc
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
    !print *,"rank,nproc,ierror: ",nrank,nproc,ierror
    !counter = counter + 1
    !write(itcount,'(i4)') counter
    !write(proc,'(i1)') nrank
#ifdef USE_NVTX
  CALL nvtxStartRange("COMP",0)
#endif
    allocate(A_h, source=A_d)
    !A_h=0.0d0
    !open(unit=10, status='old', file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)),form='unformatted')
    !read(10) A_h
    !close(10)
    l2normerr = 0.d0
    norm = 0.d0
    maxerr = 0.d0
    imax=1; jmax=1; kmax=1
    do k=lbound(A,3),ubound(A,3)
      do j=lbound(A,2),ubound(A,2)
        do i=lbound(A,1),ubound(A,1)
          if(abs(A(i,j,k)) >= 1e-10) then
            perr = abs(A(i,j,k) - A_h(i,j,k))/abs(A(i,j,k))*100.d0
           nperr = abs(A(i,j,k) + A_h(i,j,k))/abs(A(i,j,k))*100.d0
            if(nperr < perr) then 
               perr = nperr
               A_h(i,j,k) = -A_h(i,j,k)
            endif

            norm = norm + abs(A(i,j,k)*A(i,j,k));
            l2normerr = l2normerr + abs((A(i,j,k) - A_h(i,j,k))*(A(i,j,k) - A_h(i,j,k)))
          else
            perr = 0.d0
          endif
          if(perr>maxerr .and. abs(a(i,j,k))/=0.0d0 .and. abs(a_h(i,j,k))/=0.0d0) then
            maxerr = perr
            imax = i
            jmax = j
            kmax = k
          endif
        enddo
      enddo
    enddo
 call MPI_ALLREDUCE(norm,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
 norm = buf
 call MPI_ALLREDUCE(l2normerr,buf,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
 l2normerr = buf

  norm = sqrt(norm)
  l2normerr = sqrt(l2normerr)

  if(l2normerr /= 0.d0) then
    l2normerr = l2normerr/norm

  do myid = 0,nproc-1

    if(nrank == 0) then

      if(myid == 0) then
        rmaxerr = maxerr
        rimax = imax
        rjmax = jmax
        rkmax = kmax
        rAimax = A(imax,jmax,kmax)
        rA_himax = A_h(imax,jmax,kmax)
      else
        !print *,"rank: ",nrank," recving",ierror
        call MPI_RECV(rmaxerr  ,1,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rimax    ,1,MPI_INTEGER         ,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rjmax    ,1,MPI_INTEGER         ,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rkmax    ,1,MPI_INTEGER         ,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rAimax   ,2,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
        call MPI_RECV(rA_himax ,2,MPI_DOUBLE_PRECISION,myid,MPI_ANY_TAG,MPI_COMM_WORLD,mpistat,ierror)
      endif
      IF( PRESENT( filename ) ) THEN
        write(*,"(A10,I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,I8,I8,A6,2X,E20.14,1X,E20.14,2X,A6,2X,E20.14,1X,E20.14)") &
        filename,myid,"l2norm error",l2normerr,"max error",rmaxerr,"% at ",rimax,rjmax,rkmax,"cpu=",REAL(rAimax),AIMAG(rAimax),"gpu=",REAL(rA_himax),AIMAG(rA_himax)
      ELSE
        write(*,"(    I4,A16,2X,ES10.3,A12,ES10.3,A6,I8,I8,I8,A6,2X,E20.14,1X,E20.14,2X,A6,2X,E20.14,1X,E20.14)") &
                 myid,"l2norm error",l2normerr,"max error",rmaxerr,"% at ",rimax,rjmax,rkmax,"cpu=",REAL(rAimax),AIMAG(rAimax),"gpu=",REAL(rA_himax),AIMAG(rA_himax)
      END IF

    elseif(nrank == myid .and. myid /= 0) then

      !print *,"rank: ",nrank," sending"
      call MPI_SEND(maxerr        ,1,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(imax          ,1,MPI_INTEGER         ,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(jmax          ,1,MPI_INTEGER         ,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(kmax          ,1,MPI_INTEGER         ,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(A(imax,jmax,kmax)  ,2,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)
      call MPI_SEND(A_h(imax,jmax,kmax),2,MPI_DOUBLE_PRECISION,0,0,MPI_COMM_WORLD,ierror)

    end if

  end do

   else

      IF( PRESENT( filename ) ) THEN
       write(*,"(A10,I4,A16)") filename, nrank, "EXACT MATCH"
      ELSE
       write(*,"(I4,A16)") nrank, "EXACT MATCH"
      ENDIF

   endif

  deallocate(A_h)

  if(l2normerr >= 0.00001d0 ) STOP

#ifdef USE_NVTX
  CALL nvtxEndRange
#endif
  end subroutine compare_cpu_gpu_complex_3d

  subroutine read_complex_1d(A_d, filename)
    implicit none
    complex(DP), dimension(:) :: A_d
    attributes( device ) :: A_d
    character (len=*), OPTIONAL, intent(IN) :: filename
    complex(DP), pinned, dimension(:), allocatable :: A_h
    character (len=4) :: itcount
    character (len=1) :: proc
    integer :: nrank, ierror, nproc
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
    counter = counter + 1
    write(itcount,'(i4)') counter
    write(proc,'(i1)') nrank
    allocate(A_h, source=A_d)
    A_h=0.0d0
    open(unit=10, status='old', file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)),form='unformatted')
    read(10) A_h
    close(10)
    call compare( A_h, A_d, filename )
    deallocate( A_h )
  end subroutine read_complex_1d

  subroutine read_complex_2d(A_d, filename)
    implicit none
    complex(DP), dimension(:,:) :: A_d
    attributes( device ) :: A_d
    character (len=*), OPTIONAL, intent(IN) :: filename
    complex(DP), pinned, dimension(:,:), allocatable :: A_h
    character (len=4) :: itcount
    character (len=1) :: proc
    integer :: nrank, ierror, nproc
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
    counter = counter + 1
    write(itcount,'(i4)') counter
    write(proc,'(i1)') nrank
    allocate(A_h, source=A_d)
    A_h=0.0d0
    open(unit=10, status='old', file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)),form='unformatted')
    read(10) A_h
    close(10)
    call compare( A_h, A_d, filename )
    deallocate( A_h )
  end subroutine read_complex_2d

  subroutine read_complex_3d(A_d, filename)
    implicit none
    complex(DP), dimension(:,:,:) :: A_d
    attributes( device ) :: A_d
    character (len=*), OPTIONAL, intent(IN) :: filename
    complex(DP), pinned, dimension(:,:,:), allocatable :: A_h
    character (len=4) :: itcount
    character (len=1) :: proc
    integer :: nrank, ierror, nproc
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
    counter = counter + 1
    write(itcount,'(i4)') counter
    write(proc,'(i1)') nrank
    allocate(A_h, source=A_d)
    A_h=0.0d0
    open(unit=10, status='old', file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)),form='unformatted')
    read(10) A_h
    close(10)
    call compare( A_h, A_d, filename )
    deallocate( A_h )
  end subroutine read_complex_3d

  subroutine read_real_1d(A_d, filename)
    implicit none
    real(DP), dimension(:) :: A_d
    attributes( device ) :: A_d
    character (len=*), OPTIONAL, intent(IN) :: filename
    real(DP), pinned, dimension(:), allocatable :: A_h
    character (len=4) :: itcount
    character (len=1) :: proc
    integer :: nrank, ierror, nproc
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
    counter = counter + 1
    write(itcount,'(i4)') counter
    write(proc,'(i1)') nrank
    allocate(A_h, source=A_d)
    A_h=0.0d0
    open(unit=10, status='old', file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)),form='unformatted')
    read(10) A_h
    close(10)
    call compare( A_h, A_d, filename )
    deallocate( A_h )
  end subroutine read_real_1d

  subroutine read_real_2d(A_d, filename)
    implicit none
    real(DP), dimension(:,:) :: A_d
    attributes( device ) :: A_d
    character (len=*), OPTIONAL, intent(IN) :: filename
    real(DP), pinned, dimension(:,:), allocatable :: A_h
    character (len=4) :: itcount
    character (len=1) :: proc
    integer :: nrank, ierror, nproc
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
    counter = counter + 1
    write(itcount,'(i4)') counter
    write(proc,'(i1)') nrank
    allocate(A_h, source=A_d)
    A_h=0.0d0
    open(unit=10, status='old', file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)),form='unformatted')
    read(10) A_h
    close(10)
    call compare( A_h, A_d, filename )
    deallocate( A_h )
  end subroutine read_real_2d

  subroutine read_real_3d(A_d, filename)
    implicit none
    real(DP), dimension(:,:,:) :: A_d
    attributes( device ) :: A_d
    character (len=*), OPTIONAL, intent(IN) :: filename
    real(DP), pinned, dimension(:,:,:), allocatable :: A_h
    character (len=4) :: itcount
    character (len=1) :: proc
    integer :: nrank, ierror, nproc, i
    call MPI_BARRIER(MPI_COMM_WORLD,ierror)
    call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierror)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
    counter = counter + 1
    write(itcount,'(i4)') counter
    write(proc,'(i1)') nrank
    allocate(A_h, source=A_d)
    A_h=0.0d0
    open(unit=10, status='old', file=filename//TRIM(ADJUSTL(itcount))//".dbg."//TRIM(ADJUSTL(proc)),form='unformatted')
    read(10) A_h
    close(10)

    Do i=lbound(A_d,3), ubound(A_d,3)
      call compare( A_h(:,:,i), A_d(:,:,i), filename )
    end do

    deallocate( A_h )
  end subroutine read_real_3d

end module ep_debug

#endif
