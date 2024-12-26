!-------------------------------------------------------------
! Extracted from mergertree.f90
! Short program to read in progenitor binary MPI data dump.
!-------------------------------------------------------------



program read_progenitors

  use mpi
  implicit none
  ! include 'mpif.h'


  !--------------------------------------
  ! NEEDS TO BE SET MANUALLY
  !--------------------------------------

  integer :: ifout = 11 ! output_00010 is at ifout=11




  integer :: info, myid, ntasks
  integer, parameter :: dp = kind(1.d0), i8b=kind(1)
  logical :: verbose = .true., make_mock_galaxies=.true.





  !----------------------
  ! Mergertree general
  !----------------------

  logical :: make_mergertree = .true.     ! whether to make merger trees
  logical :: use_exclusive_mass = .true.  ! whether to use exclusive or inclusive halo mass definition for treemaking
  integer :: nmost_bound = 200            ! maximal number of most bound particles to track

  real(dp), allocatable, dimension(:,:) :: most_bound_energy      ! stores the energy of nmost_bound particles per peak
  integer,  allocatable, dimension(:,:) :: most_bound_pid         ! stores particle ID of nmost_bound particles per peak

  integer,  allocatable, dimension(:)   :: prog_id                ! global ID of progenitors
  integer,  allocatable, dimension(:)   :: prog_owner             ! CPU that owns the progenitor
  real(dp), allocatable, dimension(:)   :: prog_mass              ! list of progenitors masses
  integer,  allocatable, dimension(:)   :: tracers_all            ! list of progenitor tracers (global particle IDs) of all progs
  integer,  allocatable, dimension(:)   :: tracers_loc_pid        ! list of progenitor tracers (local particle IDs)
  integer,  allocatable, dimension(:)   :: tracer_loc_progids_all ! list of progenitor IDs for tracers (local prog ID) of all progs
  integer,  allocatable, dimension(:)   :: tracer_loc_progids     ! list of progenitor IDs for tracers (local prog ID)
                                                                  ! only on this CPU
  integer,  allocatable, dimension(:)   :: galaxy_tracers         ! list of active galaxy tracers 
                                                                  ! (the absolutely most bound  particle of progenitor) 
  integer,  allocatable, dimension(:)   :: main_prog              ! main progenitor of each descendant 
  integer,  allocatable, dimension(:)   :: main_desc              ! main descendant of each progenitor

  integer :: progenitorcount = 0          ! count the number of clumps that will be progenitors
  integer :: progenitorcount_written = 0  ! count the number of progenitors for output
  integer :: nprogs = 0                   ! number of progenitors read in/to work with for creating tree
  integer :: prog_free = 1                ! first free progenitor local index
  integer :: ntracers = 0                 ! number of tracers on this CPU



  real(dp), allocatable, dimension(:)   :: clmp_mass_exclusive  ! exclusive clump mass, containing only bound particles
  integer,  allocatable, dimension(:)   :: prog_outputnr        ! snapshot number of progenitor
  ! real(dp), allocatable, dimension(:,:) :: clmp_vel_exclusive
  integer  :: killed_tot, appended_tot ! count killed or appended clumps that were too small
  real(dp) :: partm_common






  !-------------------------------
  ! Progenitor-Descendant matrix
  !-------------------------------

  type prog_desc_mat


    ! dimension 1:nprogs
    integer, dimension(:), allocatable :: first ! first descendant for prog
    integer, dimension(:), allocatable :: cnt   ! number of descendants/progenitors

    ! dimension 1:10*npeaks_max
    integer, dimension(:), allocatable :: ntrace  ! number of tracers for this prog/desc pair
    integer, dimension(:), allocatable :: clmp_id ! descendant or progenitor ID. Desc: global; prog: local
    integer, dimension(:), allocatable :: next    ! next descendant for prog

    integer :: mat_free_ind = 1 ! first free index of matrix

  end type prog_desc_mat

  type(prog_desc_mat) :: p2d_links,d2p_links ! sparse matrix for progenitor/descendants linking





  !--------------------------
  ! Multi-shapshot matching
  !--------------------------


  integer, dimension(:), allocatable :: pmprogs         ! Past Merged Progenitors for multi-snapshot matching
  integer, dimension(:), allocatable :: pmprogs_galaxy  ! Past Merged Progenitors' galaxy particles
  integer, dimension(:), allocatable :: pmprogs_t       ! Time at which past progenitors have been merged (= ifout-1 at merging time)
  integer, dimension(:), allocatable :: pmprogs_owner   ! Current owner Merged Progenitors' galaxy particles
  real(dp),dimension(:), allocatable :: pmprogs_mass    ! Mass of past merged progenitors


  integer :: npastprogs = 0           ! number of past progenitors stored
  integer :: npastprogs_max           ! max number for array loops/allocation
  integer :: pmprog_free = 1          ! First free index in pmprogs* arrays
  integer :: max_past_snapshots = 0   ! maximal number of snapshots to store



  real(dp), allocatable, dimension(:) :: mpeak                      ! peak mass and expansion factor at time of mass peak
  real(dp), allocatable, dimension(:) :: prog_mpeak                 ! peak mass and expansion factor for progenitors
  real(dp), allocatable, dimension(:) :: pmprogs_mpeak              ! stellar mass of past merged progenitors

  integer, allocatable, dimension(:)  :: orphans_local_pid          ! local particle id of orphans
  integer, allocatable, dimension(:)  :: prog_galaxy_local_id       ! local particle id of progenitor galaxies











  call mpi_init(info)
  call mpi_comm_size(MPI_COMM_WORLD, ntasks, info)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, info)
  myid = myid + 1




  ! first allocate some global arrays for current peaks
  ! allocate(main_prog(1:npeaks_max))
  ! main_prog = 0
  ! allocate(prog_outputnr(1:npeaks_max))
  ! prog_outputnr = 0

  ! Read in progenitor files
  call read_progenitor_data()



  call mpi_finalize(info)



contains

!================================
subroutine make_merger_tree()
!================================

  !-----------------------------------------
  ! This subroutine is the main routine for
  ! creating merger trees that calls all
  ! others.
  !-----------------------------------------

  implicit none

  !----------------------------
  ! Create mergertrees
  !----------------------------

  if (myid==1)write(*,'(A31,x,I5)') " Calling merger tree for output", ifout

  ! create trees only if progenitors might exist
  if (ifout > 1) then 

    ! first allocate some global arrays for current peaks
    ! allocate(main_prog(1:npeaks_max))
    ! main_prog = 0
    ! allocate(prog_outputnr(1:npeaks_max))
    ! prog_outputnr = 0

    ! Read in progenitor files
    call read_progenitor_data()

    ! if (nprogs > 0) then
      ! ! Sort out progenitor data.
      ! call process_progenitor_data()
      !
      ! ! link progenitors and descendants
      ! call create_prog_desc_links()
      !
      ! ! create trees
      ! call make_trees()
    ! endif


  endif


  !--------------------
  ! Say good bye.
  !--------------------

  if(verbose) write(*,*) "Finished making merger tree."

end subroutine make_merger_tree






! !=====================================
! subroutine process_progenitor_data()
! !=====================================
!
!   !-----------------------------------------------------
!   ! This subroutine processes the progenitor data:
!   !   - It finds which tracer particles are currently
!   !     on this CPU
!   !   - It finds which CPU is owner of which progenitor,
!   !     and communicates that info across all CPUs.
!   !   - Counts how many tracer particles of any
!   !     progenitor are on this processor, create cleaned
!   !     up list
!   !------------------------------------------------------
!
!   implicit none
!
! #ifndef WITHOUTMPI
!   integer,      allocatable, dimension(:) :: local_owners_info
! #endif
!
!   integer(i8b), allocatable, dimension(:) :: idp_copy, galaxy_tracers_copy
!   integer,      allocatable, dimension(:) :: part_local_ind, sort_ind_past, dummy
!   real(dp),     allocatable, dimension(:) :: dummy_real
!   integer,      allocatable, dimension(:) :: tracers_local_pid_long, tracer_local_ids_long
!   integer :: itrace, ipart, igalaxy, ipastprog, i, iprog
!
!   if (verbose) write(*,*) " processing progenitor data."
!
!
!   !-------------------
!   ! Allocate stuff
!   !-------------------
!
!   ! Create copies of arrays you mustn't modify
!   ! allocate(idp_copy(1:npartmax))
!   ! idp_copy = idp
!
!   allocate(galaxy_tracers_copy(1:nprogs))
!   galaxy_tracers_copy = galaxy_tracers
!
!   ! allocate others
!   if (nprogs > npastprogs) then
!     allocate(dummy(1:nprogs))
!   else
!     allocate(dummy(1:npastprogs))
!   endif
!   dummy = 1
!   allocate(dummy_real(1:npastprogs))
!
!   ! allocate(part_local_ind(1:npartmax))
!   ! part_local_ind = [(i, i = 1, npartmax)]
!
!   allocate(sort_ind_past(1:npastprogs))
!   sort_ind_past = [(i, i=1, npastprogs)]
!
!   allocate(tracers_local_pid_long(1:nprogs*nmost_bound)) ! local particle   id for tracers (room enough for all tracers)
!   allocate(tracer_local_ids_long(1:nprogs*nmost_bound))  ! local progenitor id for tracers (room enough for all tracers)
!
!   if (make_mock_galaxies) then
!     allocate(orphans_local_pid(1:npastprogs+nprogs))
!     orphans_local_pid = 0
!     allocate(prog_galaxy_local_id(1:nprogs))
!     prog_galaxy_local_id = 0
!   endif
!
!
!
!   !---------------------------------
!   ! Sort arrays for quick matching
!   !---------------------------------
!
!   ! Sort arrays for matching
!   ! call quick_sort_int_int(idp_copy, part_local_ind, npartmax)
!   call quick_sort_int_int(tracers_all, tracer_loc_progids_all, nprogs*nmost_bound)
!   call quick_sort_int_int(galaxy_tracers_copy, dummy, nprogs)
!
!   ! sort past progenitors in-place by galaxy particle ID;
!   ! will be needed later for multi-snapshot progenitor search
!   call quick_sort_int_int(pmprogs_galaxy, sort_ind_past, npastprogs)
!
!   dummy(1:npastprogs) = pmprogs(1:npastprogs)
!   do i = 1, npastprogs
!     pmprogs(i) = dummy(sort_ind_past(i))
!   enddo
!
!   dummy(1:npastprogs) = pmprogs_t(1:npastprogs)
!   do i = 1, npastprogs
!     pmprogs_t(i) = dummy(sort_ind_past(i))
!   enddo
!
!   dummy_real(1:npastprogs) = pmprogs_mass(1:npastprogs)
!   do i = 1, npastprogs
!     pmprogs_mass(i) = dummy_real(sort_ind_past(i))
!   enddo
!
!   dummy_real(1:npastprogs) = pmprogs_mpeak(1:npastprogs)
!   do i = 1, npastprogs
!     pmprogs_mass(i) = dummy_real(sort_ind_past(i))
!   enddo
!
!   deallocate(dummy, dummy_real, sort_ind_past)
!
!
!
!
!   !------------------------
!   ! find starting indices
!   !------------------------
!
!   itrace = 0; ipart = 0; igalaxy = 0; iprog = 0;
!
!   ! do i = 1, npartmax
!   !   if (idp_copy(i) > 0) then
!   !     ipart = i
!   !     exit
!   !   endif
!   ! enddo
!
!   do i = 1, nprogs*nmost_bound
!     if (tracers_all(i) > 0) then
!       itrace = i
!       exit
!     endif
!   enddo
!
!   igalaxy = 1
!   ipastprog = 1
!
!   ntracers = 0
!
!
!
!
!   !--------------------------------------------------------------
!   ! Identify which progenitor tracer particles are on this CPU.
!   ! loop over all local particles and all progenitor particles.
!   ! at this point, both arrays are sorted.
!   ! Raise the array index of array which has lower particle ID
!   ! to find matches.
!   !--------------------------------------------------------------
!
!   do while ( ipart <= npartmax)
!
!     !-------------------------------------------------------------------------------------
!     ! Check for tracers and past progenitor galaxies while itrace <= nprogs*nmost_bound
!     !-------------------------------------------------------------------------------------
!
!     do while (itrace <= nprogs * nmost_bound .and. ipart <= npartmax)
!       ! if particles aren't a match, raise the index where lower ID is
!       if (tracers_all(itrace) < idp_copy(ipart)) then
!         itrace = itrace + 1
!
!       else if (tracers_all(itrace) > idp_copy(ipart)) then
!         ! before raising local particle index, check whether you own a past progenitor
!         ! no need to check whether npastprogs > 0: arrays are allocated
!         ! (1:npastprogs+nprogs) to have extra space in case new progs
!         ! need to be added to the list
!         do while (ipastprog <= npastprogs)
!           if (pmprogs_galaxy(ipastprog) < idp_copy(ipart)) then
!               ipastprog = ipastprog + 1
!           else if (pmprogs_galaxy(ipastprog) == idp_copy(ipart)) then
!             ! you found a match!
!             pmprogs_owner(ipastprog) = myid
!             if (make_mock_galaxies) orphans_local_pid(ipastprog) = part_local_ind(ipart)
!             ipastprog = ipastprog + 1
!           else
!             exit
!           endif
!         enddo
!         ipart = ipart + 1
!
!
!       else
!         !----------------
!         ! found a match!
!         !----------------
!         iprog = tracer_loc_progids_all(itrace)
!         ! count the tracer for corresponding prog
!         ! add to list of this CPU
!         ntracers = ntracers + 1                                   ! count one more tracer on this cpu
!         tracers_local_pid_long(ntracers) = part_local_ind(ipart)  ! save local particle ID
!         tracer_local_ids_long(ntracers) = iprog                   ! save local prog ID it belongs to
!
!         ! check if found tracer is also a galaxy:
!         do while (igalaxy <= nprogs)
!           if (galaxy_tracers_copy(igalaxy) < tracers_all(itrace)) then
!             igalaxy = igalaxy + 1
!           else if (galaxy_tracers_copy(igalaxy) == tracers_all(itrace)) then
!             prog_owner(iprog) = myid
!             if (make_mock_galaxies) prog_galaxy_local_id(iprog) = part_local_ind(ipart)   ! store local ID of galaxy particle
!             igalaxy = igalaxy + 1
!             exit
!           else
!             exit
!           endif
!         enddo
!
!         ipart = ipart + 1
!         itrace = itrace + 1
!       endif
!
!     enddo
!
!
!     ! reset ipart in case it reached the limit in the previous loop
!     if (ipart>npartmax) ipart = npartmax
!
!     !---------------------------------------------------------------------------
!     ! If there are no more tracers to check for, check only for past galaxies
!     !---------------------------------------------------------------------------
!     do while (ipastprog <= npastprogs)
!       if (pmprogs_galaxy(ipastprog) < idp_copy(ipart)) then
!           ipastprog = ipastprog + 1
!       else if (pmprogs_galaxy(ipastprog) == idp_copy(ipart)) then
!         ! you found a match!
!         pmprogs_owner(ipastprog) = myid
!         if (make_mock_galaxies) orphans_local_pid(ipastprog) = part_local_ind(ipart)
!         ipastprog = ipastprog + 1
!       else
!         exit
!       endif
!     enddo
!
!     if (ipastprog > npastprogs) exit
!     ipart = ipart + 1
!
!   enddo
!
!
!   deallocate(idp_copy, galaxy_tracers_copy, part_local_ind)
!
!
!
!
!   !---------------------------------------------------
!   ! Clean up local tracers information for this CPU
!   !---------------------------------------------------
!
!   allocate(tracers_loc_pid(1:ntracers))
!   allocate(tracer_loc_progids(1:ntracers))
!
!   tracers_loc_pid(1:ntracers) = tracers_local_pid_long(1:ntracers)
!   tracer_loc_progids(1:ntracers) = tracer_local_ids_long(1:ntracers)
!
!   ! Deallocate auxilliary arrays
!   deallocate(tracer_local_ids_long, tracers_local_pid_long)
!
!   ! Deallocate global arrays that aren't used anymore
!   deallocate(tracers_all)
!   deallocate(tracer_loc_progids_all)
!
!
!
!
! #ifndef WITHOUTMPI
!   !--------------------------------
!   ! communicate progenitor owners
!   !--------------------------------
!
!   ! first the actual progenitors
!   allocate(local_owners_info(1:nprogs))     ! progenitor owners array, local to each cpu
!   local_owners_info = prog_owner
!   call MPI_ALLREDUCE(local_owners_info, prog_owner, nprogs, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, i)
!   deallocate(local_owners_info)
!
!   ! then the past progenitors
!   allocate(local_owners_info(1:npastprogs)) ! progenitor owners array, local to each cpu
!   local_owners_info = pmprogs_owner(1:npastprogs)
!   call MPI_ALLREDUCE(local_owners_info, pmprogs_owner, npastprogs, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, i)
!   deallocate(local_owners_info)
! #endif
!
!
!   return
!
! end subroutine process_progenitor_data
!
!
!
!
!




!====================================
subroutine read_progenitor_data()
!====================================
  
  !---------------------------------
  ! Reads in all progenitor data
  !---------------------------------

  implicit none

  integer           :: prog_read, prog_read_local, startind, tracer_free
  integer           :: nprogs_to_read, progcount_to_read, np
  integer           :: iprog, i
  character(LEN=80) :: fileloc
  character(LEN=5)  :: output_to_string
  logical           :: exists

  integer, allocatable, dimension(:) :: read_buffer       ! temporary array for reading in data
  real(dp),allocatable, dimension(:) :: read_buffer_2     ! temporary array for reading in data
  real(dp),allocatable, dimension(:) :: read_buffer_mpeak ! temporary array for reading in mock galaxy data

#ifndef WITHOUTMPI
  integer, dimension (1:MPI_STATUS_SIZE) :: state
  integer, dimension(1:4)                :: buf
  integer                                :: mpi_err, filehandle
#endif

  if (verbose) write(*,*) " Calling read progenitor data."

  call title(ifout-1, output_to_string)
  ! ifout -1: read from previous output!
  nprogs = 0
  nprogs_to_read = 0
  progcount_to_read = 0



  !========================
  ! Read progenitor counts
  !========================

  if (myid == 1) then ! read in stuff
    
    !-------------------------------------------------------------
    ! Current progenitor counts
    ! Both of these count files need to be present in any case.
    !-------------------------------------------------------------

    fileloc=TRIM('output_'//TRIM(output_to_string)//'/progenitorcount.dat')

    ! check that file exists
    inquire(file=fileloc, exist=exists)
    if (.not.exists) then
      write(*, *) "ID", myid, "didn't find file ", fileloc
      stop
    endif


    open(unit=666,file=fileloc,form='unformatted')
    read(666) nprogs, nprogs_to_read, progcount_to_read, npastprogs
    ! open(unit=666,file=fileloc,form='formatted')
    ! read(666, '(2(I7,x))') nprogs, nprogs_to_read
    close(666)


  endif


#ifndef WITHOUTMPI
  buf = (/nprogs, progcount_to_read, nprogs_to_read, npastprogs/)
  call MPI_BCAST(buf, 4, MPI_INTEGER, 0, MPI_COMM_WORLD, mpi_err)
  nprogs = buf(1)
  progcount_to_read = buf(2)
  nprogs_to_read = buf(3)
  npastprogs = buf(4)
#endif
  

  write(*,*) "READ IN:"
  write(*,*) "nprogs", nprogs
  write(*,*) "progcount_to_read", progcount_to_read
  write(*,*) "nprogs_to_read", nprogs_to_read
  write(*,*) "npastprogs", npastprogs



  !==================================
  ! READ CURRENT PROGENITOR DATA
  !==================================

  !---------------------------
  ! Allocate arrays
  !---------------------------

  allocate(prog_id(1:nprogs))
  prog_id = 0       ! list of progenitor global IDs
  prog_free = 1     ! first free local progenitor ID

  allocate(prog_owner(1:nprogs))
  prog_owner = -1

  allocate(prog_mass(1:nprogs))
  prog_mass = 0

  allocate(tracers_all(1:nprogs*nmost_bound))
  tracers_all = 0

  allocate(galaxy_tracers(1:nprogs))

  allocate(tracer_loc_progids_all(1:nprogs*nmost_bound))
  tracer_loc_progids_all = 0
  tracer_free = 1   ! first free local tracer index

  if (make_mock_galaxies) then
    i = nprogs
  else
    ! just to prevent "may be uninitialized" warnings
    i = 1
  endif
  allocate(prog_mpeak(1:i))
  prog_mpeak = 0


  if (nprogs > 0) then

    !--------------------------------
    ! Read progenitor particles
    !--------------------------------

    fileloc=TRIM('output_'//TRIM(output_to_string)//'/progenitor_data.dat')

    inquire(file=fileloc, exist=exists)
    if (.not.exists) then
      write(*,*) "ID", myid, "didn't find file ", fileloc
      stop
    endif


    allocate(read_buffer(1:progcount_to_read))

#ifndef WITHOUTMPI
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fileloc, &
      MPI_MODE_RDONLY, MPI_INFO_NULL,filehandle, mpi_err)
    call MPI_FILE_READ(filehandle, read_buffer, &
      progcount_to_read, MPI_INTEGER, state, mpi_err)
    call MPI_FILE_CLOSE(filehandle, mpi_err)
#else
    open(unit=666,file=fileloc,form='unformatted')
    read(666) read_buffer
    close(666)
#endif





    !--------------------------------
    ! Read progenitor masses
    !--------------------------------

    fileloc=TRIM('output_'//TRIM(output_to_string)//'/progenitor_mass.dat')

    inquire(file=fileloc, exist=exists)
    if (.not.exists) then
      write(*, *) "ID", myid, "didn't find file ", fileloc
      stop
    endif


    allocate(read_buffer_2(1:nprogs_to_read))

    ! just to prevent "may be uninitialized" warnings
    if (make_mock_galaxies) then
      i = nprogs_to_read
    else
      i = 1
    endif
    allocate(read_buffer_mpeak(1:i))
    read_buffer_mpeak = 0

#ifndef WITHOUTMPI
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fileloc, MPI_MODE_RDONLY, MPI_INFO_NULL, filehandle, mpi_err)
    call MPI_FILE_READ(filehandle, read_buffer_2, nprogs_to_read, MPI_DOUBLE_PRECISION, state, mpi_err)
    if (make_mock_galaxies) then
      call MPI_FILE_READ(filehandle, read_buffer_mpeak, nprogs_to_read, MPI_DOUBLE_PRECISION, state, mpi_err)
    endif
    call MPI_FILE_CLOSE(filehandle, mpi_err)
#else
    open(unit=666,file=fileloc,form='unformatted')
    read(666) read_buffer_2
    if (make_mock_galaxies) then
      read(666) read_buffer_mpeak
    endif
    close(666)
#endif




    !----------------------------------
    ! Sort out the data you just read
    !----------------------------------

    tracer_free = 1
    iprog = 1
    startind = 1

    write(*,*) "sorting out data"
    write(*,*) "==============================="

    do while (startind <= progcount_to_read)

      ! write(*,*) "SI", startind, "iprog", iprog

      prog_read = read_buffer(startind)
      np = read_buffer(startind + 1)
      write(*,*) "Prog read", prog_read, "np", np, "iprog", iprog

      ! get local instead global ID in prog_read (past tense "read")
      call get_local_prog_id(prog_read, prog_read_local)

      ! prog_mass(prog_read_local) = read_buffer_2(iprog)
      ! if (make_mock_galaxies) then
      !   prog_mpeak(prog_read_local) = read_buffer_mpeak(iprog)
      ! endif

      write(*, '(A10)', advance='no') "particles:"
      do i = startind+2, startind+1+np
        if (read_buffer(i) > 0) then
          tracers_all(tracer_free) = read_buffer(i)               ! add new tracer particle
          tracer_loc_progids_all(tracer_free) = prog_read_local   ! write which progenitor tracer belongs to
          tracer_free = tracer_free + 1                           ! raise index for next tracer
        else 
          ! found a galaxy particle
          tracers_all(tracer_free) = -read_buffer(i)              ! add new tracer particle
          galaxy_tracers(prog_read_local) = -read_buffer(i)       ! add new galaxy tracer
          tracer_loc_progids_all(tracer_free) = prog_read_local   ! write which progenitor tracer belongs to
          tracer_free = tracer_free + 1                           ! raise index for next tracer
        endif
        write(*,'(I10)', advance='no') read_buffer(i)
      enddo
      write(*,*)

      iprog = iprog + 1
      startind = startind + 2 + np

    enddo

    deallocate(read_buffer, read_buffer_2)
    if (make_mock_galaxies) deallocate(read_buffer_mpeak)

  endif ! nprogs > 0




  !========================================
  ! READ PAST PROGENITOR DATA
  !========================================
  
  !-------------------------
  ! Allocate arrays
  !-------------------------

  ! overestimate size to fit new ones if necessary
  npastprogs_max = npastprogs + nprogs
  pmprog_free = npastprogs + 1

  ! Past Merged Progenitors for multi-snapshot matching
  allocate(pmprogs(1:npastprogs_max))
  pmprogs = 0 

  ! Current owner Merged Progenitors' galaxy particles
  allocate(pmprogs_owner(1:npastprogs_max))
  pmprogs_owner = 0

  ! Past Merged Progenitors' galaxy particles
  allocate(pmprogs_galaxy(1:npastprogs_max))
  pmprogs_galaxy = 0

  ! Time at which past progenitors have been merged (= ifout at merging time)
  allocate(pmprogs_t(1:npastprogs_max))
  pmprogs_t = 0

  ! past merged progenitor mass
  allocate(pmprogs_mass(1:npastprogs_max))
  pmprogs_mass = 0

  ! mock galaxy stuff
  if (make_mock_galaxies) then
    allocate(pmprogs_mpeak(1:npastprogs_max))
    pmprogs_mpeak = 0
  endif
  

  if (npastprogs > 0) then

    !-----------------------------
    ! Read in data
    !-----------------------------

    allocate(read_buffer(1:3*npastprogs))

      fileloc=TRIM('output_'//TRIM(output_to_string)//'/past_merged_progenitors.dat')

#ifndef WITHOUTMPI
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fileloc, &
      MPI_MODE_RDONLY, MPI_INFO_NULL,filehandle, mpi_err)
    call MPI_FILE_READ(filehandle, read_buffer, &
      npastprogs*3, MPI_INTEGER, state, mpi_err)
    call MPI_FILE_CLOSE(filehandle, mpi_err)
#else
    open(unit=666,file=fileloc,form='unformatted')
    read(666) read_buffer
    close(666)
#endif


    !---------------------------------
    ! Read past progenitor's masses
    !---------------------------------

    fileloc=TRIM('output_'//TRIM(output_to_string)//'/past_merged_progenitor_mass.dat')

    inquire(file=fileloc, exist=exists)
    if (.not.exists) then
      write(*, *) "ID", myid, "didn't find file ", fileloc
      stop
    endif

#ifndef WITHOUTMPI
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fileloc, MPI_MODE_RDONLY, MPI_INFO_NULL, filehandle, mpi_err)
    call MPI_FILE_READ(filehandle, pmprogs_mass(1:npastprogs), npastprogs, MPI_DOUBLE_PRECISION, state, mpi_err)
    if (make_mock_galaxies) then
      call MPI_FILE_READ(filehandle, pmprogs_mpeak(1:npastprogs), npastprogs, MPI_DOUBLE_PRECISION, state, mpi_err)
    endif
    call MPI_FILE_CLOSE(filehandle, mpi_err)
#else
    open(unit=666,file=fileloc,form='unformatted')
    read(666) pmprogs_mass(1:npastprogs)
    if (make_mock_galaxies) then
      read(666) pmprogs_mpeak(1:npastprogs)
    endif
    close(666)
#endif



    !----------------------------------
    ! Sort out the data you just read
    !----------------------------------

    pmprog_free = 1
    iprog = 1
    do while (iprog <= 3*npastprogs)
      pmprogs(pmprog_free) = read_buffer(iprog)
      pmprogs_galaxy(pmprog_free) = read_buffer(iprog + 1)
      pmprogs_t(pmprog_free) = read_buffer(iprog + 2)
      iprog = iprog + 3
      pmprog_free = pmprog_free + 1
    enddo

    deallocate(read_buffer)

  endif ! npastprogs > 0

end subroutine read_progenitor_data










!=================================================
subroutine get_local_prog_id(global_id, local_id)
!=================================================
  
  implicit none
  integer, intent(in)  :: global_id 
  integer, intent(out) :: local_id 
  integer              :: i

  !---------------------------------
  ! Gets the local progenitor ID
  ! given the global progenitor ID.
  !---------------------------------

  ! Check if prog is already included
  do i = 1, prog_free-1
    if (prog_id(i) == global_id) then
      local_id = i
      return
    endif
  enddo

  ! progenitor is not in there; Give him a new index
  i = prog_free ! save for later
  prog_id(i) = global_id
  prog_free = prog_free + 1
  local_id = i
  return

end subroutine get_local_prog_id 





!==================================================
subroutine fill_matrix(mat, key, tar, np, act)
!==================================================

  !----------------------------------------------------
  ! Performs given operation on matrix elements of
  ! prog-desc matrices.
  !
  ! clump key:    local id for progenitor (p2d_links),
  !               local id for descendant (d2p_links)
  ! clump target: local id for progenitor (p2d_links),
  !               global id for descendant (d2p_links)
  ! act:          what to do with matrix elements 
  !   act = 'add': add up
  !   act = 'set': don't add, just set
  !   act = 'inv': make value negative
  !   'add' and 'set' will introduce new entries if
  !   there wasn't one before, 'inv' won't.
  !----------------------------------------------------

  implicit none

  integer,             intent(in)    :: key, tar, np
  character(len=3),    intent(in)    :: act
  type(prog_desc_mat), intent(inout) :: mat

  integer :: i, itar

  ! if there is no first value, 
  ! add new first value
  if (mat%cnt(key) == 0) then
    if (act/='inv') then
      mat%first(key) = mat%mat_free_ind
      mat%clmp_id(mat%mat_free_ind) = tar
      mat%cnt(key) = mat%cnt(key) + 1

      if (act=='add') then
        mat%ntrace(mat%mat_free_ind) = mat%ntrace(mat%mat_free_ind) + np
      else if (act=='set') then
        mat%ntrace(mat%mat_free_ind) = np
      endif

      mat%mat_free_ind = mat%mat_free_ind + 1
    endif
    return

  else
    ! Try to find a match first
    itar = mat%first(key)
    do i = 1, mat%cnt(key)
      ! If you found a match:
      if (mat%clmp_id(itar) == tar) then
        if (act=='add') then
          mat%ntrace(itar) = mat%ntrace(itar) + np
        else if (act=='set') then
          mat%ntrace(itar) = np
        else if (act=='inv') then
          mat%ntrace(itar) = -mat%ntrace(itar)
        endif
        return
      endif

      if (mat%next(itar) > 0) then
        ! cycle
        itar = mat%next(itar)
      else
        exit
      endif
    enddo
  endif



  ! if you didn't find anything, add new value
  if (act /= 'inv') then  
    mat%next(itar) = mat%mat_free_ind
    mat%clmp_id(mat%mat_free_ind) = tar
    mat%cnt(key) = mat%cnt(key) + 1
    if (act=='add') then
      mat%ntrace(mat%mat_free_ind) = mat%ntrace(mat%mat_free_ind) + np
    else if (act=='set') then
      mat%ntrace(mat%mat_free_ind) = np
    endif

    mat%mat_free_ind = mat%mat_free_ind + 1
  endif
  return

end subroutine fill_matrix

SUBROUTINE quick_sort_int_int(list, order, n)

  !------------------------------------------------------------
  ! Sort array of integers (list), rearrange array of integers
  ! (order) in the same way.
  !------------------------------------------------------------
  
  IMPLICIT NONE
  ! Quick sort routine from:
  ! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
  ! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
  ! Modified to sort the second given array by the same rules.


  INTEGER :: n
  INTEGER, DIMENSION (1:n), INTENT(INOUT)  :: list
  INTEGER, DIMENSION (1:n), INTENT(INOUT)  :: order


  CALL quick_sort_1_int_int(list, order, 1, n, n)

END SUBROUTINE quick_sort_int_int

RECURSIVE SUBROUTINE quick_sort_1_int_int(list, order, left_end, right_end, n)

  INTEGER, INTENT(IN) :: left_end, right_end, n
  INTEGER, DIMENSION (1:n), INTENT(INOUT)  :: list
  INTEGER, DIMENSION (1:n), INTENT(INOUT)  :: order

  !     Local variables
  INTEGER             :: i, j, itemp
  INTEGER             :: reference, temp
  INTEGER, PARAMETER  :: max_simple_sort_size = 6

  IF (right_end < left_end + max_simple_sort_size) THEN
     ! Use interchange sort for small lists
     CALL interchange_sort_int_int(list, order, left_end, right_end, n)

  ELSE
     ! Use partition ("quick") sort
     reference = list((left_end + right_end)/2)
     i = left_end - 1; j = right_end + 1

    DO
        ! Scan list from left end until element >= reference is found
        DO
           i = i + 1
           IF (list(i) >= reference) EXIT
        END DO
        ! Scan list from right end until element <= reference is found
        DO
           j = j - 1
           IF (list(j) <= reference) EXIT
        END DO


        IF (i < j) THEN
           ! Swap two out-of-order elements
           temp = list(i); list(i) = list(j); list(j) = temp
           itemp = order(i); order(i) = order(j); order(j) = itemp
        ELSE IF (i == j) THEN
           i = i + 1
           EXIT
        ELSE
           EXIT
        END IF
     END DO
     IF (left_end < j) CALL quick_sort_1_int_int(list, order, left_end, j, n)
     IF (i < right_end) CALL quick_sort_1_int_int(list, order, i, right_end, n)
  END IF

END SUBROUTINE quick_sort_1_int_int


SUBROUTINE interchange_sort_int_int(list, order, left_end, right_end, n)

  INTEGER, INTENT(IN) :: left_end, right_end, n
  INTEGER, DIMENSION (1:n), INTENT(INOUT)  :: list
  INTEGER, DIMENSION (1:n), INTENT(INOUT)  :: order

  !     Local variables
  INTEGER             :: i, j, itemp, temp

  DO i = left_end, right_end - 1
     DO j = i+1, right_end
        IF (list(i) > list(j)) THEN
           temp = list(i); list(i) = list(j); list(j) = temp
           itemp = order(i); order(i) = order(j); order(j) = itemp
        END IF
     END DO
  END DO

END SUBROUTINE interchange_sort_int_int

!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character(LEN=5)::nchar

  character(LEN=1)::nchar1
  character(LEN=2)::nchar2
  character(LEN=3)::nchar3
  character(LEN=4)::nchar4
  character(LEN=5)::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif


end subroutine title






end program read_progenitors
