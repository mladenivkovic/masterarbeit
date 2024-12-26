!========================================================
! Find associated clump ID for all orphan galaxies.
!========================================================


!===================================
program get_clump_ID_of_orphans
!===================================

  !==============================================
  ! cmd line args:
  ! get_clump_ID_of_orphans outputnr
  ! outputnr: 67 for output_00067
  !==============================================

  use omp_lib
  use constants_and_parameters
  use io_module
  use quick_sort

  implicit none

  real(dp) :: mthresh = 1d10 ! M_Sol

  !-------------------------
  ! Arrays
  !-------------------------
  ! read-in data arrays
  integer                               :: pguess = 5*2000000 ! initial guess for number of particles. Take npartmax!

  integer , dimension(:), allocatable   :: order    ! sorted orphan order by ID
  integer , dimension(:), allocatable   :: oclumpid ! orphan clump ID
  real(dp),dimension(:, :), allocatable :: xo_sorted
  real(dp),dimension(:), allocatable    :: mo_sorted

  ! character(len=80) :: fname    ! output filename
  integer :: i, nohalos

  !---------------------------
  ! Set-up
  !---------------------------

  call get_srcdir_from_cmdlineargs() ! from io_module
  call read_info() ! from io_module
  call read_galaxy_data(ignore_centrals_satellite_distinction=.true.)
  call filter_orphans()

  ! sort orphans by their ID
  allocate(order(1:norphans))
  do i = 1, norphans
    order(i) = i
  enddo
  call quick_sort_int(ido, order, norphans)
  allocate(xo_sorted(1:norphans, 1:3))
  allocate(mo_sorted(1:norphans))
  do i = 1, norphans
    xo_sorted(i, 1) = xo(order(i), 1)
    xo_sorted(i, 2) = xo(order(i), 2)
    xo_sorted(i, 3) = xo(order(i), 3)
    mo_sorted(i) = mo(order(i))
  enddo

  call cleanup_duplicate_orphans()


  allocate(oclumpid(1:norphans))
  oclumpid = -1 

  call print_infos()
  call read_particles() ! the important part here!
  call write_clumpids() ! write results

  ! safety checks
  nohalos = 0
  do i = 1, norphans
    if (oclumpid(i) == 0) nohalos = nohalos + 1
    if (oclumpid(i) < 0) write(*, *) "Found orphan clump id < 0?", i, oclumpid(i), xo_sorted(i, 1:3)
  enddo

  write(*, *) "Found", nohalos, "/", norphans, "orphans outside of clumps"
  write(*,*) "get_clump_ID_of_orphans.f03 finished."

  
  call deallocate_io_module()
  deallocate(order, oclumpid)
  deallocate(xo_sorted, mo_sorted)



contains


  subroutine filter_orphans()
    
    !------------------------------------
    ! Apply mass threshold to orphans
    ! in order to save some time
    !------------------------------------

    implicit none
    integer :: i,ikeep, counts
    real(dp), allocatable, dimension(:, :) :: xokeep
    real(dp), allocatable, dimension(:)    :: mokeep
    integer(i8),allocatable,dimension(:)   :: idokeep

    counts = 0
    do i = 1, norphans
      if (mo(i) >= mthresh) then
        counts = counts + 1
      endif
    enddo

    allocate(xokeep(1:counts, 1:3))
    allocate(mokeep(1:counts))
    allocate(idokeep(1:counts))

    ikeep = 0
    do i = 1, norphans
      if (mo(i) >= mthresh) then
        ikeep = ikeep + 1
        xokeep(ikeep, 1:3) = xo(i, 1:3)
        mokeep(ikeep) = mo(i)
        idokeep(ikeep) = ido(i)
      endif
    enddo

    deallocate(xo)
    allocate(xo(1:counts, 1:3))
    xo(:, :) = xokeep(:, :)
    deallocate(mo)
    allocate(mo(1:counts))
    mo(:) = mokeep(:)
    deallocate(ido)
    allocate(ido(1:counts))
    ido(:) = idokeep(:)

    norphans = counts
    write(*, *) "Keeping", norphans, "after applying mass threshold"

  end subroutine filter_orphans




  subroutine cleanup_duplicate_orphans()

    !---------------------------------------------
    ! In versions of ramses before 11.05.2021, 
    ! a bug was present where orphan galaxies
    ! could be added twice to the list of orphans.
    ! Remove the duplicate with the lower mass.
    !---------------------------------------------

    implicit none
    integer :: i, j, duplicates

    i = 1
    duplicates = 0
    do while (i < norphans - 1) 
      if (ido(i) == ido(i+1)) then
        ! found duplicate
        duplicates = duplicates + 1
        if (mo_sorted(i+1) > mo_sorted(i)) then
          ! keep the more massive one
          mo_sorted(i) = mo_sorted(i+1)
        endif
        ! move everything that comes after in the array
        ! one position back
        do j = i+1, norphans - 1
          ido(j) = ido(j+1)
          xo_sorted(j, 1) = xo_sorted(j+1, 1)
          xo_sorted(j, 2) = xo_sorted(j+1, 2)
          xo_sorted(j, 3) = xo_sorted(j+1, 3)
          mo_sorted(j) = mo_sorted(j+1)
        enddo
        norphans = norphans - 1
      endif
      i = i + 1
    enddo

    write(*, *) "removed ", duplicates, "duplicates" 


  end subroutine cleanup_duplicate_orphans



  subroutine print_infos()

    implicit none
    write(*,*) "Working with parameters:"
    write(*,'(A20,x,A25)') "  srcdir:        ", srcdir
    write(*,'(A20,x,I25)') "  ncpu:          ", ncpu
    write(*,'(A20,x,I25)') "  ngalaxies:     ", ngalaxies
    write(*,'(A20,x,I25)') "  norphans:      ", norphans
    ! write(*,'(A20,x,I25)')      "  levelmin:      ", levelmin
    ! write(*,'(A20,x,I25)')      "  levelmax:      ", levelmax
    ! write(*,'(A20,x,E25.3,x,F25.3,x,A9)') "  unit_l:        ", unit_l, (unit_l/Mpc), "[Mpc]"
    ! write(*,'(A20,x,F25.3)')    "  a_exp:         ", aexp
    ! write(*,'(A20,x,F25.3)')    "  z:             ", 1.d0/aexp - 1.d0
    ! write(*,'(A20,x,F25.3)')    "  H0:            ", H0

  end subroutine print_infos






  subroutine read_particles()
    !---------------------------------------------
    ! Read in all particle data, file by file,
    ! and get clump ID of orphan particles
    !---------------------------------------------

    implicit none
    character(len=80)                     :: fname
    character(len=5)                      :: outputnr, cpunr
    integer                               :: i, j, id, ifile, nparts_loc, finished
    integer                               :: nparts_tot
    integer,  allocatable, dimension(:)   :: pclmpid_loc
    integer(i8), allocatable, dimension(:):: id_loc
    real :: fract
    
    outputnr = srcdir(8:12)
    finished = 0
    nparts_tot = (2**levelmin)**3

    write(*,*) "Reading in particle data. This may take a while."



    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(id,i,j,ifile,fname,cpunr,fract,&
    !$OMP         nparts_loc,id_loc,pclmpid_loc) 
      id = omp_get_thread_num()

      allocate(id_loc(1:pguess))
      id_loc = 0
      allocate(pclmpid_loc(1:pguess))
      pclmpid_loc = 0

      !$OMP DO
        do ifile=1, ncpu
          call title(ifile, cpunr)
          fname=TRIM(srcdir//'/part_'//outputnr//'.out'//cpunr)
          open(unit=600+id, file=fname, form='unformatted')
          read(600+id) !ncpu
          read(600+id) !ndim
          read(600+id) nparts_loc
          read(600+id) !localseed
          read(600+id) !nstar_tot
          read(600+id) !mstar_tot
          read(600+id) !mstar_lost
          read(600+id) !nsink

          read(600+id) !xp_loc(1:nparts_loc,1)
          read(600+id) !xp_loc(1:nparts_loc,2)
          read(600+id) !xp_loc(1:nparts_loc,3)
          read(600+id) !vx
          read(600+id) !vy
          read(600+id) !vz
          read(600+id) !mass
          read(600+id) id_loc(1:nparts_loc)
          close(600+id)

          fname=TRIM(srcdir//'/unbinding_'//outputnr//'.out'//cpunr)
          open(unit=600+id, file=fname, form='unformatted')
          read(600+id) pclmpid_loc(1:nparts_loc)
          close(600+id)

          do i=1, nparts_loc
            do j=1, norphans
              if (ido(j) > id_loc(i)) then
                exit 
              else if (ido(j) == id_loc(i)) then
                oclumpid(j) = abs(pclmpid_loc(i))
                exit
              endif
            enddo
          enddo

          !$OMP CRITICAL
            finished = finished + nparts_loc
            fract = real(finished) / real(nparts_tot) * 100
            write(*, '(F6.2,x,A)') fract, "%"
          !$OMP END CRITICAL

        enddo !loop over files
      !$OMP END DO

      deallocate(id_loc, pclmpid_loc)
    !$OMP END PARALLEL

  end subroutine read_particles



  subroutine write_clumpids()

    !---------------------------------------------------
    ! write the results (clump IDs of orphans) to file
    !---------------------------------------------------

    implicit none

    character(len=80) :: fname
    integer :: i

    fname = 'orphan_clump_IDs-'//srcdir(8:12)//'.txt'
    open(unit=666, file=fname, form='formatted')
    write(666, '(6(A15,x))') "# orphan_ID", "clump_ID", "x", "y", "z", "mass"
    do i = 1, norphans
      write(666, '(2(I15,x),4(E15.6E3,x))') ido(i), oclumpid(i), &
          xo_sorted(i, 1), xo_sorted(i, 2), xo_sorted(i, 3), mo_sorted(i)
    enddo
    close(666)

  end subroutine write_clumpids


end program
