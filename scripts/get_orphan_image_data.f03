!===========================================
! Get an image of a halo of your choosing
! and all the orphans and galaxies within it
!===========================================


!===================================
program get_orphan_image_data
!===================================

  !========================================================
  ! cmd line args:
  ! get_halostats.o outputnr halo halo_x halo_y halo_z
  ! outputnr: 67 for output_00067
  ! halo:     halo ID
  ! halo_x, halo_y, halo_z : peak positions of halo, used
  !                          to center image
  !========================================================

  use omp_lib
  use constants_and_parameters
  use io_module
  use density_projection
  implicit none

  !-------------------------------------------
  ! Global variables and parameters
  !-------------------------------------------

  real(dp) :: image_width_user = 5. ! Mpc
  real(dp) :: mthresh = 1.d10 ! M_Sol mass threshold


  !------------------
  ! Read-in data
  !------------------
  integer             :: nparts, nparts_tot
  integer             :: ng_write, no_write ! number of galaxies/orphans in halo
  integer             :: halo               ! the halo to work with


  !-------------------------
  ! Arrays
  !-------------------------
  ! read-in data arrays
  integer                               :: pguess = 5*2000000 ! initial guess for number of particles. Take npartmax!
  real(dp), dimension(:,:), allocatable :: xow, xgw           ! positions of orphans/galaxies in halo  (w: write)
  real(dp), dimension(:), allocatable   :: mow, mgw           ! masses of orphans/galaxies in halo  (w: write)

  real(dp), dimension(:,:), allocatable :: xp     ! positions of particles [code units]
  real(dp), dimension(:), allocatable   :: mp     ! masses of particles [code units]


  ! real(dp) :: xmin, xmax, ymin, ymax, zmin, zmax ! image borders; defined in density_projection
  ! whether x,y, or z have been moved for periodic boundary corrections
  logical   :: moved_x = .false.
  logical   :: moved_y = .false.
  logical   :: moved_z = .false.
  

  !---------------------------
  ! Set-up
  !---------------------------

  ! imported from density_projection module
  ! interpolation = 'ngp'    ! which interpolation to use: cic for cloud in cell, ngp for nearest grid point
  interpolation = 'cic'   ! which interpolation to use: cic for cloud in cell, ngp for nearest grid point
  nc = 1000                ! number of cells per dimension for image


  call get_cmdlineargs(halo) ! from io_module
  call read_info() ! from io_module
  write(*, *) "Working for halo", halo
  call find_image_borders()

  ! read galaxy data first in case we can exit early
  call read_galaxy_data(.true.)
  call filter_galaxy_data()

  call read_particles()

  call allocate_density_fields()
  allocate(mp(1:nparts))
  mp = 1.d0 / (2**(levelmin)) ** 3
  call get_density_projection(xp, mp, nparts)


  call write_image()

  write(*,*) "get_orphan_image_data.f03 finished."

  call deallocate_density_fields()
  call deallocate_io_module()
  deallocate(xp, mp)
  deallocate(xow, xgw, mow, mgw)



contains


  !========================================
  subroutine find_image_borders()
  !========================================
    
    !---------------------------------
    ! Find xmin, xmax, ... etc 
    !---------------------------------

    implicit none
    real(dp) :: dx_code_units

    dx_code_units = 0.5 * image_width_user / unit_l_Mpc

    xmin = CoM(1) - dx_code_units
    xmax = CoM(1) + dx_code_units
    ymin = CoM(2) - dx_code_units
    ymax = CoM(2) + dx_code_units
    zmin = CoM(3) - dx_code_units
    zmax = CoM(3) + dx_code_units

    ! if periodicity correction necessary, 
    ! always move from negative to positive side
    if (xmin < 0.d0) then
      moved_x = .true.
      xmin = 1.d0 + xmin
      xmax = 1.d0 + xmax
    endif
    if (xmax > 1.d0) then
      moved_x = .true.
    endif
    if (ymin < 0.d0) then
      moved_y = .true.
      ymin = 1.d0 + ymin
      ymax = 1.d0 + ymax
    endif
    if (ymax > 1.d0) then
      moved_y = .true.
    endif
    if (zmin < 0.d0) then
      moved_z = .true.
      zmin = 1.d0 + zmin
      zmax = 1.d0 + zmax
    endif
    if (zmax > 1.d0) then
      moved_z = .true.
    endif

    ! write(*, '(A,x,6F12.8)') "limits:", &
    !     xmin, xmax, ymin, ymax, zmin, zmax


  end subroutine find_image_borders



  
  !========================================
  subroutine filter_galaxy_data()
  !========================================

    !---------------------------------------------
    ! Find which galaxies and orphans to keep
    !---------------------------------------------

    implicit none

    integer :: i, id, ikeep, shift
    real(dp) :: x, y, z
    real(dp), allocatable, dimension(:, :) :: xgloc, xoloc
    real(dp), allocatable, dimension(:) :: mgloc, moloc
    integer,  allocatable, dimension(:) :: pcounts_o, pcounts_g
    logical :: check

    ng_write = 1
    no_write = 1

    allocate(pcounts_o(0:omp_get_max_threads()))
    pcounts_o = 0
    allocate(pcounts_g(0:omp_get_max_threads()))
    pcounts_g = 0

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP   PRIVATE(i, x, y, z, id, xgloc, xoloc, &
    !$OMP           mgloc, moloc, check, ikeep, shift)
      id = omp_get_thread_num()
      allocate(xgloc(1:ngalaxies, 1:3))
      allocate(mgloc(1:ngalaxies))
      ikeep = 0
      !$OMP DO
        do i = 1, ngalaxies
          x = xg(i, 1) 
          y = xg(i, 2)
          z = xg(i, 3)
          if (moved_x .and. x < 0.5d0) x = x + 1.d0
          if (moved_y .and. y < 0.5d0) y = y + 1.d0
          if (moved_z .and. z < 0.5d0) z = z + 1.d0

          check = (x >= xmin) .and. (x <= xmax)
          check = check .and. (y >= ymin) .and. (y <= ymax)
          check = check .and. (z >= zmin) .and. (z <= zmax)
          check = check .and. (mg(i) >= mthresh)

          if (check) then
            ikeep = ikeep + 1
            xgloc(ikeep,1) = x
            xgloc(ikeep,2) = y
            xgloc(ikeep,3) = z
            mgloc(ikeep) = mg(i)
          endif
        enddo
      !$OMP END DO
      pcounts_g(id+1) = ikeep

      allocate(xoloc(1:norphans, 1:3))
      allocate(moloc(1:norphans))
      ikeep = 0
      !$OMP DO
        do i = 1, norphans
          x = xo(i, 1) 
          y = xo(i, 2)
          z = xo(i, 3)
          if (moved_x .and. x < 0.5d0) x = x + 1.d0
          if (moved_y .and. y < 0.5d0) y = y + 1.d0
          if (moved_z .and. z < 0.5d0) z = z + 1.d0

          check = (x >= xmin) .and. (x <= xmax)
          check = check .and. (y >= ymin) .and. (y <= ymax)
          check = check .and. (z >= zmin) .and. (z <= zmax)
          check = check .and. (mo(i) >= mthresh)

          if (check) then
            ikeep = ikeep + 1
            xoloc(ikeep,1) = x
            xoloc(ikeep,2) = y
            xoloc(ikeep,3) = z
            moloc(ikeep) = mo(i)
          endif
        enddo
      !$OMP END DO
      pcounts_o(id+1) = ikeep

      !$OMP BARRIER
      !$OMP SINGLE
        do i = 0, omp_get_num_threads()-1
          pcounts_g(i+1) = pcounts_g(i+1)+pcounts_g(i)
        enddo
        allocate(xgw(1:pcounts_g(omp_get_num_threads()), 1:3))
        allocate(mgw(1:pcounts_g(omp_get_num_threads())))
        ng_write = pcounts_g(omp_get_num_threads())
      !$OMP END SINGLE
      !$OMP SINGLE
        do i = 0, omp_get_num_threads()-1
          pcounts_o(i+1) = pcounts_o(i+1)+pcounts_o(i)
        enddo
        allocate(xow(1:pcounts_o(omp_get_num_threads()), 1:3))
        allocate(mow(1:pcounts_o(omp_get_num_threads())))
        no_write = pcounts_o(omp_get_num_threads())
      !$OMP END SINGLE
      !$OMP BARRIER


      ! write the particles you want to keep into an array in parallel
      ikeep = pcounts_g(id+1)-pcounts_g(id)
      if (ikeep > 0) then
        shift = pcounts_g(id)
        do i=1, ikeep
          xgw(shift+i,1) = xgloc(i,1)
          xgw(shift+i,2) = xgloc(i,2)
          xgw(shift+i,3) = xgloc(i,3)
          mgw(shift+i) = mgloc(i)
        enddo
      endif

      ikeep = pcounts_o(id+1)-pcounts_o(id)
      if (ikeep > 0) then
        shift = pcounts_o(id)
        do i=1, ikeep
          xow(shift+i,1) = xoloc(i,1)
          xow(shift+i,2) = xoloc(i,2)
          xow(shift+i,3) = xoloc(i,3)
          mow(shift+i) = moloc(i)
        enddo
      endif

      deallocate(xgloc, xoloc)
      deallocate(mgloc, moloc)
    !$OMP END PARALLEL

  write(*, *) "Galaxies kept:", ng_write, "/", ngalaxies
  write(*, *) "Orphans kept:", no_write, "/", norphans
  if (ng_write == 0 .and. no_write==0) then
    write(*, *) "No galaxies and no orphans within boundary and above mass threshold"
    stop
  endif

  end subroutine filter_galaxy_data






  !========================================
  subroutine read_particles()
  !========================================
    !---------------------------------------------
    ! Reads in particle data, file by file, and
    ! only keeps necessary data.
    ! (Particles that are within image borders)
    !---------------------------------------------

    implicit none
    character(len=80)                     :: fname
    character(len=5)                      :: outputnr, cpunr
    integer                               :: i, j, id, ifile, nparts_loc
    integer                               :: ikeep, shift
    real(dp), allocatable, dimension(:,:) :: xp_loc, xkeep
    integer,  allocatable, dimension(:)   :: pcounts
    real(dp) :: x, y, z
    logical :: check

    allocate(pcounts(0:omp_get_max_threads()))
    pcounts = 0

    nparts = 0
    nparts_tot = 0
    outputnr = srcdir(8:12)

    write(*,*) "Reading in particle data. This may take a while."

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(id,i,j,ifile,fname,cpunr,nparts_loc,xp_loc,&
    !$OMP         ikeep,xkeep,shift,x,y,z,check)
      id = omp_get_thread_num()

      allocate(xp_loc(1:pguess, 1:3))
      allocate(xkeep(1:pguess, 1:3))

      ikeep = 0;

      !$OMP DO REDUCTION(+:nparts_tot)
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

          read(600+id) xp_loc(1:nparts_loc,1)
          read(600+id) xp_loc(1:nparts_loc,2)
          read(600+id) xp_loc(1:nparts_loc,3)
          close(600+id)
          ! read(600+id) !vx
          ! read(600+id) !vy
          ! read(600+id) !vz
          ! read(600+id) !mass_loc(1:nparts_loc)
          ! read(600+id) !id_loc(1:nparts_loc)

          do i=1, nparts_loc

            x = xp_loc(i, 1) 
            y = xp_loc(i, 2)
            z = xp_loc(i, 3)
            if (moved_x .and. x < 0.5d0) x = x + 1.d0
            if (moved_y .and. y < 0.5d0) y = y + 1.d0
            if (moved_z .and. z < 0.5d0) z = z + 1.d0

            check = (x >= xmin) .and. (x <= xmax)
            check = check .and. (y >= ymin) .and. (y <= ymax)
            check = check .and. (z >= zmin) .and. (z <= zmax)

            if (check) then
              ikeep = ikeep + 1
              xkeep(ikeep,1) = x
              xkeep(ikeep,2) = y
              xkeep(ikeep,3) = z
            endif
          enddo
          nparts_tot = nparts_tot + nparts_loc
        enddo !loop over files
      !$OMP END DO

      pcounts(id+1) = ikeep
      !$OMP BARRIER
      !$OMP SINGLE
        do i = 0, omp_get_num_threads()-1
          pcounts(i+1) = pcounts(i+1)+pcounts(i)
        enddo
        allocate(xp(1:pcounts(omp_get_num_threads()), 1:3))
      !$OMP END SINGLE
      !$OMP BARRIER


      ! write the particles you want to keep into an array in parallel
      if (ikeep > 0) then
        shift = pcounts(id)
        do i=1, ikeep
          xp(shift+i,1) = xkeep(i,1)
          xp(shift+i,2) = xkeep(i,2)
          xp(shift+i,3) = xkeep(i,3)
        enddo
        !$OMP ATOMIC
        nparts = nparts + ikeep
      endif


      deallocate(xp_loc)
      deallocate(xkeep)
    !$OMP END PARALLEL

    write(*,*) "Nparts kept:", nparts, "out of", nparts_tot
    deallocate(pcounts)

  end subroutine read_particles






  !========================================
  subroutine get_cmdlineargs(halo)
  !========================================
    !-----------------------------------------------------------------  
    ! reads halo ID and output_XXXXX directory from cmdline args.
    !-----------------------------------------------------------------

    implicit none

    integer, intent(out) :: halo
    ! real(dp), intent(out) :: CoM ! in io_module

    character(len=12) :: arg
    character(len=5)  :: dirnr_str
    integer           :: halo_temp, dirnr, i
    logical           :: exists

    do i = 1, iargc()
      call getarg(i, arg)
      if (i==1) then
        read(arg, *) dirnr
        call title(dirnr, dirnr_str)
        srcdir = TRIM('output_'//dirnr_str)
        inquire(file=srcdir, exist=exists)
        if (.not. exists) then
          write(*,*) "Couldn't find directory ", srcdir
          stop
        endif
      else if (i==2) then
        read(arg, "(I12)") halo_temp
        if (halo_temp>=0) halo = halo_temp
      else if (i==3) then
        read(arg, "(F12.8)") CoM(1) 
      else if (i==4) then
        read(arg, "(F12.8)") CoM(2) 
      else if (i==5) then
        read(arg, "(F12.8)") CoM(3) 
      else
        exit
      endif
    end do
  end subroutine get_cmdlineargs






  !========================================
  subroutine write_image()
  !========================================

    !----------------------------------
    ! Write the image data to file.
    !----------------------------------

    implicit none
    character(len=100) :: fname    ! output filename
    character(len=20) :: halo_str ! halo ID as a string
    character(80) :: cmnd

    ! first make subdirectory
    cmnd = 'mkdir -p orphan_plots'
    call system(cmnd)
    cmnd = 'mkdir -p orphan_plots/orphan_plots-'//srcdir(8:12)
    call system(cmnd)
    
    call get_halo_string(halo_str, halo)
    fname = TRIM('orphan_plots/orphan_plots-'//&
                srcdir(8:12)//'/orphan_plot-halo-'//&
                TRIM(halo_str)//'.dat')

    if (.not.infofile_read) then
      write(*, *) "Missing data from info file to dump density projections :/"
      stop
    endif

    open(unit=667, file=fname, form='unformatted')
    write(667) nc
    write(667) aexp
    write(667) unit_l
    write(667) xmax
    write(667) xmin
    write(667) ymax
    write(667) ymin
    write(667) zmax
    write(667) zmin
    write(667) density_field_xy
    write(667) density_field_yz
    write(667) density_field_xz
    write(667) ng_write
    write(667) xgw(:, 1)
    write(667) xgw(:, 2)
    write(667) xgw(:, 3)
    write(667) mgw(:)
    write(667) no_write
    write(667) xow(:, 1)
    write(667) xow(:, 2)
    write(667) xow(:, 3)
    write(667) mow(:)
    close(667)

  end subroutine write_image

end program
