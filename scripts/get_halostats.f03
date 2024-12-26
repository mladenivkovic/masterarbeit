!===========================================
! Get an image of a halo of your choosing
! and write out galaxies and orphans that
! are inside the halo
!===========================================


!===================================
program get_halostats
!===================================

  !==============================================
  ! cmd line args:
  ! get_halostats.o outputnr halo
  ! outputnr: 67 for output_00067
  ! halo:     halo ID
  !==============================================

  use omp_lib
  use constants_and_parameters
  use io_module
  use density_projection
  implicit none

  !-------------------------------------------
  ! Global variables and parameters
  !-------------------------------------------

  ! Manual parameters
  integer, parameter  :: nsamples            = 200    ! number of samples for output; Number of bins for histograms

  logical :: halo_density_projection_plots = .true. ! make plots
  ! logical :: halo_density_projection_plots = .false.

  ! logical :: print_consistency_checks = .true. ! print some diagnostics
  logical :: print_consistency_checks = .false. ! print some diagnostics

  character(len=7) :: which_CoM = "pcenter" ! use partcenter.f90 output. See run_partcenter.sh script for output
  ! character(len=7) :: which_CoM = "peakpos" ! use halo finder's peak position as center of mass
  ! character(len=7) :: which_CoM = "cgalaxy" ! use central galaxy as center of mass
  ! character(len=7) :: which_CoM = "compute" ! compute center of mass based on particles



  !------------------
  ! Read-in data
  !------------------
  integer             :: nchildren, nparts, nparts_tot, nparts_tot_expect
  integer             :: ng_halo, no_halo              ! number of galaxies/orphans in halo
  integer             :: halo                          ! the halo to work with


  !-------------------------
  ! Arrays
  !-------------------------
  ! read-in data arrays
  integer                               :: pguess = 5*2000000 ! initial guess for number of particles. Take npartmax!
  integer                               :: childguess = 10000 ! initial guess for number of children
  real(dp), dimension(:,:), allocatable :: xow, xgw           ! positions of orphans/galaxies in halo  (w: write)
  real(dp), dimension(:), allocatable   :: mow, mgw           ! masses of orphans/galaxies in halo  (w: write)
  integer,  dimension(:), allocatable   :: children           ! list of children ID's

  real(dp), dimension(:,:), allocatable :: xp     ! positions of particles [code units]
  real(dp), dimension(:), allocatable   :: mp     ! masses of particles [code units]
  integer , dimension(:), allocatable   :: idp    ! id of particles


  !-----------------------------
  ! Other variables
  !-----------------------------
  real(dp)  :: R200                          ! R200 wrt critical density
  real(dp)  :: M200                          ! M200 wrt critical density
  real(dp)  :: mpart                         ! mass of single particle
  real(dp)  :: rmax                          ! Particle furthest away from center
  real(dp)  :: partmass, partmass_tot        ! read-in particle masses for consistency checks
  real(dp)  :: rho_crit                      ! critical density

  logical   :: moved_x, moved_y, moved_z ! whether x,y, or z have been moved for periodic boundary corrections
  

  character(len=80) :: fname    ! output filename
  character(len=20) :: halo_str ! halo ID as a string

  !---------------------------
  ! Set-up
  !---------------------------

  ! imported from density_projection module
  ! interpolation = 'ngp'    ! which interpolation to use: cic for cloud in cell, ngp for nearest grid point
  interpolation = 'cic'   ! which interpolation to use: cic for cloud in cell, ngp for nearest grid point
  nc = 1000                ! number of cells per dimension for image


  call get_srcdir_and_halo_from_cmdlineargs(halo) ! from io_module
  call read_info() ! from io_module
  call read_clump_data(halo)
  call print_infos()

  call find_all_children()
  call read_particles()
  call check_periodicity_particles()
  call compute_derived_quantities_and_add_units()

  if (halo_density_projection_plots)  then
    call allocate_density_fields()
    allocate(mp(1:nparts))
    mp = 1.d0 / (2**(levelmin)) ** 3
    call get_density_projection(xp, mp, nparts)
    call get_halo_string(halo_str, halo)
    fname = TRIM(srcdir//'/density_image_new-halo-'//TRIM(halo_str)//'.dat')
    call write_density_projections(fname)
    call deallocate_density_fields()
  endif

  call read_galaxy_data()
  call check_periodicity_galaxies()
  call write_galaxy_data()
  call get_center()
  call get_R200()
  call write_halo_data()
  ! call get_radial_profile() ! deprecated, python does this now



  write(*,*) "get_halostats.f03 finished."

  
  call deallocate_io_module()
  deallocate(xp, mp)
  deallocate(xow, xgw, mow, mgw)



contains




  !=======================================
  subroutine check_periodicity_galaxies()
  !=======================================
    !-----------------------------------------------
    ! Move galaxies the same way you moved particles
    ! with check_periodicity_particles()
    !-----------------------------------------------
    
    implicit none
    integer :: i

    if (moved_x .or. moved_y .or. moved_z) then
      write(*,*) "Correcting galaxy positions for periodicity."
    else
      return
    endif

    !$OMP PARALLEL PRIVATE(i)
      if (moved_x) then
        !$OMP DO 
          do i=1, ngalaxies
            if (xg(i,1)<=0.5d0) xg(i,1) = xg(i,1)+1.d0
          enddo
        !$OMP END DO NOWAIT
        !$OMP DO 
          do i=1, norphans
            if (xo(i,1)<=0.5d0) xo(i,1) = xo(i,1)+1.d0
          enddo
        !$OMP END DO NOWAIT
      endif
      if (moved_y) then
        !$OMP DO 
          do i=1, ngalaxies
            if (xg(i,2)<=0.5d0) xg(i,2) = xg(i,2)+1.d0
          enddo
        !$OMP END DO NOWAIT
        !$OMP DO
          do i=1, norphans
            if (xo(i,2)<=0.5d0) xo(i,2) = xo(i,2)+1.d0
          enddo
        !$OMP END DO NOWAIT
      endif
      if (moved_z) then
        !$OMP DO 
          do i=1, ngalaxies
            if (xg(i,3)<=0.5d0) xg(i,3) = xg(i,3)+1.d0
          enddo
        !$OMP END DO NOWAIT
        !$OMP DO 
          do i=1, norphans
            if (xo(i,3)<=0.5d0) xo(i,3) = xo(i,3)+1.d0
          enddo
        !$OMP END DO
      endif

    !$OMP END PARALLEL

  end subroutine check_periodicity_galaxies






  !=======================================
  subroutine check_periodicity_particles()
  !=======================================
    !-----------------------------------------------
    ! Check if clumps are cut in half by boundary,
    ! and correct if that's the case.
    !-----------------------------------------------
    
    implicit none
    integer :: i

    write(*,*) "Checking for periodicity corrections."

    xmax = -1
    xmin = 2
    ymax = -1
    ymin = 2
    zmax = -1
    zmin = 2
    moved_x = .false.
    moved_y = .false.
    moved_z = .false.

    !$OMP PARALLEL
      !$OMP DO REDUCTION(max:xmax,ymax,zmax) REDUCTION(min:xmin,ymin,zmin)
        do i = 1, nparts
          xmax = max(xmax, xp(i,1))
          xmin = min(xmin, xp(i,1))
          ymax = max(ymax, xp(i,2))
          ymin = min(ymin, xp(i,2))
          zmax = max(zmax, xp(i,3))
          zmin = min(zmin, xp(i,3))
        enddo
      !$OMP END DO

      if (xmax-xmin>0.5d0) then
        !$OMP SINGLE
          write(*,*) "Need to move some x-values", xmin, xmax
          moved_x = .true.
        !$OMP END SINGLE
        !$OMP SINGLE
          xmax = -1
          xmin = 2
        !$OMP END SINGLE
        !$OMP BARRIER
        !$OMP DO REDUCTION(max:xmax), REDUCTION(min:xmin)
          do i=1, nparts
            if (xp(i,1)<=0.5d0) xp(i,1) = xp(i,1)+1.d0
            xmax = max(xmax, xp(i,1))
            xmin = min(xmin, xp(i,1))
          enddo
        !$OMP END DO
      endif
      if (ymax-ymin>0.5d0) then
        !$OMP SINGLE
          write(*,*) "Need to move some y-values", ymin, ymax
          moved_y=.true.
        !$OMP END SINGLE
        !$OMP SINGLE
          ymax = -1
          ymin = 2
        !$OMP END SINGLE
        !$OMP BARRIER
        !$OMP DO REDUCTION(max:ymax), REDUCTION(min:ymin)
          do i=1, nparts
            if (xp(i,2)<=0.5d0) xp(i,2) = xp(i,2)+1.d0
            ymax = max(ymax, xp(i,2))
            ymin = min(ymin, xp(i,2))
          enddo
        !$OMP END DO
      endif
      if (zmax-zmin>0.5d0) then
        !$OMP SINGLE
          write(*,*) "Need to move some z-values", zmin, zmax
          moved_z=.true.
        !$OMP END SINGLE
        !$OMP SINGLE
          zmax = -1
          zmin = 2
        !$OMP END SINGLE
        !$OMP BARRIER
        !$OMP DO REDUCTION(max:zmax), REDUCTION(min:zmin)
          do i=1, nparts
            if (xp(i,3)<=0.5d0) xp(i,3) = xp(i,3)+1.d0
            zmax = max(zmax, xp(i,3))
            zmin = min(zmin, xp(i,3))
          enddo
        !$OMP END DO
      endif

    !$OMP END PARALLEL

  end subroutine check_periodicity_particles




  !=======================================
  subroutine find_all_children()
  !=======================================
    use quick_sort

    !------------------------------------------------
    ! Find all children of given halo recusively.
    !------------------------------------------------

    implicit none
    integer      :: child_start, child_end, i, j
    logical      :: found_new
    ! integer, allocatable, dimension(:) :: order

    ! first find all direct children

    nchildren = 0
    allocate(children(1:childguess))
    children = 0

    !$OMP PARALLEL PRIVATE(i,j) DEFAULT(shared)
      !$OMP DO
        do i = 1, nclumps
          if (parent_id(i)==halo .and. clmp_id(i)/=parent_id(i)) then
            !$OMP CRITICAL
              nchildren = nchildren + 1
              children(nchildren) = clmp_id(i)
            !$OMP END CRITICAL
          endif
        enddo
      !$OMP ENDDO

      !$OMP SINGLE
        child_start = 1
        child_end = nchildren
        found_new = .true.
      !$OMP END SINGLE
      !$OMP BARRIER

      ! Now look for children's children
      do while (found_new) 
        !$OMP DO
          do i=1, nclumps
            do j=child_start, child_end
              if(children(j)==parent_id(i)) then
                if (clmp_id(i)/=parent_id(i)) then
                  !$OMP CRITICAL
                    nchildren = nchildren + 1
                    children(nchildren) = clmp_id(i)
                  !$OMP END CRITICAL
                  exit
                endif
              endif
            enddo
          enddo
        !$OMP ENDDO
        !$OMP SINGLE
          found_new = .false.
          if (child_end<nchildren) then
            found_new=.true.
            child_start=child_end+1
            child_end=nchildren
          endif
        !$OMP END SINGLE
        !$OMP BARRIER
      enddo
    !$OMP END PARALLEL

    ! for checks with other scripts
    ! allocate(order(1:nchildren))
    ! order = 0
    ! call quick_sort_int(children, order, nchildren)
    !
    ! do i = 1, nchildren
    !   write(*, '(I12)') children(i)
    ! enddo
    ! write(*, '(A5,x,I12)') "Halo", halo
    ! write(*, '(A5,x,I12)') "Total", nchildren
    ! stop

    ! Add halo to list of children for simplicity
    nchildren = nchildren + 1
    children(nchildren) = halo

    write(*,*) "Found", nchildren-1, "child clumps out of ", nclumps, "total clumps in snapshot."

  end subroutine find_all_children



  !========================================
  subroutine get_center()
  !========================================

    !------------------------------------------------
    ! Determine which center to use, and do 
    ! additional work if necessary.
    ! Needs to be called AFTER particle periodicity
    ! checks, and AFTER galaxy periodicity checks
    !------------------------------------------------

    implicit none
    real(dp)    :: xc, yc, zc, rmax, dx, dy, dz, dc
    integer     :: i
    character(len=80) :: partcenterfile
    character(len=20) :: halo_str
    logical :: file_exists


    ! central galaxy is set at index 0 of xgw array
    if (which_CoM == "cgalaxy") then
      CoM(1) = xgw(0,1)
      CoM(2) = xgw(0,2)
      CoM(3) = xgw(0,3)

    else if (which_CoM == "peakpos") then
      ! if you use peak position as center of mass, apply periodicity correction
      ! if necessary
      CoM(1) = peakpos(1)
      CoM(2) = peakpos(2)
      CoM(3) = peakpos(3)
      write(*, *) "Checking CoM periodicity"
      if (moved_x) then
        if (CoM(1) <= 0.5) CoM(1) = CoM(1) + 1.d0
        write(*,*) "Correcting CoM for x"
      endif
      if (moved_y) then
        if (CoM(2) <= 0.5) CoM(2) = CoM(2) + 1.d0
        write(*,*) "Correcting CoM for y"
      endif
      if (moved_z) then
        if (CoM(3) <= 0.5) CoM(3) = CoM(3) + 1.d0
        write(*,*) "Correcting CoM for z"
      endif

    else if (which_CoM == "compute") then
      xc = 0; yc = 0; zc = 0;
      !$OMP PARALLEL PRIVATE(i) DEFAULT(shared)
        ! Compute center of mass
        !$OMP DO REDUCTION(+:xc, yc, zc)
          do i = 1, nparts
            xc = xc + xp(i,1)
            yc = yc + xp(i,2)
            zc = zc + xp(i,3)
          enddo
        !$OMP ENDDO
      !$OMP END PARALLEL
      CoM(1) = xc/nparts
      CoM(2) = yc/nparts
      CoM(3) = zc/nparts


    else if (which_CoM == "pcenter") then
      write(*, *) "Using partcenter results!"

      call get_halo_string(halo_str, halo)
      partcenterfile=TRIM(srcdir//'/partcenter-output-halo-'//TRIM(halo_str)//'.dat')

      inquire(file=partcenterfile, exist=file_exists)
      if (.not.file_exists) then
        write(*, *) "Couldn't find file", partcenterfile
        stop
      endif

      open(unit=666, file=partcenterfile)
      ! leave 1PE here, that's how partcenter.f90 does it
      read(666,'(4(1X, 1PE14.7))') rmax, CoM(1), CoM(2), CoM(3)
      close(666)

      write(*, *) "Checking CoM periodicity"
      if (moved_x) then
        if (CoM(1) <= 0.5) CoM(1) = CoM(1) + 1.d0
        write(*,*) "Correcting CoM for x"
      endif
      if (moved_y) then
        if (CoM(2) <= 0.5) CoM(2) = CoM(2) + 1.d0
        write(*,*) "Correcting CoM for y"
      endif
      if (moved_z) then
        if (CoM(3) <= 0.5) CoM(3) = CoM(3) + 1.d0
        write(*,*) "Correcting CoM for z"
      endif

      dx = CoM(1) - peakpos(1)
      dy = CoM(2) - peakpos(2)
      dz = CoM(3) - peakpos(3)
      dc = sqrt(dx**2 + dy**2 + dz**2)
      ! write(*, '(A10, 4(x, E14.7))') "Diff:" dc
      ! write(*, '(A,4(x,E14.7))') " Diff Partcenter - Peakpos [Mpc]:", &
      !                               dc * unit_l, dx*unit_l, dy*unit_l, dz*unit_l
      !
    end if
  
  end subroutine get_center






  !=================================
  subroutine get_R200()
  !=================================

    !------------------------------------------------
    ! Computes the Center of Mass of the particles
    ! and the R200 with the CoM as the origin
    !------------------------------------------------

    use quick_sort

    implicit none
    integer :: i, j, ikeep
    real(dp), allocatable, dimension(:) :: r
    integer, allocatable, dimension(:)  :: order_r
    real(dp)    :: f

    real(dp) :: M200_inner, M200_outer, R200_inner, R200_outer


    allocate(r(1:nparts))
    allocate(order_r(1:nparts))

    !$OMP PARALLEL PRIVATE(i, j) DEFAULT(shared)
      ! compute distances from center of mass
      !$OMP DO
        do i = 1, nparts
          r(i) = 0.d0
          do j = 1, 3
            r(i) = r(i) + (xp(i,j) - CoM(j))**2
          enddo
          r(i) = sqrt(r(i))
          order_r(i) = i
        enddo
      !$OMP ENDDO
    !$OMP END PARALLEL

    call quick_sort_real(r, order_r, nparts)

    write(*, '(A, 4(x,E12.3))') "Particle at Rmax", xp(order_r(nparts), :), &
        sqrt((xp(order_r(nparts),1) - CoM(1))**2 + (xp(order_r(nparts),2) - CoM(2))**2 + (xp(order_r(nparts),3) - CoM(3))**2) 

    ! Sorting checks
    ! do i = 1, nparts - 1
    !   if (r(i) > r(i+1)) then
    !     write(*, *) "SOMETHING'S WRONG WITH THE SORTING!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    !   endif
    ! enddo
 

    rmax = r(nparts)

    ! Get M200 and R200
    ! Start from the outside, and work your way towards the center
    f = 4./3. * pi * 200 * rho_crit ! this is in M_Sol / Mpc^3
    f = f * M_Sol / Mpc**3

    M200_outer = nparts * mpart 
    ikeep = nparts
    do  i = nparts, 1, -1
      ! write(*, "(I10, 2E18.6, F12.6, E18.6)") i, M200, f*r(i)**3, r(i)/rmax, r(i)
      ! if (M200_outer * unit_l**3 / M_Sol >= f * r(i)**3 * unit_l**3 / Mpc**3) then
      if (M200_outer >= f * r(i)**3) then
        ikeep = i
        exit
      endif
      M200_outer = M200_outer - mpart
    end do
    R200_outer = r(ikeep)

    M200_inner = 100*mpart
    ikeep = 100
    do  i = ikeep+1, nparts
      if (M200_inner < f * r(i)**3) then
        ikeep = i
        exit
      endif
      M200_inner = M200_inner + mpart
    end do
    R200_inner = r(ikeep)

    write(*, '(A, 2(x,E12.3))') "M200 inner/outer", &
              M200_inner * unit_l**3 / M_Sol , M200_outer * unit_l**3 / M_Sol 
    write(*, '(A, 2(x,E12.3))') "R200 inner/outer", &
              R200_inner * unit_l / Mpc, R200_outer * unit_l / Mpc

    if (M200_inner > M200_outer) then
      M200 = M200_inner
      R200 = R200_inner
    else
      M200 = M200_outer
      R200 = R200_outer
    endif

    write(*,'(A,2(x,E14.6), A)') "M200:", M200, M200*unit_l**3/M_Sol, " [MSol]"
    write(*,'(A,2(x,E14.6), A)') "R200:", R200, R200*unit_l/Mpc, " [Mpc]"
    write(*,'(A,2(x,E14.6), A)') "Rmax:", rmax, rmax*unit_l/Mpc, " [Mpc]"

    deallocate(r, order_r)

  end subroutine get_R200





  !====================================
  subroutine print_infos()
  !====================================

    implicit none
    write(*,*) "============================"
    write(*,*) "get_halostats.f03 started."
    write(*,*) "============================"
    write(*,*)                "Working with parameters:"
    write(*,'(A20,x,A25)')      "  srcdir:        ", srcdir
    write(*,'(A20,x,I25)')      "  halo:          ", halo
    write(*,'(A20,x,I25)')      "  ncpu:          ", ncpu
    write(*,'(A20,x,I25)')      "  nc:            ", nc
    write(*,'(A20,x,I25)')      "  nclumps:       ", nclumps
    ! write(*,'(A20,x,I25)')      "  ngalaxies:     ", ngalaxies
    write(*,'(A20,x,A25)')      "  interpolation: ", interpolation
    write(*,'(A20,x,I25)')      "  levelmin:      ", levelmin
    write(*,'(A20,x,I25)')      "  levelmax:      ", levelmax

    write(*,'(A20,x,E25.3,x,F25.3,x,A9)') "  unit_l:        ", unit_l, (unit_l/Mpc), "[Mpc]"
    write(*,'(A20,x,F25.3)')    "  a_exp:         ", aexp
    write(*,'(A20,x,F25.3)')    "  z:             ", 1.d0/aexp - 1.d0
    write(*,'(A20,x,F25.3)')    "  H0:            ", H0
    write(*,'(A20,x,F25.3)')    "  x_peak:        ", peakpos(1)
    write(*,'(A20,x,F25.3)')    "  y_peak:        ", peakpos(2)
    write(*,'(A20,x,F25.3)')    "  z_peak:        ", peakpos(3)

  end subroutine print_infos






  !========================================
  subroutine read_particles()
  !========================================
    !---------------------------------------------
    ! Reads in particle data, file by file, and
    ! only keeps necessary data.
    ! (Particles that are in halo or its children)
    !---------------------------------------------

    implicit none
    character(len=80)                     :: fname
    character(len=5)                      :: outputnr, cpunr
    integer                               :: i, j, id, ifile, nparts_loc
    integer                               :: ikeep, shift, ikeep_old
    integer                               :: temp
    real(dp), allocatable, dimension(:,:) :: xp_loc, xkeep
    real(dp), allocatable, dimension(:)   :: mass_loc
    integer,  allocatable, dimension(:)   :: pclmpid_loc, pcounts
    integer(i8), allocatable, dimension(:) :: id_loc, idkeep
    real(dp) :: partmass_read, partmass_tot_read

    allocate(pcounts(0:omp_get_max_threads()))
    pcounts = 0

    nparts = 0
    nparts_tot = 0
    outputnr = srcdir(8:12)

    write(*,*) "Reading in particle data. This may take a while."

    !$OMP PARALLEL DEFAULT(SHARED) &
    !$OMP PRIVATE(id,i,j,ifile,fname,cpunr, &
    !$OMP         nparts_loc,xp_loc,id_loc,pclmpid_loc, &
    !$OMP         ikeep,xkeep,idkeep,temp,ikeep_old,shift,mass_loc) 
      id = omp_get_thread_num()

      allocate(xp_loc(1:pguess, 1:3))
      allocate(id_loc(1:pguess))
      allocate(pclmpid_loc(1:pguess))
      allocate(mass_loc(1:pguess))

      allocate(xkeep(1:pguess, 1:3))
      allocate(idkeep(1:pguess))

      ikeep = 0; ikeep_old = 0;

      !$OMP DO REDUCTION(+:nparts_tot, partmass_tot_read)
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
          read(600+id) !vx
          read(600+id) !vy
          read(600+id) !vz
          read(600+id) mass_loc(1:nparts_loc)
          read(600+id) id_loc(1:nparts_loc)
          close(600+id)

          if (ifile == 1) then
            partmass_read = mass_loc(1)
          endif
          partmass_tot_read = partmass_tot_read + sum(mass_loc(1:nparts_loc))

          fname=TRIM(srcdir//'/unbinding_'//outputnr//'.out'//cpunr)
          open(unit=600+id, file=fname, form='unformatted')
          read(600+id) pclmpid_loc(1:nparts_loc)
          close(600+id)

          ikeep_old = ikeep

          do i=1, nparts_loc
            temp = abs(pclmpid_loc(i))
            do j=1, nchildren
              if (children(j)==temp) then
                ikeep = ikeep + 1
                xkeep(ikeep,1) = xp_loc(i,1)
                xkeep(ikeep,2) = xp_loc(i,2)
                xkeep(ikeep,3) = xp_loc(i,3)
                idkeep(ikeep)  = id_loc(i)
                exit
              endif
            enddo
          enddo
          nparts_tot = nparts_tot + nparts_loc
          ! write(*,'(A2,x,I5,x,A34,x,I10,x,A6,x,I10,x,A10)') &
            ! "ID", id, "read part_"//outputnr//'.out'//cpunr//'; Keeping',&
            ! ikeep-ikeep_old, "out of", nparts_loc, "particles."
        enddo !loop over files
      !$OMP END DO

      pcounts(id+1) = ikeep
      !$OMP BARRIER
      !$OMP SINGLE
        do i = 0, omp_get_num_threads()-1
          pcounts(i+1) = pcounts(i+1)+pcounts(i)
        enddo
        allocate(xp(1:pcounts(omp_get_num_threads()), 1:3))
        allocate(idp(1:pcounts(omp_get_num_threads())))
      !$OMP END SINGLE
      !$OMP BARRIER


      ! write the particles you want to keep into an array in parallel
      if (ikeep > 0) then
        shift = pcounts(id)
        do i=1, ikeep
          xp(shift+i,1) = xkeep(i,1)
          xp(shift+i,2) = xkeep(i,2)
          xp(shift+i,3) = xkeep(i,3)
          idp(shift+i)  = idkeep(i)
        enddo
        !$OMP ATOMIC
        nparts = nparts + ikeep
      endif


      deallocate(xp_loc, id_loc, pclmpid_loc, mass_loc)
      deallocate(xkeep, idkeep)
    !$OMP END PARALLEL

    write(*,*) "Nparts kept:", nparts, "out of", nparts_tot
    partmass = partmass_read
    partmass_tot = partmass_tot_read
    deallocate(pcounts)

  end subroutine read_particles




  !=======================================================
  subroutine compute_derived_quantities_and_add_units()
  !=======================================================

    !--------------------------------------------------------
    ! After reading in particle data and info data, compute
    ! derived quantities and check them for consistency
    ! Also transform from code units into physical units.
    !--------------------------------------------------------

    implicit none

    real(dp) :: t1, t2, t3, t4, t5, t6, m0, m1

    ! unit_m = unit_rho * unit_l**3 / M_Sol

    nparts_tot_expect = (2**levelmin)**3

    rho_crit = 3.d0 * H0**2 / 8.d0 / pi / G
    mpart = unit_rho / nparts_tot_expect ! per unit_l^3

    if (print_consistency_checks) then
      write(*, '(A)')  "--- Consistency Checks"

      write(*, '(A,2(x,E18.6E2),E12.3,A)') &
            "--- rho_crit", rho_crit, &
            2.7754d11*(H0/100)**2, &
            rho_crit /( 2.7754d11*(H0/100)**2), "/1.0"

      write(*, '(A,2(I16),2(I4,A2))') &
            "--- Total Particles", nparts_tot, &
            nparts_tot_expect, nparts_tot / nparts_tot_expect, "/1", &
            modulo(nparts_tot, nparts_tot_expect), "/0"

      m0 = mpart ! unit density at this point
      m1 = mpart*unit_l**3/M_Sol
      write(*, '(A,2(x,E12.6),A)') &
            "--- Particle mass Version 1: ", m0, m1, " [M_Sol]"

      ! rho_crit is in units of M_Sol/Mpc^3 already!
      t1 = Omega_m*rho_crit/Mpc**3*M_Sol / dble(nparts_tot_expect)
      t2 = Omega_m*rho_crit*(unit_l/Mpc)**3 / dble(nparts_tot_expect)
      write(*, '(A,2(x,E12.6),A)') &
            "--- Particle mass Version 2: ", t1, t2, " [M_Sol]"

      ! rho_crit is in units of M_Sol/Mpc^3 already!
      ! hardcoded omega_b used for ICs
      t3 = (Omega_m + 0.0486d0)*rho_crit/Mpc**3*M_Sol / dble(nparts_tot_expect)
      t4 = (Omega_m + 0.0486d0)*rho_crit*(unit_l/Mpc)**3 / dble(nparts_tot_expect)
      write(*, '(A,2(x,E12.6),A)') &
            "--- Particle mass Version 3: ", t3, t4, " [M_Sol]"

      ! rho_crit is in units of M_Sol/Mpc^3 already!
      ! hardcoded omega_b used for ICs
      t5 = (Omega_m - 0.0486d0)*rho_crit/Mpc**3*M_Sol / dble(nparts_tot_expect)
      t6 = (Omega_m - 0.0486d0)*rho_crit*(unit_l/Mpc)**3 / dble(nparts_tot_expect)
      write(*, '(A,2(x,E12.6),A)') &
            "--- Particle mass Version 4: ", t5, t6, " [M_Sol]"


      write(*, '(A,4(x,E12.4))') &
            "--- Ratios: ", m0/t1, m1/t2, m0/t3, m1/t4


      write(*,'(A, 4(x,E12.3),A)') &
            "--- Particle mass array check 1:", &
            partmass, &
            partmass_tot, &
            partmass_tot / dble(nparts_tot), &
            partmass_tot / dble(nparts_tot) / partmass, "/1.0"

      write(*,'(A, 3(x,E12.4))') &
            "--- Particle mass array check 2:", &
            partmass*unit_rho*unit_l**3, &
            partmass_tot*unit_rho*unit_l**3, &
            partmass_tot*unit_rho*unit_l**3 / (rho_crit * Omega_m * unit_l**3)

      write(*, '(A,4(x,E12.4))') "--- Redshift Checks", aexp, aexp**3, 1d0/(aexp), (1d0/aexp)**3
    endif

  end subroutine compute_derived_quantities_and_add_units










  !===============================================
  subroutine write_galaxy_data()
  !===============================================
    !--------------------------------------------------------
    ! Find out which galaxies and oprhans should be included 
    ! and write them to file.
    !--------------------------------------------------------

    use quick_sort

    implicit none

    integer :: i, ip, io, ig, ic
    integer,  dimension(:), allocatable  :: order_p
    integer,  dimension(:), allocatable  :: order_o
    integer,  dimension(:), allocatable  :: order_g
    integer,  dimension(:), allocatable  :: order_c
    real(dp) :: l
    character(len=80) :: fname
    character(len=20) :: halo_str



    !------------------------------
    ! Sort out galaxies
    !------------------------------

    call get_halo_string(halo_str, halo)

    write(*,*) "Sorting out Clump Galaxies"


    allocate(order_g(1:ngalaxies))
    allocate(order_c(1:nchildren))
    allocate(xgw(0:ngalaxies,1:3)) !store central galaxy at index 0
    allocate(mgw(0:ngalaxies)) !store central galaxy at index 0


    order_g = [(i, i=1, ngalaxies)]
    order_c = [(i, i=1, nchildren)]
    call quick_sort_int(children, order_c, nchildren)
    call quick_sort_int(gal_id, order_g, ngalaxies)


    ig = 1
    ic = 1
    ng_halo = 0
    do
      if (gal_id(ig)<children(ic)) then
        ig = ig + 1
      else if (gal_id(ig)>children(ic)) then
        ic = ic + 1
      else !found match!
        if (gal_id(ig)==halo) then !central
          xgw(0,1) = xg(order_g(ig),1)
          xgw(0,2) = xg(order_g(ig),2)
          xgw(0,3) = xg(order_g(ig),3)
          mgw(0) = mg(order_g(ig))
        else
          ng_halo = ng_halo + 1
          xgw(ng_halo,1) = xg(order_g(ig),1)
          xgw(ng_halo,2) = xg(order_g(ig),2)
          xgw(ng_halo,3) = xg(order_g(ig),3)
          mgw(ng_halo) = mg(order_g(ig))
        endif
        ig = ig + 1
        ic = ic + 1
      endif
      if (ig>ngalaxies) exit
      if (ic>nchildren) exit
    enddo

    write(*,'(A6,I10,A45,I10,A25)') "Found", ng_halo, "satellite galaxies with clumps out of", ngalaxies, "galaxies in snapshot."

    l = unit_l / Mpc

    fname = TRIM(srcdir//'/profile-galaxies-'//TRIM(halo_str)//'.dat')
    open(unit=666, file=fname, form='formatted')
    write(666, '(A,A14,A15,A15,A15)') "#","   m","   x", "   y", "   z"
    ! write central first
    do ig = 0, ng_halo
      write(666,'(4(E15.6))') mgw(ig), xgw(ig,1)*l, xgw(ig,2)*l, xgw(ig,3)*l
    enddo
    close(666)
  
    deallocate(order_c, order_g)





    !------------------------------
    ! Sort out orphans
    !------------------------------

    write(*,*) "Sorting out Orphan Galaxies"

    allocate(order_p(1:nparts))
    allocate(order_o(1:norphans))
    allocate(xow(1:norphans,1:3))
    allocate(mow(1:norphans))

    order_p = [(i, i=1, nparts)]
    order_o = [(i, i=1, norphans)]

    call quick_sort_int(idp, order_p, nparts)
    call quick_sort_int(ido, order_o, norphans)


    ip = 1
    io = 1
    no_halo = 0
    do
      if (idp(ip)<ido(io)) then
        ip = ip + 1
      else if (idp(ip)>ido(io)) then
        io = io + 1
      else !found match!
        no_halo = no_halo + 1
        xow(no_halo,1) = xo(order_o(io),1)
        xow(no_halo,2) = xo(order_o(io),2)
        xow(no_halo,3) = xo(order_o(io),3)
        mow(no_halo) = mo(order_o(io))
        ip = ip + 1
        io = io + 1
      endif
      if (ip>nparts) exit
      if (io>norphans) exit
    enddo 


    fname = TRIM(srcdir//'/profile-orphans-'//TRIM(halo_str)//'.dat')
    open(unit=666, file=fname, form='formatted')
    write(666, '(A,A14,A15,A15,A15)') "#","  m","   x", "   y", "   z"
    do ig = 1, no_halo
      write(666,'(4(E15.6))') mow(ig), xow(ig,1)*l, xow(ig,2)*l, xow(ig,3)*l
    enddo
    close(666)


    write(*,'(A6,I10,A45,I10,A25)') "Found", no_halo, &
                                    "orphan galaxies in halo out of", norphans, &
                                    "orphans in snapshot."



    deallocate(order_o, order_p)
    deallocate(children)
    deallocate(idp)

  end subroutine write_galaxy_data





  !=======================================
  subroutine write_halo_data()
  !=======================================

    !------------------------------------------------
    ! Writes halo metadata, like center of mass and
    ! R200, M200....
    ! Writes it in outut_XXXX/radial-profile-metadata
    ! file
    !------------------------------------------------

    implicit none

    character(len=80) :: fname
    character(len=20) :: halo_str


    call get_halo_string(halo_str, halo)
    if (which_CoM == 'pcenter') then
      fname = TRIM(srcdir//'/radial-profile-partcenter-metadata-'//TRIM(halo_str)//'.dat')
    else
      fname = TRIM(srcdir//'/radial-profile-metadata-'//TRIM(halo_str)//'.dat')
    endif

    open(unit=666, file=fname, form='formatted')
    write(666, "(A)") "# center of mass [internal units]"
    write(666, "(E14.6)") CoM(1)
    write(666, "(E14.6)") CoM(2)
    write(666, "(E14.6)") CoM(3)
    write(666, "(A)") "# R200 [internal units]"
    write(666, "(E14.6)") R200
    write(666, "(A)") "# Rmax [internal units]"
    write(666, "(E14.6)") rmax
    write(666, "(A)") "# M200 [Msol]"
    write(666, "(E14.6)") M200 * unit_l**3 / M_Sol 
    ! TODO: check whether sum of whole array is correct
    write(666, "(A)") "# Unit_l"
    write(666, "(E14.6)") unit_l
    write(666, "(A)") "# Unit_l [Mpc]"
    write(666, "(E14.6)") unit_l / Mpc

    close(666)

  end subroutine write_halo_data



    ! !=================================
    ! subroutine get_radial_profile()
    ! !=================================
    !
    !   DEPRECATED, PYTHON DOES THIS NOW
    !
    !   !--------------------------------------------
    !   ! Computes and writes to file the cumulative
    !   ! radial number density of galaxies and
    !   ! orphans
    !   !--------------------------------------------
    !
    !   implicit none
    !   integer  :: i, ind
    !   real(dp) :: d
    !   real(dp), dimension(0:nsamples) :: prof_ng, prof_no
    !   real(dp), dimension(0:nsamples) :: prof_mg, prof_mo
    !   ! character(len=80) :: fname
    !   ! character(len=20) :: halo_str
    !   real(dp) :: dr
    !
    !   ! rmax = max(xgw(0)-xmin, xmax-xgw(0), ygw(0)-ymin, ymax-ygw(0))
    !   ! find and overestimated maximal radius for bin edgte
    !   ! rmax = sqrt(4*(xmax-xmin)**2 + (zmax-zmin)**2)*unit_l ! xmax=ymax, xmin=ymin for box
    !   ! write(*,*) "Rmax Radial Profile:", rmax, xgw(0,1), xgw(0,2)
    !
    !   prof_mg = 0
    !   prof_mo = 0
    !   prof_ng = 0
    !   prof_no = 0
    !
    !   !$OMP PARALLEL PRIVATE(i,d,ind)
    !
    !     !$OMP DO
    !       do i = 1, ng_halo
    !         d = sqrt((xgw(i,1)-CoM(1))**2+(xgw(i,2)-CoM(2))**2+(xgw(i,3)-CoM(3))**2) * unit_l
    !         ind = int(d/rmax*nsamples)+1
    !         !$OMP ATOMIC
    !         prof_ng(ind) = prof_ng(ind)+1
    !         prof_mg(ind) = prof_mg(ind)+mgw(i)
    !       enddo
    !     !$OMP ENDDO NOWAIT
    !
    !     !$OMP DO
    !       do i = 1, no_halo
    !         d = sqrt((xow(i,1)-CoM(1))**2+(xow(i,2)-CoM(2))**2+(xow(i,3)-CoM(3))**2) * unit_l
    !         ind = int(d/rmax*nsamples)+1
    !         !$OMP ATOMIC
    !         prof_no(ind) = prof_no(ind)+1
    !         prof_mo(ind) = prof_mo(ind)+mow(i)
    !       enddo
    !     !$OMP ENDDO
    !   !$OMP END PARALLEL
    !
    !   dr = rmax/nsamples
    !
    !   do i=1, nsamples
    !     prof_ng(nsamples - i) = prof_ng(nsamples - i) + prof_ng(nsamples - i + 1)
    !     prof_mg(nsamples - i) = prof_mg(nsamples - i) + prof_mg(nsamples - i + 1)
    !     prof_no(nsamples - i) = prof_no(nsamples - i) + prof_no(nsamples - i + 1)
    !     prof_mo(nsamples - i) = prof_mo(nsamples - i) + prof_mo(nsamples - i + 1)
    !   enddo
    !
    !   do i = 1, nsamples
    !     prof_ng(i) = prof_ng(i) / (pi * ((i+0.5) * dr)**2 )
    !     prof_mg(i) = prof_mg(i) / (pi * ((i+0.5) * dr)**2 )
    !     prof_no(i) = prof_no(i) / (pi * ((i+0.5) * dr)**2 )
    !     prof_mo(i) = prof_mo(i) / (pi * ((i+0.5) * dr)**2 )
    !   enddo
    !
    !   ! if you want to do this, first you need to fix the units properly
    !   ! the computed density is wrong
    !
    !   ! call get_halo_string(halo_str, halo)
    !   ! fname = TRIM(srcdir//'/radial_galaxies_halo-'//TRIM(halo_str)//'.dat')
    !   ! open(unit=666, file=fname, form='formatted')
    !   ! write(666, '(5A15)') "#Distance [R200]", "N_no_orphans", "M_no_orphans", "N_with_orphans", "M_with_orphans"
    !   ! write(666, '(A,E15.6)') "#R200", R200
    !   ! do i = 0, nsamples-1
    !   !   write(666,'(5E15.6)') i*rmax/nsamples/R200, prof_ng(i), prof_mg(i), prof_no(i), prof_mo(i)
    !   ! enddo
    !   ! close(666)
    !
    ! end subroutine get_radial_profile

end program
