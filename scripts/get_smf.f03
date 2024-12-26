program get_smf

  use omp_lib
  use constants_and_parameters
  use io_module
  use density_projection
  implicit none

  !----------------------
  ! Manual parameters
  !----------------------
  integer             :: output_start = 21          ! first snapshot to use
  integer             :: output_end   = 67          ! last snapshot to use
  integer, parameter  :: nsamples     = 100         ! number of samples for output; Number of bins for histograms
  real(dp), parameter :: minmass = 4.d0
  real(dp), parameter :: maxmass = 13.0d0
  logical             :: calculate_masses = .false. ! set =.true. if input is halo masses instead of galaxy masses
  character(len=4)    :: SMHM_model = 'behr'        ! if calculate_masses: which model to use: behroozi 2013
  ! character(len=4)    :: SMHM_model = 'SMHM_most'      ! if calculate_masses: which model to use: moster 2013
  logical             :: apply_filter = .true.     ! whether to exclude orphans that are closer than 
                                                    ! filter_distance to a non-orphan galaxy
  real(dp), parameter :: filter_distance = 20.      ! kpc; exclude orphans that are closer than this distance from a real galaxy


  !---------------------
  ! Arrays
  !---------------------

  ! computed global arrays
  real(dp), allocatable, dimension(:) :: logmass         ! log10(stellar mass) bins
  logical :: pic = .true.               ! whether to create density projection image
  


  !-----------------------------
  ! Other variables
  !-----------------------------
  integer                          :: i, outp
  real(dp)                         :: dm
  character(len=5)                 :: dirnr
  character(len=80)                :: fname




  !----------------------------
  ! Preparations
  !----------------------------
  nc = 1000        ! number of cells for density field
  interpolation = 'cic'

  ! read necessary info
  call title(output_end, dirnr)
  srcdir = TRIM('output_'//dirnr)
  call read_info()        

  ! Allocate and fill constant global arrays
  allocate(logmass(0:nsamples));  logmass =0;
  dm = (maxmass-minmass)/nsamples
  do i = 0, nsamples 
    logmass(i) = minmass+i*dm
  enddo
  





  write(*,*) "============================"
  write(*,*) "get_smf.f03 started."
  write(*,*) "============================"
  write(*,*)                "Working with parameters:"
  write(*,'(A20,I25)')      "  ncpu:          ", ncpu
  write(*,'(A20,I25)')      "  nc:            ", nc
  write(*,'(A20,1PE25.8)')  "  log10(M_S)_min:", minmass
  write(*,'(A20,1PE25.8)')  "  log10(M_S)_max:", maxmass
  write(*,'(A20,A25)')      "  interpolation: ", interpolation
  if (calculate_masses) then
  write(*,'(A20,A25)')      "  SM model:      ", SMHM_model
  endif
  if (pic) then
  write(*,'(A20,A25)')      "  make picture:  ", "yes"
  endif




  !--------------------------
  ! MAIN LOOP
  !--------------------------

  do outp = output_start, output_end

    ! Reset stuff
    ngalaxies = 0

    call title(outp, dirnr)
    srcdir = TRIM('output_'//dirnr)

    call read_info()  ! from io_module
    call read_clump_data(-1) ! -1: don't look for particular data of a specific halo
    call read_galaxy_data()

    write(*,'(A60)')          "-------------------------------------------------------------"
    write(*,'(A20,A25)')      "  srcdir:        ", srcdir
    write(*,'(A20,I25)')      "  ngalaxies:     ", ngalaxies
    write(*,'(A20,F25.3,A)')  "  unit_l:        ", unit_l_cMpc, " comoving Mpc"
    write(*,'(A20,F25.3)')    "  a_exp:         ", aexp
    write(*,'(A20,F25.3)')    "  z:             ", 1.d0/aexp - 1.d0

    ! if (skip_this) then
    !   write(*,'(A45)') "   No data found. Skipping this snapshot."
    !   cycle
    ! endif


    !-----------------------------------------------
    write(*,*) "Calculating Stellar Mass Function"
    !-----------------------------------------------

    if (calculate_masses) call calculate_sm(SMHM_model)
    call smf()


    if (pic) then
      call allocate_density_fields()

      ! fist galaxies only
      call get_density_projection(xg, mg, ngalaxies)
      fname = TRIM(srcdir//'/density_image_galaxies-sub.dat')
      call write_density_projections(fname)

      ! now add orphans
      call get_density_projection(xo, mo, norphans)
      fname = TRIM(srcdir//'/density_image_galaxies-all.dat')
      call write_density_projections(fname)

      call deallocate_density_fields()
    endif

    call deallocate_io_module()

  enddo ! loop over all directories

 

  write(*,*) "get_smf.f03 finished."
  deallocate(logmass)


contains



  ! !===========================================
  ! subroutine get_density_field(which, interp)
  ! !===========================================
  !   !------------------------------------------------------------
  !   ! Computes the density field and overdensity field delta
  !   ! which = 1: main + satellites
  !   ! which = 2: main + satellites + orphans
  !   ! interp = 'ngp' : nearest grid point interpolation
  !   ! interp = 'cic' : cloud in cell interpolation
  !   !------------------------------------------------------------
  !
  !   implicit none
  !   integer, intent(in)  :: which
  !   character(len=3)     :: interp
  !   integer              :: i, j, gal
  !   integer              :: iup, idown, jup, jdown, kup, kdown
  !   real(dp)             :: rho, xup, yup, zup, hdc, cv
  !   logical              :: found_it
  !
  !
  !   hdc = dc/2
  !   cv = dc**3
  !
  !   !$OMP PARALLEL DEFAULT(shared) PRIVATE(i,j,gal,found_it,iup,idown,jup,jdown,kup,kdown,rho,xup,yup,zup)
  !     !$OMP DO COLLAPSE(2)
  !       ! reset density
  !       do i=1, nc
  !         do j=1, nc
  !           density_field(i,j) = 0.d0
  !         enddo
  !       enddo
  !     !$OMP END DO
  !
  !     if (interp=='ngp') then
  !       !$OMP DO
  !         do gal=1, ngalaxies
  !           if (which==2 .or. gal_id(gal)>0) then
  !             i = int(x(gal)/boxlen*nc)+1
  !             j = int(y(gal)/boxlen*nc)+1
  !             !$OMP ATOMIC
  !             density_field(i,j) = density_field(i,j)+stellar_mass(gal)
  !           endif
  !         enddo
  !       !$OMP END DO
  !
  !     else if (interp=='cic') then
  !       !$OMP DO
  !         do gal=1, ngalaxies
  !           if (which==2 .or. gal_id(gal)>0) then
  !             iup   = int((x(gal)+hdc)/boxlen*nc)+1
  !             idown = iup-1
  !             jup   = int((y(gal)+hdc)/boxlen*nc)+1
  !             jdown = jup-1
  !
  !             rho = stellar_mass(gal)/cv
  !             xup = x(gal) + hdc - (iup-1)*dc
  !             yup = y(gal) + hdc - (jup-1)*dc
  !
  !             if (iup>nc)  iup   = iup-nc
  !             if (idown<1) idown = nc+idown
  !             if (jup>nc)  jup   = jup-nc
  !             if (jdown<1) jdown = nc+jdown
  !
  !             !$OMP CRITICAL
  !               density_field(iup,   jup  )   = density_field(iup,   jup  )  + xup      * yup      * rho
  !               density_field(idown, jup  )   = density_field(idown, jup  )  + (dc-xup) * yup      * rho
  !               density_field(iup,   jdown)   = density_field(iup,   jdown)  + xup      * (dc-yup) * rho
  !               density_field(idown, jdown)   = density_field(idown, jdown)  + (dc-xup) * (dc-yup) * rho
  !             !$OMP END CRITICAL
  !
  !           endif
  !         enddo
  !       !$OMP END DO
  !
  !     else
  !       write(*,*) "Didn't recognize interpolation method ", interp
  !       stop
  !     endif
  !
  !
  !     !$OMP DO COLLAPSE(2)
  !       do i=1, nc
  !         do j=1, nc
  !           density_field(i,j) = density_field(i,j)/volume_cMpc
  !         enddo
  !       enddo
  !     !$OMP END DO
  !   !$OMP END PARALLEL
  ! end subroutine get_density_field






  !========================================
  subroutine smf()
  !========================================
    !-----------------------------------------------
    ! Computes and writes stellar mass function
    !-----------------------------------------------

    real(dp), allocatable, dimension(:) :: smf_main
    real(dp), allocatable, dimension(:) :: smf_sub
    real(dp), allocatable, dimension(:) :: smf_orph
    real(dp), allocatable, dimension(:) :: smf_orph_filtered
    real(dp)                            :: lm, f, r, dx, dy, dz
    integer                             :: i,guess,g,nskipped_orph
    character(len=80)                   :: fname
    logical                             :: skip_orphan
    real(dp)                            :: r_in_kpc, cutoff

    ! get factor for physical distance between two particles
    ! factor aexp: unit_l is in comoving distance
    ! factor 1000: filter distance is given in kpc
    r_in_kpc = unit_l_Mpc * 1000.
    cutoff = filter_distance / r_in_kpc

    allocate(smf_main(0:nsamples)); smf_main=0;
    allocate(smf_sub(0:nsamples)); smf_sub=0;
    allocate(smf_orph(0:nsamples)); smf_orph=0;
    allocate(smf_orph_filtered(0:nsamples)); smf_orph_filtered=0;


    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP PRIVATE(i,guess,lm,f,g,r,skip_orphan,dx,dy,dz)
      ! "real" galaxies first
      !$OMP DO
        do i = 1, ngalaxies
          ! determine index
          lm = log10(mg(i))
          guess = int((lm-minmass)/dm) + 1
          if (guess<1) then
            write(*,'(A9,x,E12.5, A42, E12.5)') "Got mass", lm, " which is smaller than your set minmass of", minmass
            write(*,*) "Consider lowering the threshold."
            cycle
          endif

          if (is_central(i)) then
            ! found main
            !$OMP ATOMIC
            smf_main(guess) = smf_main(guess) + 1
          else
            ! found satellite
            !$OMP ATOMIC
            smf_sub(guess) = smf_sub(guess) + 1
          endif
        enddo
      !$OMP END DO NOWAIT

      !$OMP DO
        do i = 1, norphans
          ! determine index
          lm = log10(mo(i))
          guess = int((lm-minmass)/dm) + 1
          !$OMP ATOMIC
          smf_orph(guess) = smf_orph(guess) + 1

          ! filtered orphans
          if (apply_filter) then
            skip_orphan = .false.
            do g = 1, ngalaxies
              dx = xo(i, 1) - xg(g, 1)
              dy = xo(i, 2) - xg(g, 2)
              dz = xo(i, 3) - xg(g, 3)
              r = sqrt(dx**2 + dy**2 + dz**2)
              if (r <= cutoff) then
                skip_orphan = .true.
                exit
              endif
            enddo
            if (.not.skip_orphan) then
              !$OMP ATOMIC
              smf_orph_filtered(guess) = smf_orph_filtered(guess) + 1
            else
              !$OMP ATOMIC
              nskipped_orph = nskipped_orph + 1
            endif
          endif
        enddo
      !$OMP ENDDO

      ! !$OMP DO
      !   do i=1, nsamples
      !     f = (logmass(i)-logmass(i-1))*volume_cMpc
      !     smf_all(i) = smf_all(i)/f
      !     smf_sub(i) = smf_sub(i)/f
      !     smf_main(i) = smf_main(i)/f
      !   enddo
      ! !$OMP ENDDO

    !$OMP END PARALLEL

    if (apply_filter) then
      write(*, '(A,I12,A,I12)') "Finished computing SMFs. Skipped", &
                                nskipped_orph, "/", norphans
    endif


    ! old version; changed so I don't overwrite things on accident
    ! old version = pre 13.04.2021
    ! fname=TRIM(TRIM(srcdir)//"/smf-new.txt")
    ! open(unit=666, file=fname, form='formatted')
    ! write(666, '(2A20)') "Mass bins", "SMF counts"

    fname=TRIM(TRIM(srcdir)//"/smf-new.txt")
    open(unit=666, file=fname, form='formatted')
    write(666, '(A2, 2A20)') "# ", "Mass_bins", "SMF_counts"
    do i=0, nsamples
      write(666, '(F20.6, E20.12)') logmass(i), smf_main(i)
    enddo
    close(666)

    write(*,*) "SMF for main haloes stuff finished, results written to ", fname

    fname=TRIM(TRIM(srcdir)//"/smf-including-satellites.txt")
    open(unit=666, file=fname, form='formatted')
    if (apply_filter) then
      write(666, '(A,x,F12.3,x,A)') "# Cutoff Radius for filter:", filter_distance, "kpc"
      write(666, '(A2,4A20,A27)') "# ", "Mass_bins", "SMF_counts_main", &
                                "SMF_counts_sub", "SMF_counts_orphan", "SMF_counts_orphan_filtered"
      do i=0, nsamples
        write(666, '(F20.6, 4E20.12)') logmass(i), smf_main(i), smf_sub(i), smf_orph(i), smf_orph_filtered(i)
      enddo
    else
      write(666, '(A,x,F12.3,x,A)') "# Cutoff Radius for filter: NONE; NO FILTER USED"
      write(666, '(A2,4A20,A27)') "# ", "Mass_bins", "SMF_counts_main", &
                                "SMF_counts_sub", "SMF_counts_orphan", "SMF_counts_orphan_filtered"
      do i=0, nsamples
        write(666, '(F20.6, 4E20.12)') logmass(i), smf_main(i), smf_sub(i), smf_orph(i), 0.
      enddo
    endif
    close(666)

    write(*,*) "SMF for all haloes stuff finished, results written to ", fname


    deallocate(smf_main)
  end subroutine smf


end program
