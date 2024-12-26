program get_correlation

  !==============================================
  ! cmd line args:
  ! get_correlation.o outputnr nc
  ! outputnr: 67 for output_00067
  ! nc is optional
  !==============================================

  use,intrinsic :: iso_c_binding
  use omp_lib
  use constants_and_parameters
  use io_module
  use SMHM_relations
  implicit none
  include 'fftw3.f03'

  !-------------------------------------------
  ! Global variables and parameters
  !-------------------------------------------

  ! Manual parameters
  integer             :: nc                  = 1024    ! number of cells per dimension
  integer, parameter  :: nsamples            = 200     ! number of samples for output; Number of bins for histograms
  character(len=3)    :: interpolation       = 'cic'   ! which interpolation to use: cic for cloud in cell, ngp for nearest grid point
  ! character(len=3)    :: interpolation      = 'ngp'    ! which interpolation to use: cic for cloud in cell, ngp for nearest grid point
  character(len=4)    :: SMHM_model           = 'behr'  ! if calculate_masses: which model to use: behroozi 2013
  ! character(len=4)    :: SMHM_model         = 'most' ! if calculate_masses: which model to use: moster 2013
  logical             :: check_normalisation = .false. ! print normalisation check via Parseval's theorem
  logical             :: calculate_masses    = .false. ! set =.true. if input is halo masses instead of galaxy masses




  !------------------
  ! Grid/Plot stuff
  !------------------
  real(dp)            :: dk, dx                        ! physical values of cell size in real/k space
  real(dp)            :: dc                            ! cell size in code units
  real(dp)            :: nccube                        ! (nc^3)
  integer             :: nchalf                        ! nchalf


  !-------------------------
  ! Arrays
  !-------------------------
  ! computed global arrays
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:,:)  :: Pk_field
  real(C_DOUBLE),            allocatable, dimension(:,:,:)  :: density_field
  real(C_DOUBLE),            allocatable, dimension(:)      :: Pk, correlation, projected_corr
  integer,                   allocatable, dimension(:)      :: Pk_counts, correlation_counts, projected_corr_counts
  real(dp),                  allocatable, dimension(:)      :: distances, distances_k, distances_p

  ! for wp
  ! real(C_DOUBLE),            allocatable, dimension(:,:)    :: density_field2d
  ! complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:)    :: Pk_field2d

  !-----------------------------
  ! Other variables
  !-----------------------------
  integer                          :: i, void, c, m
  real(dp)                         :: kmax, lmax, d
  character(len=4)                 :: add_name
  character(len=3), dimension(1:2) :: cases=(/'sub','all'/)
  real(dp), dimension(1:4)         :: thresholds=(/0.d0, 1d9, 1d10, 1d11/)

  ! FFTW stuff
  type(C_PTR)                      :: plan_r2c_3d, plan_c2r_3d
  ! type(C_PTR)                      :: plan_r2c_2d, plan_c2r_2d ! for wp


  ! checking normalisation
  real(dp) :: deltasq, delta_ksq, corrsq, pksq

  

  !---------------------------
  ! Set-up
  !---------------------------

  ! read necessary info
  call read_cmdlineargs(srcdir, nc) ! defined in this file, not io_module
  call read_info()        
  call read_clump_data(-1) ! -1: don't look for data of a specific halo
  call read_galaxy_data()


  ! Some definitions now that nc and unit_l are available
  dk         = 2*pi/unit_l_Mpc
  dx         = unit_l_Mpc/nc
  dc         = 1.d0/nc
  nchalf     = nc/2
  nccube     = nc**3.d0
  lmax       = 0.5d0*unit_l_Mpc
  kmax       = (nchalf+1) * dk
  ! lmax       = sqrt(3.d0)*unit_l_Mpc*1.000001
  ! kmax       = sqrt(3.d0)*(nchalf+1) * dk



  write(*,*) "============================"
  write(*,*) "get_correlation.f03 started."
  write(*,*) "============================"
  write(*,*)                "Working with parameters:"
  write(*,'(A20,A25)')      "  srcdir:        ", srcdir
  write(*,'(A20,I25)')      "  ncpu:          ", ncpu
  write(*,'(A20,I25)')      "  nc:            ", nc
  write(*,'(A20,I25)')      "  ngalaxies:     ", ngalaxies
  write(*,'(A20,A25)')      "  interpolation: ", interpolation
  write(*,'(A20,I25)')      "  levelmax:      ", levelmax
  write(*,'(A20,F25.3,A9)') "  unit_l:        ", unit_l_Mpc,  "[Mpc]"
  write(*,'(A20,F25.3,A9)') "  lmax:          ", lmax,    "[Mpc]"
  write(*,'(A20,F25.3,A9)') "  kmax:          ", kmax,    "[Mpc-1]"
  write(*,'(A20,F25.3)')    "  a_exp:         ", aexp
  write(*,'(A20,F25.3)')    "  z:             ", 1.d0/aexp - 1.d0



  if (calculate_masses) call calculate_sm(SMHM_model)



  !------------------------------
  ! Allocate arrays
  !------------------------------

  allocate(distances_k(0:nsamples))
  allocate(Pk(0:nsamples))
  allocate(Pk_counts(0:nsamples))
  allocate(distances(0:nsamples))
  allocate(correlation(0:nsamples))
  allocate(correlation_counts(0:nsamples))
  allocate(distances_p(0:nchalf))
  allocate(projected_corr(0:nchalf))
  allocate(projected_corr_counts(0:nchalf))


  allocate(density_field(1:nc, 1:nc, 1:nc))
  allocate(Pk_field(1:nchalf+1, 1:nc, 1:nc))

  ! allocate(density_field2d(1:nc, 1:nc))
  ! allocate(Pk_field2d(1:nchalf+1, 1:nc))




  !----------------------------------
  ! Initiate fftw stuff
  !----------------------------------
  void = fftw_init_threads()
  if (void==0) then
    write(*,*) "Error in fftw_init_threads, quitting"
    stop
  endif
  call fftw_plan_with_nthreads(omp_get_max_threads())
  plan_r2c_3d = fftw_plan_dft_r2c_3d(nc, nc, nc, density_field, Pk_field, FFTW_ESTIMATE)
  plan_c2r_3d = fftw_plan_dft_c2r_3d(nc, nc, nc, Pk_field, density_field, FFTW_ESTIMATE)
  ! plan_r2c_2d = fftw_plan_dft_r2c_2d(nc, nc, density_field2d, Pk_field2d, FFTW_ESTIMATE)
  ! plan_c2r_2d = fftw_plan_dft_c2r_2d(nc, nc, Pk_field2d, density_field2d, FFTW_ESTIMATE)


  !----------------------------------
  ! Get bin distances, reset values.
  !----------------------------------
  do i = 0, nsamples
    d = real(i,dp)/real(nsamples,dp)
    distances_k(i) = d*kmax
    distances(i) = d*lmax
  enddo
  ! do i = 0, nchalf
  !   distances_p(i) = real(i,dp)/real(nchalf)*lmax
  ! enddo

  !-------------------------------------------
  ! Do everything twice:
  ! If c = 1, don't include orphan galaxies
  ! If c = 2, include orphans
  !-------------------------------------------

  do m = 1, 4
    do c = 1, 2

      write(*,*)
      add_name = cases(c)

      if (check_normalisation) then
        deltasq = 0; delta_ksq = 0 ; pksq = 0; corrsq = 0
      endif


      ! reset arrays for this loop iteration
      call reset_arrays()
      call get_density_field(c, interpolation, thresholds(m))
      call get_Pk()
      call get_xi()
      ! call get_wp() ! deprecated, done in python now

      call write_results(thresholds(m))

    enddo
  enddo

  deallocate(Pk_field, density_field, distances, distances_k, Pk, correlation, Pk_counts, correlation_counts)
  ! deallocate(Pk_field2d, density_field2d)
  call fftw_destroy_plan(plan_r2c_3d)
  call fftw_destroy_plan(plan_c2r_3d)
  ! call fftw_destroy_plan(plan_r2c_2d)
  ! call fftw_destroy_plan(plan_c2r_2d)
  write(*,*) "get_correlation.f03 finished."

contains





  !================================
  subroutine get_Pk()
  !================================
    !---------------------------------------
    ! Computes and histogramms P(k)
    !---------------------------------------
    implicit none
    integer :: i,j,k,ix,iy,iz,ik,ig,f
    real(dp):: d

    !-----------------------------------------------------------------------
    write(*,*) "Computing and histogramming P(k) FFT for case ", add_name
    !-----------------------------------------------------------------------

    !-------------------------------
    ! Get P(k)
    !-------------------------------
    call fftw_execute_dft_r2c(plan_r2c_3d, density_field, Pk_field)

    !$OMP PARALLEL PRIVATE(i,j,k,ix,iy,iz,ik,ig,d,f)
      !$OMP DO COLLAPSE(3)
        do i=1, nchalf+1
          do j=1, nc
            do k=1, nc
              ! Note: With this normalisation, the "normalisation check" doesn't exactly
              ! work out, but you do less loops.
              ! Pk_field(i,j,k) = abs(Pk_field(i,j,k))**2.d0
              ! Pk_field(i,j,k) = abs(Pk_field(i,j,k)/nccube)**2.d0
              !todo: testing
              Pk_field(i,j,k) = abs(Pk_field(i,j,k)*dx**3)**2/unit_l_Mpc**3
            enddo
          enddo
        enddo
      !$OMP ENDDO

      if (check_normalisation) then
        !$OMP DO REDUCTION (+:delta_ksq, pksq)
          do i = 1, nchalf+1
            f=2
            if (i==1 .or. i==nchalf+1) f=1
            do j = 1, nc
              do k = 1, nc
                delta_ksq = delta_ksq + f*real(Pk_field(i,j,k))
                pksq = pksq + f*real(Pk_field(i,j,k))**2
              enddo
            enddo
          enddo
        !$OMP END DO
      endif


      !-------------------------------
      ! Histogram P(k)
      !-------------------------------
      !$OMP DO
        do i = 1, nchalf+1
          ix = i-1
          f = 2
          if (i==1 .or. i==nchalf+1) f=1
          do j = 1, nc
            iy = j-1
            if (iy>nchalf) iy = iy-nc
            do k = 1, nc
              iz = k-1
              if (iz>nchalf) iz = iz-nc

              d = sqrt(ix**2.d0+iy**2.d0+iz**2.d0)*dk
              if (d<=kmax .and. d>0.d0) then
                ! ig = int( log(d/dk)*nsamples/log(kmax/dk) )
                ig = int(d/kmax*nsamples)
                !$OMP ATOMIC
                Pk(ig) = Pk(ig)+f*real(Pk_field(i,j,k))
                !$OMP ATOMIC
                Pk_counts(ig) = Pk_counts(ig) + f
              endif
            enddo
          enddo
        enddo
      !$OMP END DO
    !$OMP END PARALLEL
  end subroutine get_Pk





  !================================
  subroutine get_xi()
  !================================
    !----------------------------------------------
    ! Computes and histogrammises 3d correlation
    !----------------------------------------------
    
    implicit none
    integer :: i,j,k,ix,iy,iz,ik,ig
    real(dp):: d

    !-----------------------------------------------------------------------
    write(*,*) "Computing and histogramming xi(r) FFT for case ", add_name
    !-----------------------------------------------------------------------

    !-------------------------------
    ! Get correlation
    !-------------------------------
    call fftw_execute_dft_c2r(plan_c2r_3d, Pk_field, density_field)

    !$OMP PARALLEL PRIVATE(i,j,k,ix,iy,iz,ik,ig,d)

      if (check_normalisation) then
        !$OMP DO COLLAPSE(3) REDUCTION (+:corrsq)
          do i = 1, nc
            do j = 1, nc
              do k = 1, nc
                corrsq=corrsq + density_field(i,j,k)**2
              enddo
            enddo
          enddo
        !$OMP END DO
      endif

      ! todo: temp
      !$OMP DO COLLAPSE(3)
        do i = 1, nc
          do j = 1, nc
            do k = 1, nc
              density_field(i,j,k) = density_field(i,j,k)*dk**3/twopicube
            enddo
          enddo
        enddo
      !$OMP END DO


      !-------------------------------
      ! Histogram correlation
      !-------------------------------

      !$OMP DO
        do i = 1, nc
          ix = i-1
          if (ix>nchalf) ix = ix-nc
          do j = 1, nc
            iy = j-1
            if (iy>nchalf) iy = iy-nc
            do k = 1,nc
              iz = k-1
              if (iz>nchalf) iz = iz-nc
              d = sqrt(ix**2.d0+iy**2.d0+iz**2.d0)*dx
              if (d<=lmax .and. d>0) then
                ig = int(d*nsamples/lmax)
                !$OMP ATOMIC
                correlation(ig) = correlation(ig)+density_field(i,j,k)
                !$OMP ATOMIC
                correlation_counts(ig) = correlation_counts(ig) + 1
              endif
            enddo
          enddo
        enddo
      !$OMP END DO



    !   !$OMP DO
        ! do i = 1, nc
        ! ! do i = 1, nc
        !   ix = i-1
        !   if (ix>nchalf) ix = ix-nc
        !   do j = 1, nc
        !     iy = j-1
        !     if (iy>nchalf) iy = iy-nc
        !     d = sqrt(ix**2.d0+iy**2.d0)*dx ! get perpendicular distance
        !     do k = 1, nc
        !       ! iz = j-1
        !       ! if (iz>nchalf) iz = iz-nc
        !       ! d = abs(ix)*dx
        !       if (d<=lmax .and. d>0) then
        !         ig = int(d*nchalf/lmax)
        !         !$OMP ATOMIC
        !         projected_corr(ig) = projected_corr(ig)+density_field(i,j,k)*dx ! this is an integral, not just a sum!
        !         !$OMP ATOMIC
        !         projected_corr_counts(ig) = projected_corr_counts(ig) + 1
        !       endif
        !     enddo
        !   enddo
        ! enddo
    !   !$OMP END DO

      ! !$OMP DO
      ! do ig = 0, nchalf
      !   projected_corr(ig) = projected_corr(ig)/nc
      ! enddo
      ! !$OMP END DO

    !$OMP END PARALLEL
  end subroutine get_xi







  !======================================================
  subroutine get_density_field(which, interp, mthresh)
  !======================================================
    !------------------------------------------------------------
    ! Computes the density field and overdensity field delta
    ! which = 1: main + satellites
    ! which = 2: main + satellites + orphans
    ! interp = 'ngp' : nearest grid point interpolation
    ! interp = 'cic' : cloud in cell interpolation
    ! mthresh: mass threshold (real)
    !------------------------------------------------------------

    implicit none
    integer, intent(in)  :: which
    real(dp), intent(in) :: mthresh
    character(len=3)     :: interp
    integer              :: i, j, k, gal, kept1, kept2
    integer              :: iup, idown, jup, jdown, kup, kdown
    real(dp)             :: rho, xup, yup, zup, hdc, cv, mean_density, mean_density_2d
    real(dp)             :: totdens, totmass
    logical              :: found_it

    hdc = dc/2
    cv = dc**3
    mean_density = 0
    mean_density_2d = 0
    kept1 = 0
    kept2 = 0
    totmass = 0
    totdens = 0

    !-----------------------------------------------------------------------------------
    write(*,*) "Computing density field for case ", add_name, "threshold =", mthresh
    !-----------------------------------------------------------------------------------


    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP PRIVATE(i,j,k,gal,found_it,iup,idown,jup,jdown,kup,kdown,rho,xup,yup,zup)
      if (interp=='ngp') then
        !$OMP DO REDUCTION(+:kept1)
          do gal=1, ngalaxies
            if (mg(gal) < mthresh) cycle
            i = int(xg(gal, 1)/boxlen*nc)+1
            j = int(xg(gal, 2)/boxlen*nc)+1
            k = int(xg(gal, 3)/boxlen*nc)+1
            !$OMP ATOMIC
            density_field(i,j,k) = density_field(i,j,k)+mg(gal)
            kept1 = kept1 + 1
          enddo
        !$OMP END DO
        if (which == 2) then
          !$OMP DO REDUCTION(+:kept1)
            do gal=1, norphans
              if (mo(gal) < mthresh) cycle
              i = int(xo(gal, 1)/boxlen*nc)+1
              j = int(xo(gal, 2)/boxlen*nc)+1
              k = int(xo(gal, 3)/boxlen*nc)+1
              !$OMP ATOMIC
              density_field(i,j,k) = density_field(i,j,k)+mo(gal)
              kept2 = kept2 + 1
            enddo
          !$OMP END DO
        endif

      else if (interp=='cic') then
        !$OMP DO REDUCTION(+:kept1)
          do gal=1, ngalaxies

            if (mg(gal) < mthresh) cycle

            iup   = int((xg(gal,1)+hdc)/boxlen*nc)+1
            idown = iup-1
            jup   = int((xg(gal,2)+hdc)/boxlen*nc)+1
            jdown = jup-1
            kup   = int((xg(gal,3)+hdc)/boxlen*nc)+1
            kdown = kup-1

            rho = mg(gal)/cv
            xup = xg(gal,1) + hdc - (iup-1)*dc
            yup = xg(gal,2) + hdc - (jup-1)*dc
            zup = xg(gal,3) + hdc - (kup-1)*dc

            if (iup>nc)  iup   = iup-nc
            if (idown<1) idown = nc+idown
            if (jup>nc)  jup   = jup-nc
            if (jdown<1) jdown = nc+jdown
            if (kup>nc)  kup   = kup-nc
            if (kdown<1) kdown = nc+kdown

            !$OMP CRITICAL
              density_field(iup,   jup,   kup)   = density_field(iup,   jup,   kup)  + xup      * yup      * zup      * rho
              density_field(idown, jup,   kup)   = density_field(idown, jup,   kup)  + (dc-xup) * yup      * zup      * rho
              density_field(iup,   jdown, kup)   = density_field(iup,   jdown, kup)  + xup      * (dc-yup) * zup      * rho
              density_field(idown, jdown, kup)   = density_field(idown, jdown, kup)  + (dc-xup) * (dc-yup) * zup      * rho
              density_field(iup,   jup,   kdown) = density_field(iup,   jup,   kdown)+ xup      * yup      * (dc-zup) * rho
              density_field(idown, jup,   kdown) = density_field(idown, jup,   kdown)+ (dc-xup) * yup      * (dc-zup) * rho
              density_field(iup,   jdown, kdown) = density_field(iup,   jdown, kdown)+ xup      * (dc-yup) * (dc-zup) * rho    
              density_field(idown, jdown, kdown) = density_field(idown, jdown, kdown)+ (dc-xup) * (dc-yup) * (dc-zup) * rho
            !$OMP END CRITICAL
            kept1 = kept1 + 1
          enddo
        !$OMP END DO

        if (which==2) then
          !$OMP DO REDUCTION(+:kept1)
            do gal=1, norphans

              if (mo(gal) < mthresh) cycle

              iup   = int((xo(gal,1)+hdc)/boxlen*nc)+1
              idown = iup-1
              jup   = int((xo(gal,2)+hdc)/boxlen*nc)+1
              jdown = jup-1
              kup   = int((xo(gal,3)+hdc)/boxlen*nc)+1
              kdown = kup-1

              rho = mo(gal)/cv
              xup = xo(gal,1) + hdc - (iup-1)*dc
              yup = xo(gal,2) + hdc - (jup-1)*dc
              zup = xo(gal,3) + hdc - (kup-1)*dc

              if (iup>nc)  iup   = iup-nc
              if (idown<1) idown = nc+idown
              if (jup>nc)  jup   = jup-nc
              if (jdown<1) jdown = nc+jdown
              if (kup>nc)  kup   = kup-nc
              if (kdown<1) kdown = nc+kdown

              !$OMP CRITICAL
                density_field(iup,   jup,   kup)   = density_field(iup,   jup,   kup)  + xup      * yup      * zup      * rho
                density_field(idown, jup,   kup)   = density_field(idown, jup,   kup)  + (dc-xup) * yup      * zup      * rho
                density_field(iup,   jdown, kup)   = density_field(iup,   jdown, kup)  + xup      * (dc-yup) * zup      * rho
                density_field(idown, jdown, kup)   = density_field(idown, jdown, kup)  + (dc-xup) * (dc-yup) * zup      * rho
                density_field(iup,   jup,   kdown) = density_field(iup,   jup,   kdown)+ xup      * yup      * (dc-zup) * rho
                density_field(idown, jup,   kdown) = density_field(idown, jup,   kdown)+ (dc-xup) * yup      * (dc-zup) * rho
                density_field(iup,   jdown, kdown) = density_field(iup,   jdown, kdown)+ xup      * (dc-yup) * (dc-zup) * rho    
                density_field(idown, jdown, kdown) = density_field(idown, jdown, kdown)+ (dc-xup) * (dc-yup) * (dc-zup) * rho
              !$OMP END CRITICAL
              kept2 = kept2 + 1
            enddo
          !$OMP END DO
        endif

      else
        write(*,*) "Didn't recognize interpolation method ", interp
        stop
      endif
        

      !$OMP DO COLLAPSE(3) REDUCTION(+: mean_density)
        do i=1, nc
          do j=1, nc
            do k=1, nc
              density_field(i,j,k) = density_field(i,j,k)/volume_Mpc
              mean_density = mean_density + density_field(i,j,k)
            enddo
          enddo
        enddo
      !$OMP END DO

      ! k=20
      ! !$OMP DO COLLAPSE(2) REDUCTION(+: mean_density_2d)
      ! do i=1, nc
      !   do j=1, nc
      !     ! density_field2d(i,j) = sum(density_field(i,j,:))
      !     density_field2d(i,j) = density_field(i,j,k)
      !     ! mean_density_2d = mean_density_2d + density_field2d(i,j)
      !   enddo
      ! enddo
      ! !$OMP ENDDO

      !$OMP SINGLE
        ! mean_density_2d = mean_density_2d/nc**2
        mean_density = mean_density/nccube
      !$OMP END SINGLE

      !$OMP FLUSH(mean_density, mean_density_2d)

      !$OMP DO REDUCTION(+:totdens, totmass)
        do i=1, nc
          do j = 1, nc
            do k = 1, nc
              totdens = totdens + density_field(i,j,k)
              totmass = totmass + density_field(i,j,k) * cv
            enddo
          ! density_field2d(i,j) = density_field2d(i,j)/mean_density_2d -1.d0
          enddo
        enddo
      !$OMP END DO


      ! Compute delta(x)
      !$OMP DO COLLAPSE(3)
        do i=1, nc
          do j = 1, nc
            do k = 1, nc
              density_field(i,j,k) = (density_field(i,j,k)/mean_density - 1.d0)
            enddo
          ! density_field2d(i,j) = density_field2d(i,j)/mean_density_2d -1.d0
          enddo
        enddo
      !$OMP END DO

      if (check_normalisation) then
        !$OMP DO COLLAPSE(3) REDUCTION (+:deltasq)
          do i = 1, nc
            do j = 1, nc
              do k = 1, nc
                deltasq = deltasq + density_field(i,j,k)**2
              enddo
            enddo
          enddo
        !$OMP END DO
      endif
    !$OMP END PARALLEL

    write(*, '(A, I8, A, I8)') "Finished density field. Kept ", &
                                kept1 + kept2, " galaxies out of ", ngalaxies
    write(*, '(A, 1PE12.3E2, A, 1PE12.3E2)') "total density ", &
                                totdens, " total mass ", totmass

  end subroutine get_density_field




  !===========================================
  subroutine read_cmdlineargs(srcdir, ncells)
  !===========================================
    !-----------------------
    ! reads cmd line args
    !-----------------------
    character(len=12), intent(out) :: srcdir
    integer, intent(out) :: ncells

    character(len=12) :: arg
    character(len=5)  :: dirnr_str
    integer           :: ncells_temp, dirnr
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
        read(arg, *) ncells_temp
        if (ncells_temp>=0) ncells = ncells_temp
      endif
    end do
  end subroutine read_cmdlineargs



  
  !===========================
  subroutine reset_arrays()
  !===========================
    !----------------------------
    ! Resets array values to 0.
    !----------------------------

    implicit none
    integer :: i,j,k

    !$OMP PARALLEL DEFAULT(shared) PRIVATE(i,j,k)
      !$OMP DO
      do i=0, nsamples
        Pk_counts(i)         =0
        Pk(i)                =0
        correlation_counts(i)=0
        correlation(i)       =0
      enddo
      !$OMP ENDDO NOWAIT
      !$OMP DO
      do i=0, nchalf
        projected_corr(i)       =0
        projected_corr_counts(i)=0
      enddo
      !$OMP ENDDO NOWAIT

      !$OMP DO
        do i=1, nc
          do j=1, nc
            ! density_field2d(i,j) = 0.d0
            do k=1, nc
              density_field(i,j,k) = 0.d0
            enddo
            do k=1, nc/2+1
              Pk_field(k,i,j) = 0.d0
            enddo
          enddo
          ! do j=1, nchalf
          !   Pk_field2d(j,i) = 0.d0
          ! enddo
        enddo
      !$OMP END DO

    !$OMP END PARALLEL
  end subroutine reset_arrays



  !==================================
  subroutine write_results(threshold)
  !==================================
    implicit none
    real(dp), intent(in) :: threshold
    character(len=80) :: outfname
    character(len=8) :: thresh_str
    integer :: i

    write(thresh_str, '(1PE8.2E2)') threshold

    !$OMP PARALLEL DEFAULT(shared) PRIVATE(i,outfname)
      !$OMP SINGLE
        outfname=TRIM(TRIM(srcdir)//"/Pk_"//TRIM(add_name)//"-"//TRIM(thresh_str)//".txt")
        open(unit=666, file=outfname, form='formatted')
        write(666,*) "Results generated by get_correlation.f03"
        write(666,'(3A22,x)') "k", "P(k)", "Pk(k) counts"
        do i = 0, nsamples
          write(666, '(2(E22.10,x), I22)') distances_k(i), Pk(i), Pk_counts(i)
        enddo
        write(*,*) "Finished P(k), written results to "//TRIM(outfname)
        close(666)
      !$OMP END SINGLE

      !$OMP SINGLE
        outfname=TRIM(TRIM(srcdir)//"/correlation_"//TRIM(add_name)//".txt")
        outfname=TRIM(TRIM(srcdir)//"/correlation_"//TRIM(add_name)//"-"//TRIM(thresh_str)//".txt")
        open(unit=667, file=outfname, form='formatted')
        write(667,*) "Results generated by get_correlation.f03"
        write(667,'(3A22,x)') "r", "correlation(r)", "correlation(r) counts"
        do i = 0, nsamples
          write(667, '(2(E22.10,x), I22)') distances(i), correlation(i), correlation_counts(i)
        enddo
        write(*,*) "Finished xi(r), written results to "//TRIM(outfname)
        close(667)
      !$OMP END SINGLE

      ! !$OMP SINGLE
      !   outfname=TRIM(TRIM(srcdir)//"/projected_correlation_"//TRIM(add_name)//".txt")
      !   open(unit=667, file=outfname, form='formatted')
      !   write(667,*) "Results generated by get_correlation.f03"
      !   write(667,'(3A22,x)') "r_p", "w_p(r_p)", "w_p(r_p) counts"
      !   do i = 1, nchalf
      !     write(667, '(2(E22.10,x), I22)') 0.5*(distances_p(i-1)+distances_p(i)), projected_corr(i), projected_corr_counts(i)
      !   enddo
      !   write(*,*) "Finished wp(rp), written results to "//TRIM(outfname)
      !   close(667)
      ! !$OMP END SINGLE

      !$OMP SINGLE
        if (check_normalisation) then
          write(*,'(A25,3E18.10)') "Normalisation check: R->F", deltasq, delta_ksq, deltasq/delta_ksq
          write(*,'(A25,3E18.10)') "Normalisation check: F->R", pksq, corrsq, pksq/corrsq
        endif
      !$OMP END SINGLE
    !$OMP END PARALLEL
  end subroutine write_results



    ! !================================
    ! subroutine get_wp()
    ! !================================
    !
    !   DEPRECATED, DONE IN PYTHON NOW
    !
    !   !----------------------------------------
    !   ! Compute projected correlation
    !   !----------------------------------------
    !
    !   implicit none
    !
    !   integer :: i,j,ix,iy,ig
    !   real(dp):: d
    !
    !   call fftw_execute_dft_r2c(plan_r2c_2d, density_field2d, Pk_field2d)
    !
    !   !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j)
    !     do i=1, nchalf+1
    !       do j=1, nc
    !         Pk_field2d(i,j) = abs(Pk_field2d(i,j)*dx**2)**2/unit_l_Mpc**2
    !       enddo
    !     enddo
    !   !$OMP END PARALLEL DO
    !
    !   call fftw_execute_dft_c2r(plan_c2r_2d, Pk_field2d, density_field2d)
    !
    !   !$OMP PARALLEL PRIVATE(i,j, ix, iy, d, ig)
    !
    !     !$OMP DO COLLAPSE(2)
    !       do i = 1, nc
    !         do j = 1, nc
    !           density_field2d(i,j) = density_field2d(i,j)*dk**2/(2*pi)**2
    !         enddo
    !       enddo
    !     !$OMP END DO
    !
    !     !----------------------------------
    !     ! Histogram projected correlation
    !     !----------------------------------
    !
    !     !$OMP DO
    !       do i = 1, nc
    !         ix = i-1
    !         if (ix>nchalf) ix = ix-nc
    !         do j = 1, nc
    !           iy = j-1
    !           if (iy>nchalf) iy = iy-nc
    !           ! d = abs(ix)*dx
    !           d = sqrt(ix**2.d0+iy**2.d0)*dx ! get perpendicular distance
    !           if (d<=lmax .and. d>0) then
    !             ig = int(d*nchalf/lmax)
    !             !$OMP ATOMIC
    !             projected_corr(ig) = projected_corr(ig)+density_field2d(i,j)*dx! this is an integral, not just a sum!
    !             !$OMP ATOMIC
    !             projected_corr_counts(ig) = projected_corr_counts(ig) + 1
    !           endif
    !         enddo
    !       enddo
    !     !$OMP ENDDO
    !
    !   !$OMP END PARALLEL
    ! end subroutine get_wp

end program
