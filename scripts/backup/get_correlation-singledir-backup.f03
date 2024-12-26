program eval_galaxies

  use,intrinsic :: iso_c_binding
  use omp_lib
  implicit none
  include 'fftw3.f03'

  

  integer, parameter  :: dp = kind(1.0d0) 

  !-------------------------
  ! Manual parameters
  !-------------------------
  integer, parameter  :: nsamples = 100           ! number of samples for output; Number of bins for histograms
  character(len=3)    :: interpolation = 'cic'    ! which interpolation to use: cic for cloud in cell, ngp for nearest grid point
  ! character(len=3)    :: interpolation = 'ngp'    ! which interpolation to use: cic for cloud in cell, ngp for nearest grid point



  !-------------------------
  ! Constants
  !-------------------------
  real(dp), parameter :: pi = 3.14159265359
  real(dp), parameter :: twopicube = (2.d0*pi)**3
  real(dp), parameter :: Mpc = 3.086e+24
  real(dp), parameter :: boxlen = 1.d0



  !-----------------------
  ! Read-in data
  !-----------------------
  real(dp)            :: aexp                     ! expansion factor
  real(dp)            :: H0                       ! Hubble factor [km s-1 Mpc-1]
  real(dp)            :: unit_l                   ! length of box [Mpc]
  real(dp)            :: volume                   ! current physical volume
  integer             :: nhalos, ngalaxies, ncpu
  character(len=12)   :: srcdir                   ! output_XXXXX 


  !------------------------
  ! Grid/plot stuff
  !------------------------
  real(dp)            :: dc                       ! cell size in code units
  real(dp)            :: dk, dx                   ! physical values of cell size in real/k space
  real(dp)            :: nccube, sqrtnccube, nccubesq

  !---------------------
  ! Arrays
  !---------------------
  ! read-in data arrays
  real(dp), dimension(:), allocatable :: x, y, z          ! x,y,z positions of galaxies [code units]
  real(dp), dimension(:), allocatable :: stellar_mass     ! stellar mass of galaxies [M_Sol]
  integer,  dimension(:), allocatable :: halo_id, clmp_id ! ID's of halos in snapshot and associated clump IDs

  ! computed global arrays
  integer  :: nc
  complex(C_DOUBLE_COMPLEX), allocatable, dimension(:,:,:) :: Pk_field
  real(C_DOUBLE),            allocatable, dimension(:,:,:) :: density_field
  real(C_DOUBLE),            allocatable, dimension(:)     :: Pk, correlation
  integer,                   allocatable, dimension(:)     :: Pk_counts, correlation_counts
  real(dp),                  allocatable, dimension(:)     :: distances, distances_k

  ! FFTW stuff
  type(C_PTR)                                              :: plan_forward, plan_backward


  !-----------------------------
  ! Other variables
  !-----------------------------
  integer                          :: i,j,k, void, ix, iy, iz, ik, ig, c, f
  real(dp)                         :: mean_density, kmax, lmax, d
  character(len=80)                :: outfname
  character(len=4)                 :: add_name
  character(len=3), dimension(1:2) :: cases=(/'sub','all'/)

  ! checking normalisation
  real(dp) :: deltasq, delta_ksq, corrsq, pksq
  logical  :: check_normalisation = .true.





  write(*,*) "============================"
  write(*,*) "get_correlation.f03 started."
  write(*,*) "============================"


  !----------------------------------
  write(*,*) "Reading data."
  !----------------------------------

  nc = 128
  call read_cmdlineargs() ! read command line arguments
  call read_info()        ! read data from info files
  call read_data()        ! read actual data and nhalos, ngalaxies


  ! Some definitions
  lmax       = 0.5d0*unit_l
  dc         = 1.d0/nc
  dx         = unit_l/nc
  dk         = 2*pi/unit_l
  kmax       = (nc/2) * dk  ! choose max k to be at half "box size" so you have full spheric samples only
  volume     = unit_l**3
  nccube     = nc**3
  nccubesq   = nccube**2
  sqrtnccube = sqrt(nccube)


  write(*,*)                "Working with parameters:"
  write(*,'(A20,A25)')      "  srcdir:        ", srcdir
  write(*,'(A20,I25)')      "  ncpu:          ", ncpu
  write(*,'(A20,I25)')      "  nc:            ", nc
  write(*,'(A20,I25)')      "  nhaloes:       ", nhalos
  write(*,'(A20,I25)')      "  ngalaxies:     ", ngalaxies
  write(*,'(A20,E25.8,A9)') "  lmax:          ", lmax, "[Mpc]"
  write(*,'(A20,E25.8,A9)') "  kmax:          ", kmax, "[Mpc-1]"
  write(*,'(A20,A25)')      "  interpolation: ", interpolation




  !-----------------------------------------------
  ! write(*,*) "Calculating Stellar Mass Function"
  !-----------------------------------------------
  ! call smf()




  allocate(distances(0:nsamples));                distances=0
  allocate(distances_k(0:nsamples));              distances_k=0;
  allocate(Pk_counts(0:nsamples));                Pk_counts=0;
  allocate(correlation_counts(0:nsamples));       correlation_counts=0;
  allocate(Pk(1:nsamples));                       Pk=0;
  allocate(correlation(1:nsamples));              correlation=0;

  allocate(density_field(1:nc, 1:nc, 1:nc))
  allocate(Pk_field(1:nc/2+1, 1:nc, 1:nc))




  !----------------------------------
  ! Initiate fftw stuff
  !----------------------------------
  void = fftw_init_threads()
  if (void==0) then
    write(*,*) "Error in fftw_init_threads, quitting"
    stop
  endif
  call fftw_plan_with_nthreads(omp_get_max_threads())
  plan_forward  = fftw_plan_dft_r2c_3d(nc, nc, nc, density_field, Pk_field, FFTW_ESTIMATE)
  plan_backward = fftw_plan_dft_c2r_3d(nc, nc, nc, Pk_field, density_field, FFTW_ESTIMATE)





  !-------------------------------------------
  ! Do everything twice:
  ! If c = 1, include orphans
  ! If c = 2, don't include orphan galaxies
  !-------------------------------------------

  ! do c = 1,1
  do c = 1, 2
    
    write(*,*) 
    add_name = cases(c)

                                                                    deltasq = 0
                                                                    delta_ksq = 0
                                                                    pksq = 0
                                                                    corrsq = 0

    !---------------------------------------------------------
    write(*,*) "Computing density field for case ", add_name
    !---------------------------------------------------------

    call get_density_field(c, interpolation)
    ! call write_galaxy_images() ! call before changing density


    mean_density = 0
    !$OMP PARALLEL DEFAULT(shared) PRIVATE(i,j,k)
      !$OMP DO COLLAPSE(3) REDUCTION (+:mean_density)
        do i = 1, nc
          do j = 1, nc
            do k = 1, nc
              mean_density = mean_density + density_field(i,j,k)
            enddo
          enddo
        enddo
      !$OMP END DO

      !$OMP SINGLE
        mean_density = mean_density/nccube
      !$OMP END SINGLE

      !$OMP FLUSH(mean_density)

      !$OMP DO
      ! !$OMP DO COLLAPSE(3)
        do i=1, nc
          do j = 1, nc
            do k = 1, nc
              density_field(i,j,k) = (density_field(i,j,k)/mean_density - 1.d0)
            enddo
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




    !-----------------------------------------------------------------------
    write(*,*) "Computing and histogramming P(k) FFT for case ", add_name
    !-----------------------------------------------------------------------

    call fftw_execute_dft_r2c(plan_forward, density_field, Pk_field)

    !$OMP PARALLEL PRIVATE(i,j,k,ix,iy,iz,ik,ig,d,f)
      !$OMP DO COLLAPSE(3)
        do i=1, nc/2+1
          do j=1, nc
            do k=1, nc
              ! Note: With this normalisation, the "normalisation check" doesn't exactly 
              ! work out, but you do less loops.
              Pk_field(i,j,k) = (real(Pk_field(i,j,k))**2 + aimag(Pk_field(i,j,k))**2)/nccube
            enddo
          enddo
        enddo
      !$OMP ENDDO

      if (check_normalisation) then
        !$OMP DO REDUCTION (+:delta_ksq, pksq)
          do i = 1, nc/2+1
            f=2
            if (i==1 .or. i==nc/2+1) f=1
            do j = 1, nc
              do k = 1, nc
                delta_ksq = delta_ksq + f*real(Pk_field(i,j,k))
                pksq = pksq + f*real(Pk_field(i,j,k))**2
              enddo
            enddo
          enddo
        !$OMP END DO
      endif


      !$OMP DO
        ! Get bin distances, reset values. Use logarithmic bins
        do i = 1, nsamples
          distances_k(i) = dk * (kmax/dk)**(real(i,dp)/real(nsamples,dp))
          Pk_counts(i)=0;           
          Pk(i) = 0;
        enddo
      !$OMP ENDDO

      !$OMP DO
        do i = 1, nc/2+1
          ix = i-1
          do j = 1, nc
            iy = j-1
            if (iy>nc/2+1) iy = iy-nc
            do k = 1, nc
              iz = k-1
              if (iz>nc/2+1) iz = iz-nc

              d = sqrt(real(ix**2+iy**2+iz**2,kind=dp))*dk
              if (d<=kmax .and. d>0.d0) then
                ig = int( log(d/dk)*nsamples/log(kmax/dk) ) -1
                ig = max(ig, 1)
                do ik=ig, nsamples
                  if (d<=distances_k(ik)) exit
                enddo
                f = 2
                if (i==1 .or. i==nc/2+1) f=1
                !$OMP ATOMIC
                Pk(ik) = Pk(ik)+f*real(Pk_field(i,j,k))
                !$OMP ATOMIC
                Pk_counts(ik) = Pk_counts(ik) + f
              endif
            enddo
          enddo
        enddo
      !$OMP END DO

    !$OMP END PARALLEL




    !-----------------------------------------------------------------------
    write(*,*) "Computing and histogramming xi(r) FFT for case ", add_name
    !-----------------------------------------------------------------------

    call fftw_execute_dft_c2r(plan_backward, Pk_field, density_field)
    
    !$OMP PARALLEL PRIVATE(i,j,k,ix,iy,iz,ik,ig,d)

      !$OMP DO COLLAPSE(3)
        do i=1, nc
          do j=1, nc
            do k=1, nc
              density_field(i,j,k) = density_field(i,j,k)/sqrtnccube
            enddo
          enddo
        enddo
      !$OMP END DO

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

      !$OMP DO
        ! Get bin distances, reset values
        do i = 1, nsamples
          distances(i) = i*lmax/nsamples
          correlation_counts(i)=0;  
          correlation(i)=0;
        enddo
      !$OMP END DO

      !$OMP DO
        do i = 1, nc
          ix = i-1
          if (ix>nc/2+1) ix = ix-nc
          do j = 1, nc
            iy = j-1
            if (iy>nc/2+1) iy = iy-nc
            do k = 1,nc
              iz = k-1
              if (iz>nc/2+1) iz = iz-nc
              d = sqrt(real(ix**2+iy**2+iz**2,kind=dp))*dx
              if (d<=lmax) then
                ig = int(d*nsamples/lmax) - 1
                ig = max(ig, 1)
                do ik=ig, nsamples
                  if (d<=distances(ik)) exit
                enddo
                !$OMP ATOMIC
                correlation(ik) = correlation(ik)+density_field(i,j,k)
                !$OMP ATOMIC
                correlation_counts(ik) = correlation_counts(ik) + 1
              endif
            enddo
          enddo
        enddo
      !$OMP END DO




      !---------------------------------
      ! Write results to file
      !---------------------------------

!       !$OMP SINGLE
      !   outfname=TRIM(TRIM(srcdir)//"/Pk_"//TRIM(add_name)//".txt")
      !   open(unit=666, file=outfname, form='formatted')
      !   write(666,*) "Results generated by eval_galaxies.f03"
      !   write(666,'(2A22,x)') "k", "P(k)"
      !   do i = 1, nsamples
      !     if (Pk_counts(i)>0) then
      !       write(666, '(2E22.10,x)') 0.5*(distances_k(i-1)+distances_k(i)), Pk(i)/Pk_counts(i)
      !     endif
      !   enddo
      !   write(*,*) "Finished P(k), written results to "//TRIM(outfname)
      !   close(666)
      ! !$OMP END SINGLE
      !
      ! !$OMP SINGLE
      !   outfname=TRIM(TRIM(srcdir)//"/correlation_"//TRIM(add_name)//".txt")
      !   open(unit=667, file=outfname, form='formatted')
      !   write(667,*) "Results generated by eval_galaxies.f03"
      !   write(667,'(2A22,x)') "r", "correlation(r)"
      !   do i = 1, nsamples
      !     if (correlation_counts(i) > 0 ) then
      !       write(667, '(2E22.10,x)') 0.5*(distances(i-1)+distances(i)), correlation(i)/correlation_counts(i)
      !     endif
      !   enddo
      !   write(*,*) "Finished xi(r), written results to "//TRIM(outfname)
      !   close(667)
      ! !$OMP END SINGLE
      ! !$OMP SINGLE
      !   if (check_normalisation) then
      !     write(*,'(A25,3E18.10)') "Normalisation check: R->F", deltasq, delta_ksq, deltasq/delta_ksq
      !     write(*,'(A25,3E18.10)') "Normalisation check: F->R", pksq, corrsq, pksq/corrsq
      !   endif
      ! !$OMP END SINGLE

      !$OMP SINGLE
        outfname=TRIM(TRIM(srcdir)//"/Pk_"//TRIM(add_name)//".txt")
        open(unit=666, file=outfname, form='formatted')
        write(666,*) "Results generated by get_correlation.f03"
        write(666,'(3A22,x)') "k", "P(k)", "Pk(k) counts"
        do i = 1, nsamples
          write(666, '(2(E22.10,x), I22)') distances_k(i), Pk(i), Pk_counts(i)
        enddo
        write(*,*) "Finished P(k), written results to "//TRIM(outfname)
        close(666)
      !$OMP END SINGLE

      !$OMP SINGLE
        outfname=TRIM(TRIM(srcdir)//"/correlation_"//TRIM(add_name)//".txt")
        open(unit=667, file=outfname, form='formatted')
        write(667,*) "Results generated by get_correlation.f03"
        write(667,'(3A22,x)') "r", "correlation(r)", "correlation(r) counts"
        do i = 1, nsamples
          write(667, '(2(E22.10,x), I22)') distances(i), correlation(i), correlation_counts(i)
        enddo
        write(*,*) "Finished xi(r), written results to "//TRIM(outfname)
        close(667)
      !$OMP END SINGLE

      !$OMP SINGLE
        if (check_normalisation) then
          write(*,'(A25,3E18.10)') "Normalisation check: R->F", deltasq, delta_ksq, deltasq/delta_ksq
          write(*,'(A25,3E18.10)') "Normalisation check: F->R", pksq, corrsq, pksq/corrsq
        endif
      !$OMP END SINGLE
    !$OMP END PARALLEL


  enddo


  deallocate(x,y,z,stellar_mass,halo_id,clmp_id)
  deallocate(Pk_field, density_field, distances, distances_k, Pk, correlation, Pk_counts, correlation_counts)
  call fftw_destroy_plan(plan_forward)
  call fftw_destroy_plan(plan_backward)

  


  write(*,*) "eval_galaxies.f03 finished."

contains


  !===========================================
  subroutine get_density_field(which, interp)
  !===========================================
    !------------------------------------------------------------
    ! Computes the density field and overdensity field delta
    ! which = 1: main + satellites
    ! which = 2: main + satellites + orphans
    ! interp = 'ngp' : nearest grid point interpolation
    ! interp = 'cic' : cloud in cell interpolation
    !------------------------------------------------------------

    implicit none
    integer, intent(in)  :: which
    character(len=3)     :: interp
    integer              :: i, j, k, gal
    integer              :: iup, idown, jup, jdown, kup, kdown
    real(dp)             :: rho, xup, yup, zup, hdc, cv
    logical              :: found_it

    real(dp) :: s

    hdc = dc/2
    cv = dc**3

    !$OMP PARALLEL DEFAULT(shared) PRIVATE(i,j,k,gal,found_it,iup,idown,jup,jdown,kup,kdown,rho,xup,yup,zup)

      !$OMP DO COLLAPSE(3)
        ! reset density
        do i=1, nc
          do j=1, nc
            do k=1, nc
              density_field(i,j,k) = 0.d0
            enddo
          enddo
        enddo
      !$OMP END DO

      if (interp=='ngp') then
        !$OMP DO
          do gal=1, ngalaxies
            if (which==2 .or. clmp_id(gal)>0) then
              i = int(x(gal)/boxlen*nc)+1
              j = int(y(gal)/boxlen*nc)+1
              k = int(z(gal)/boxlen*nc)+1
              !$OMP ATOMIC
              density_field(i,j,k) = density_field(i,j,k)+stellar_mass(gal)
            endif
          enddo
        !$OMP END DO

      else if (interp=='cic') then
        !$OMP DO
          do gal=1, ngalaxies
            if (which==2 .or. clmp_id(gal)>0) then
              iup   = int((x(gal)+hdc)/boxlen*nc)+1
              idown = iup-1
              jup   = int((y(gal)+hdc)/boxlen*nc)+1
              jdown = jup-1
              kup   = int((z(gal)+hdc)/boxlen*nc)+1
              kdown = kup-1

              rho = stellar_mass(gal)/cv
              xup = x(gal) + hdc - (iup-1)*dc
              yup = y(gal) + hdc - (jup-1)*dc
              zup = z(gal) + hdc - (kup-1)*dc

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

            endif
          enddo
        !$OMP END DO

      else
        write(*,*) "Didn't recognize interpolation method ", interp
        stop
      endif
        

      !$OMP DO COLLAPSE(3)
        do i=1, nc
          do j=1, nc
            do k=1, nc
              density_field(i,j,k) = density_field(i,j,k)/volume
            enddo
          enddo
        enddo
      !$OMP END DO
    !$OMP END PARALLEL

    s = 0
    do gal=1, ngalaxies
      if (clmp_id(gal)>0) s = s + stellar_mass(gal)
    enddo
    write(*,'(3(A15,E18.10))') "TOTAL MASS:", sum(density_field)*unit_l**3, "SHOULD BE", sum(stellar_mass), "or", s

  end subroutine get_density_field





  !========================================
  subroutine read_cmdlineargs()
  !========================================
    !-----------------------
    ! reads cmd line args
    !-----------------------
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
        if (ncells_temp>=0) nc = ncells_temp
      endif
    end do
  end subroutine read_cmdlineargs



  
  !==================================
  subroutine read_data()
  !==================================
    !------------------------------------------
    ! reads in file sizes so you can allocate
    ! arrays of appropriate sizes.
    ! Essentially, get nhaloes and ngalaxies
    !------------------------------------------

    use omp_lib
    character(len=5)            :: outputnr, cpunr
    character(len=80)           :: fname
    character(len=125)          :: junk
    integer, dimension(0:ncpu)  :: counts_g, counts_h
    integer                     :: i, ind, id, io, ifile
    logical                     :: file_exists

    counts_g = 0; counts_g(0)=1;
    counts_h = 0; counts_h(0)=1;
    outputnr = srcdir(8:12)

    !$OMP PARALLEL PRIVATE(i, id, io, ifile, file_exists, fname, cpunr, ind, junk)
      id = omp_get_thread_num()

      !----------------------------------
      ! Read file sizes
      !----------------------------------

      !$OMP DO
        do ifile=1, ncpu
          call title(ifile, cpunr)
          fname=TRIM(srcdir//'/galaxies_'//outputnr//'.txt'//cpunr)
          inquire(file=fname, exist=file_exists)
          if (file_exists) then
            open(unit=600+id, file=fname)
            read(600+id, *) ! skip header
            do 
              read(600+id,*, iostat=io)
              if (io/=0) exit
              counts_g(ifile) = counts_g(ifile)+1
            enddo
            close(600+id)
          else
            write(*,*) "Didn't find file ", fname
          endif
        enddo
      !$OMP ENDDO NOWAIT

      !$OMP DO
        do ifile=1, ncpu
          call title(ifile, cpunr)
          fname=TRIM(srcdir//'/halo_'//outputnr//'.txt'//cpunr)
          inquire(file=fname, exist=file_exists)
          if (file_exists) then
            open(unit=600+id, file=fname)
            read(600+id,*) ! skip header
            do 
              read(600+id,*, iostat=io)
              if (io/=0) exit
              counts_h(ifile) = counts_h(ifile)+1
            enddo
            close(600+id)
          else
            write(*,*) "Didn't find file ", fname
          endif
        enddo
      !$OMP ENDDO

    
      !$OMP DO ORDERED
        do i=0, ncpu-1
          !$OMP ORDERED
            counts_g(i+1) = counts_g(i+1)+counts_g(i)
            counts_h(i+1) = counts_h(i+1)+counts_h(i)
          !$OMP END ORDERED
        enddo
      !$OMP ENDDO



      !----------------------------------
      ! Allocate Arrays
      !----------------------------------

      !$OMP SINGLE
        nhalos = counts_h(ncpu)-1
        ngalaxies = counts_g(ncpu)-1
        allocate(x(1:ngalaxies));             x=0;
        allocate(y(1:ngalaxies));             y=0;
        allocate(z(1:ngalaxies));             z=0;
        allocate(stellar_mass(1:ngalaxies));  stellar_mass=0;
        allocate(clmp_id(1:ngalaxies));       clmp_id=0;
        allocate(halo_id(1:nhalos));          halo_id=0;
      !$OMP END SINGLE
      !$OMP BARRIER



      !----------------------------------
      ! Read galaxy data
      !----------------------------------

      !$OMP DO
        do ifile=1, ncpu
          call title(ifile, cpunr)
          fname=TRIM(srcdir//'/galaxies_'//outputnr//'.txt'//cpunr)
          open(unit=600+id, file=fname)
          read(600+id,*) ! skip header
          do ind=counts_g(ifile-1), counts_g(ifile)-1
            read(600+id,'(I20,x,4(E20.12,x))') clmp_id(ind), stellar_mass(ind), x(ind), y(ind), z(ind)
          enddo
          close(600+id)
        enddo
      !$OMP END DO


      !----------------------------------
      ! Read halo data
      !----------------------------------

      !$OMP DO
        do ifile=1, ncpu
          call title(ifile, cpunr)
          fname=TRIM(srcdir//'/halo_'//outputnr//'.txt'//cpunr)
          open(unit=600+id, file=fname)
          read(600+id,*) ! skip header
          do ind=counts_h(ifile-1), counts_h(ifile)-1
            read(600+id,'(I10,x,A125)') halo_id(ind), junk
          enddo
          close(600+id)
        enddo
      !$OMP END DO

    !$OMP END PARALLEL
  end subroutine read_data





  !========================================
  subroutine read_info()
  !========================================
    !-------------------------------
    ! Read in galaxy data.
    !-------------------------------
    implicit none
    character(len=5)  :: outputnr
    character(len=13) :: junk
    character(len=80) :: fname
    logical :: file_exists
    integer :: i

    outputnr = srcdir(8:12)
    fname=TRIM(srcdir//'/info_'//outputnr//'.txt')
    
    inquire(file=fname, exist=file_exists)
    if (file_exists) then

      open(unit=666, file=fname)
      read(666, '(A13,I11)') junk, ncpu
      do i = 1, 8
        ! skip lines
        read(666,*)
      enddo

      read(666, '(A13,E23.15)') junk, aexp
      read(666, '(A13,E23.15)') junk, H0

      do i = 1, 4
        ! skip lines
        read(666,*)
      enddo

      read(666, '(A13,E23.15)') junk, unit_l
      unit_l = unit_l / Mpc

      close(666)

    else
      write(*,*) "Didn't find file", fname
      stop
    endif
  end subroutine read_info



  !========================================
  subroutine smf()
  !========================================
    !-----------------------------------------------
    ! Computes and writes stellar mass function
    !-----------------------------------------------

    ! log10(M_galaxy) boundaries:
    real(dp), parameter                 :: minmass = 5.5d0
    real(dp), parameter                 :: maxmass = 13.0d0
    real(dp), allocatable, dimension(:) :: smf_all, smf_sub, smf_main, logmass
    real(dp)                            :: dm, lm, f
    integer                             :: i,j,guess
    character(len=80)                   :: fname


    allocate(smf_all(1:nsamples));  smf_all =0;
    allocate(smf_sub(1:nsamples));  smf_sub =0;
    allocate(smf_main(1:nsamples)); smf_main=0;
    allocate(logmass(0:nsamples));  logmass =0;


    dm = (maxmass-minmass)/nsamples


    !$OMP PARALLEL PRIVATE(i,j,guess,lm,f)
      !$OMP DO
        do i = 0, nsamples 
          logmass(i) = minmass+i*dm
        enddo
      !$OMP ENDDO

      !$OMP DO
        do i = 1, ngalaxies
          ! determine index
          lm = log10(stellar_mass(i))
          guess = int((lm-minmass)/dm) - 1
          guess = min(guess, 1)
          do
            if (logmass(guess)>lm) then
              smf_all(guess) = smf_all(guess) + 1
              exit
            endif
            guess = guess + 1
          enddo
          
          ! check if subhalo or main
          !$OMP ATOMIC
          smf_all(guess) = smf_all(guess) + 1
          if (clmp_id(i) > 0) then
            ! found non-orphan
            !$OMP ATOMIC
            smf_sub(guess) = smf_sub(guess) + 1
            do j=1, nhalos
              if (clmp_id(i) == halo_id(j)) then
                ! found main
                !$OMP ATOMIC
                smf_main(guess) = smf_main(guess) + 1
              endif
            enddo
          endif
        enddo
      !$OMP ENDDO

      !$OMP DO
        do i=1, nsamples
          f = (logmass(i)-logmass(i-1))*volume
          smf_all(i) = smf_all(i)/f
          smf_sub(i) = smf_sub(i)/f
          smf_main(i) = smf_main(i)/f
        enddo
      !$OMP ENDDO

    !$OMP END PARALLEL


    fname=TRIM(TRIM(srcdir)//"/smf.txt")
    open(unit=666, file=fname, form='formatted')
    write(666, '(4A20)') "Mass bins", "SMF subhalos+orphans", "SMF subhalos", "SMF centrals only"
    do i=1, nsamples
      write(666, '(F20.6, 3E20.12)') (logmass(i-1)+logmass(i))*0.5, smf_all(i), smf_sub(i), smf_main(i)
    enddo
    close(666)

    write(*,*) "SMF stuff finished, results written to ", fname


    deallocate(smf_main, smf_all, smf_sub, logmass)
  end subroutine smf





  !==================================
  subroutine title(n,nchar)
  !==================================
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





  !========================================
  subroutine write_galaxy_images()
  !========================================
    !--------------------------------------------------------
    ! write images to file that can be
    ! plotted directly with plot_fortran_galaxies.py
    !--------------------------------------------------------

    real(dp), dimension(1:nc, 1:nc) :: density_image
    real(dp) :: pixel
    integer :: i, j, k

    character(len=80) :: fname

    write(*,*) "Creating density field image for plotting."

    !$OMP DO
    do i=1, nc
      do j=1, nc
        pixel = 0.d0
        do k = 1, nc
          pixel = pixel + density_field(i, j, k) 
        enddo
        density_image(i,j) = pixel
      enddo
    enddo
    !$OMP END DO


    write(*,*) "Writing density field image to file."

    fname = TRIM('density_image.txt')
    open(unit=666, file=fname, form='formatted')
    ! note switched j/i for pyplot.imshow
    do j=1, nc
      do i=1, nc
        write(666, '(E14.7)', advance='no') density_image(i,j)
      enddo
      write(666,*) ! write newline
    enddo
    close(666)

    write(*,*) fname, " written."
 
    fname = TRIM('galaxy_field.txt')
    open(unit=666, file=fname, form='formatted')
    do i=1, ngalaxies
      write(666,*) x(i), y(i)
    enddo
    close(666)

    write(*,*) fname, " written."
  end subroutine write_galaxy_images 

end program
