program eval_galaxies

  !======================================
  ! Compute SMF for centrals, satellites
  ! and orphans individually
  !======================================

  use omp_lib
  implicit none

  integer, parameter  :: dp = kind(1.0d0) 

  !----------------------
  ! Manual parameters
  !----------------------
  integer             :: output_start = 39          ! first snapshot to use
  integer             :: output_end   = 41          ! last snapshot to use
  integer             :: nc           = 2000        ! number of cells for density field
  integer, parameter  :: nsamples     = 100         ! number of samples for output; Number of bins for histograms
  character(len=3)    :: interpolation = 'cic'      ! which interpolation to use: cic for cloud in cell, ngp for nearest grid point
  ! character(len=3)    :: interpolation = 'ngp'    ! which interpolation to use: cic for cloud in cell, ngp for nearest grid point
  real(dp), parameter :: minmass = 4.d0
  real(dp), parameter :: maxmass = 13.0d0
  logical             :: calculate_masses = .false. ! set =.true. if input is halo masses instead of galaxy masses
  character(len=4)    :: model = 'behr'             ! if calculate_masses: which model to use: behroozi 2013
  ! character(len=4)    :: model = 'most'           ! if calculate_masses: which model to use: moster 2013
  logical             :: pic = .false.              ! whether to create density projection image


  !------------------
  ! Constants
  !------------------
  real(dp), parameter :: Mpc = 3.086d+24            ! Mpc in cm
  real(dp), parameter :: M_Sol = 1.998d33           ! solar mass in g
  real(dp), parameter :: boxlen = 1.d0


  !-----------------------
  ! Read-in data
  !-----------------------
  real(dp)            :: aexp                       ! expansion factor
  real(dp)            :: unit_l                     ! length of box [Mpc]
  real(dp)            :: volume                     ! current physical volume
  character(len=12)   :: srcdir                     ! output_XXXXX 
  integer             :: nhalos, ngalaxies, ncpu
  logical             :: skip_this                  ! whether to skip empty output

  
  !--------------------
  ! Grid/Plot stuff
  !--------------------
  real(dp)            :: dc                         ! cell size in code units



  !---------------------
  ! Arrays
  !---------------------
  ! read-in data arrays
  real(dp), dimension(:), allocatable :: x, y, z              ! x,y,z positions of galaxies [code units]
  real(dp), dimension(:), allocatable :: stellar_mass         ! stellar mass of galaxies [M_Sol]
  integer,  dimension(:), allocatable :: halo_id, clmp_id     ! ID's of halos in snapshot and associated clump IDs

  ! computed global arrays
  real(dp), allocatable, dimension(:)      :: logmass         ! log10(stellar mass) bins
  real(dp), allocatable, dimension(:,:)    :: density_field   ! density field
  character(len=3), dimension(1:2)         :: cases=(/'sub','all'/)



  !-----------------------------
  ! Other variables
  !-----------------------------
  integer                          :: i, outp, c
  real(dp)                         :: dm
  character(len=5)                 :: dirnr

  logical  :: seed_set = .false.
  




  !----------------------------
  ! Preparations
  !----------------------------

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
  
  dc = 1.d0/nc





  write(*,*) "============================"
  write(*,*) "get_smf.f03 started."
  write(*,*) "============================"
  write(*,*)                "Working with parameters:"
  write(*,'(A20,I25)')      "  ncpu:          ", ncpu
  write(*,'(A20,I25)')      "  nc:            ", nc
  write(*,'(A20,E25.8)')    "  log10(M_S)_min:", minmass
  write(*,'(A20,E25.8)')    "  log10(M_S)_max:", maxmass
  write(*,'(A20,A25)')      "  interpolation: ", interpolation
  if (calculate_masses) then
  write(*,'(A20,A25)')      "  SM model:      ", model
  endif
  if (pic) then
  write(*,'(A20,A25)')      "  make picture:  ", "yes"
  endif




  !--------------------------
  ! MAIN LOOP
  !--------------------------

  do outp = output_start, output_end

    ! Reset stuff
    skip_this = .false.
    nhalos = 0
    ngalaxies = 0

    call title(outp, dirnr)
    srcdir = TRIM('output_'//dirnr)

    call read_info()        ! read data from info files
    call read_data()        ! read actual data and nhalos, ngalaxies
    write(*,'(A60)')          "-------------------------------------------------------------"
    write(*,'(A20,A25)')      "  srcdir:        ", srcdir
    write(*,'(A20,I25)')      "  nhaloes:       ", nhalos
    write(*,'(A20,I25)')      "  ngalaxies:     ", ngalaxies
    write(*,'(A20,F25.3)')    "  unit_l:        ", unit_l
    write(*,'(A20,F25.3)')    "  a_exp:         ", aexp
    write(*,'(A20,F25.3)')    "  z:             ", 1.d0/aexp - 1.d0

    if (skip_this) then
      write(*,'(A45)') "   No data found. Skipping this snapshot."
      cycle
    endif


    !-----------------------------------------------
    write(*,*) "Calculating Stellar Mass Function"
    !-----------------------------------------------

    if (calculate_masses) call calculate_sm()
    call smf()


    if (pic) then 
      volume = unit_l**3
      allocate(density_field(1:nc, 1:nc))
      do c=1, 2
        write(*,'(A33,x,A3)') " Computing density field for case ", cases(c)
        call get_density_field(c, interpolation)
        call write_galaxy_images(cases(c))
      enddo
      deallocate(density_field)
    endif

    deallocate(x,y,z,stellar_mass,halo_id,clmp_id)

  enddo

 

  write(*,*) "get_smf.f03 finished."
  deallocate(logmass)

contains


  !===================================
  subroutine calculate_sm()
  !===================================
    !------------------------------------------
    ! Subroutine to compute the stellar mass
    ! assuming the read-in 'stellar_mass' array
    ! doesn't actually contain galaxy masses,
    ! but halo masses to be used for the SMHM 
    ! relation.
    !------------------------------------------
    
    implicit none
    integer :: i

    if (model=='behr') then
      !$OMP PARALLEL DO PRIVATE(i)
        do i=1, ngalaxies
          stellar_mass(i) = smf_behroozi(stellar_mass(i), aexp)
        enddo
      !$OMP END PARALLEL DO

    else if (model=='most') then
      !$OMP PARALLEL DO PRIVATE(i)
        do i=1, ngalaxies
          stellar_mass(i) = smf_moster(stellar_mass(i), aexp)
        enddo
      !$OMP END PARALLEL DO
    endif
  end subroutine calculate_sm




  !===================================
  real(dp) function smf_behroozi(m,a)
  !===================================
    !-----------------------------------------------------------
    ! Computes stellar mass given peak mass and expansion
    ! factor a at time when clump had peak mass using a 
    ! parametric SHAM relation taken from
    ! Behroozi, Wechsler & Conroy 2013
    ! DOI:	                10.1088/0004-637X/770/1/57
    ! Bibliographic Code:	  2013ApJ...770...57B
    ! http://adsabs.harvard.edu/abs/2013ApJ...770...57B 
    !-----------------------------------------------------------
    
    use omp_lib

    implicit none
    real(dp), intent(in) :: m,a ! mass

    real(dp) :: euler = 2.7182818284590
    real(dp) :: log10_2 = 0.301029995663981 ! log_10(2)

    real(dp) :: M_10    =   11.514 
    real(dp) :: M_1a    = -  1.793
    real(dp) :: M_1z    = -  0.251
    real(dp) :: e_0     = -  1.777
    real(dp) :: e_a     = -  0.006
    real(dp) :: e_z     =    0.000
    real(dp) :: e_a2    = -  0.119
    real(dp) :: alpha_0 = -  1.412
    real(dp) :: alpha_a =    0.731
    real(dp) :: delta_0 =    3.508 
    real(dp) :: delta_a =    2.608
    real(dp) :: delta_z = -  0.043
    real(dp) :: gamma_0 =    0.316 
    real(dp) :: gamma_a =    1.319 
    real(dp) :: gamma_z =    0.279
    real(dp) :: xi_0    =    0.218
    real(dp) :: xi_a    = -  0.023

    real(dp) :: nu
    real(dp) :: loge
    real(dp) :: alpha
    real(dp) :: delta
    real(dp) :: gam
    real(dp) :: xi, sig
    real(dp) :: logM1

    real(dp) :: z, f0

    integer :: n, clock, i
    integer, dimension(:), allocatable:: seed

    z = 1.d0/a - 1.d0

    nu = exp(-4.d0*a**2)
    logM1 = M_10    + nu*(M_1a   *(a-1)  + M_1z*z)
    loge  = e_0     + nu*(e_a    *(a-1)  + e_z*z ) + e_a2*(a-1)
    alpha = alpha_0 + nu*(alpha_a*(a-1))
    delta = delta_0 + nu*(delta_a*(a-1)            + delta_z*z)
    gam   = gamma_0 + nu*(gamma_a*(a-1)            + gamma_z*z)

    !-----------------
    ! get scatter
    !-----------------
  
    ! set unique seed for every task if necessary
    if (.not.seed_set) then
      call random_seed(size=n)
      allocate(seed(n))
      call system_clock(count=clock)
      seed = clock + omp_get_thread_num() + 53*[(n-i, i=1, n)]
      call random_seed(put=seed)
      deallocate(seed)
      seed_set = .true.
    endif


    call random_number(xi)      ! temporarily store random number in xi
    call random_number(sig)     ! get sign
    if (sig > 0.5d0) then
      sig = -1.d0
    else
      sig = 1.d0
    endif

    xi = sig*xi*(xi_0 + xi_a * (a-1))


    !-------------------------------
    ! finally, get stellar mass
    !-------------------------------
    f0 = -log10_2 + delta*(log10_2**gam)/(1.d0 + euler)
    smf_behroozi = 10.d0**(loge+logM1 + behr_f(log10(m)-logM1, alpha, delta, gam)- f0 + xi)
  end function smf_behroozi



  !==================================================
  real(dp) function behr_f(x, alpha, delta, gam)
  !==================================================
    implicit none
    real(dp), intent(in) :: x, alpha, delta, gam
    behr_f = -log10(10.d0**(alpha*x) + 1) + delta*log10(1.d0+exp(x))**gam/(1.d0 + exp(10.d0**(-x)))
  end function behr_f





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
    integer              :: i, j, gal
    integer              :: iup, idown, jup, jdown, kup, kdown
    real(dp)             :: rho, xup, yup, zup, hdc, cv
    logical              :: found_it


    hdc = dc/2
    cv = dc**3

    !$OMP PARALLEL DEFAULT(shared) PRIVATE(i,j,gal,found_it,iup,idown,jup,jdown,kup,kdown,rho,xup,yup,zup)
      !$OMP DO COLLAPSE(2)
        ! reset density
        do i=1, nc
          do j=1, nc
            density_field(i,j) = 0.d0
          enddo
        enddo
      !$OMP END DO

      if (interp=='ngp') then
        !$OMP DO
          do gal=1, ngalaxies
            if (which==2 .or. clmp_id(gal)>0) then
              i = int(x(gal)/boxlen*nc)+1
              j = int(y(gal)/boxlen*nc)+1
              !$OMP ATOMIC
              density_field(i,j) = density_field(i,j)+stellar_mass(gal)
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

              rho = stellar_mass(gal)/cv
              xup = x(gal) + hdc - (iup-1)*dc
              yup = y(gal) + hdc - (jup-1)*dc

              if (iup>nc)  iup   = iup-nc
              if (idown<1) idown = nc+idown
              if (jup>nc)  jup   = jup-nc
              if (jdown<1) jdown = nc+jdown

              !$OMP CRITICAL
                density_field(iup,   jup  )   = density_field(iup,   jup  )  + xup      * yup      * rho
                density_field(idown, jup  )   = density_field(idown, jup  )  + (dc-xup) * yup      * rho
                density_field(iup,   jdown)   = density_field(iup,   jdown)  + xup      * (dc-yup) * rho
                density_field(idown, jdown)   = density_field(idown, jdown)  + (dc-xup) * (dc-yup) * rho
              !$OMP END CRITICAL

            endif
          enddo
        !$OMP END DO

      else
        write(*,*) "Didn't recognize interpolation method ", interp
        stop
      endif
        

      !$OMP DO COLLAPSE(2)
        do i=1, nc
          do j=1, nc
            density_field(i,j) = density_field(i,j)/volume
          enddo
        enddo
      !$OMP END DO
    !$OMP END PARALLEL
  end subroutine get_density_field






  !=======================================
  real(dp) function smf_moster(m,a)
  !=======================================
    !-----------------------------------------------------------
    ! Computes stellar mass given peak mass and expansion
    ! factor a at time when clump had peak mass using a 
    ! parametric SHAM relation taken from
    ! Benjamin P. Moster, Thorsten Naab, Simon D. M. White
    ! DOI: 	                10.1093/mnras/sts261
    ! arXiv:                1205.5807
    ! IMPORTANT NOTE: moster et al 2013
    ! Actually use m= M(m_inf, a_inf) for subhaloes
    ! But just to check smfs of central galaxies, 
    ! it should be fine
    !-----------------------------------------------------------

    implicit none
    real(dp), intent(in) :: m,a

    real(dp) :: logM1, M1
    real(dp) :: N   
    real(dp) :: beta
    real(dp) :: gam 
    real(dp) :: M_10     =  11.590
    real(dp) :: M_11     =   1.195
    real(dp) :: N_10     =   0.0351
    real(dp) :: N_11     = - 0.0247
    real(dp) :: beta_10  =   1.376
    real(dp) :: beta_11  = - 0.826
    real(dp) :: gamma_10 =   0.608
    real(dp) :: gamma_11 =   0.329
    real(dp) :: am1

    am1 = a - 1

    logM1 = M_10     + M_11*am1
    M1    = 10**logM1
    N     = N_10     + N_11*am1
    beta  = beta_10  + beta_11*am1
    gam   = gamma_10 + gamma_11*am1

    smf_moster = m*2.d0*N/((m/M1)**(-beta)+(m/M1)**gam)
  end function smf_moster





  
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
    integer                     :: i, ind, id, io, ifile, s1, s2
    logical                     :: file_exists

    counts_g = 0; counts_g(0)=1;
    counts_h = 0; counts_h(0)=1;
    outputnr = srcdir(8:12)

    s1=0
    s2=0

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
          ! else
          !   write(*,*) "Didn't find file ", fname
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
          ! else
          !   write(*,*) "Didn't find file ", fname
          endif
        enddo
      !$OMP ENDDO

      !$OMP DO REDUCTION(+:s1, s2)
        do ifile=0, ncpu
          s1 = s1 + counts_g(ifile)
          s2 = s2 + counts_h(ifile)
        enddo
      !$OMP ENDDO

      if (s1==1 .or. s2==1) then
        !$OMP MASTER
          if (s1==1) then
            write(*,*) "No galaxy files found."
          endif
          if (s2==1) then
            write(*,*) "No halo files found."
          endif
        skip_this = .true.
        !$OMP END MASTER
      endif
    !$OMP END PARALLEL

    if (skip_this) return


    do i=0, ncpu-1
        counts_g(i+1) = counts_g(i+1)+counts_g(i)
        counts_h(i+1) = counts_h(i+1)+counts_h(i)
    enddo



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

    !$OMP PARALLEL PRIVATE(i, id, io, ifile, file_exists, fname, cpunr, ind, junk)
      id = omp_get_thread_num()
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
      !$OMP END DO NOWAIT


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
    !--------------------------------------
    ! Read in galaxy snapshot info data.
    !--------------------------------------
    implicit none
    character(len=5)   :: outnr_str
    character(len=13)  :: junk
    character(len=80)  :: fname
    logical            :: file_exists
    integer            :: i

    outnr_str = srcdir(8:12)
    fname=TRIM(srcdir//'/info_'//outnr_str//'.txt')
    
    inquire(file=fname, exist=file_exists)
    if (file_exists) then

      open(unit=666, file=fname)
      read(666, '(A13,I11)') junk, ncpu
      do i = 1, 8
        ! skip lines
        read(666,*)
      enddo

      read(666, '(A13,E23.15)') junk, aexp

      do i = 1, 5
        ! skip lines
        read(666,*)
      enddo

      read(666, '(A13,E23.15)') junk, unit_l
      unit_l = unit_l / Mpc / aexp

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

    real(dp), allocatable, dimension(:) :: smf_main, smf_sub, smf_orph
    real(dp)                            :: lm, f
    integer                             :: i,j,guess
    character(len=80)                   :: fname
    logical                             :: found


    allocate(smf_main(0:nsamples)); smf_main=0;
    allocate(smf_sub(0:nsamples)); smf_sub=0;
    allocate(smf_orph(0:nsamples)); smf_orph=0;


    !$OMP PARALLEL PRIVATE(i,j,guess,lm,f, found)
      !$OMP DO
        do i = 1, ngalaxies
          ! determine index
          lm = log10(stellar_mass(i))
          guess = int((lm-minmass)/dm) + 1
          if (guess<1) then
            write(*,'(A9,x,E12.5, A42, E12.5)') "Got mass", lm, " which is smaller than your set minmass of", minmass
            write(*,*) "Consider lowering the threshold."
            cycle
          endif

          if (clmp_id(i) > 0) then
            ! found non-orphan
            found = .false.
            do j=1, nhalos
              if (clmp_id(i) == halo_id(j)) then
                ! found main
                !$OMP ATOMIC
                smf_main(guess) = smf_main(guess) + 1
                found = .true.
                exit
              endif
            enddo
            if (.not.found) then
              ! if you didn't exit at this point, you have a subhalo
              !$OMP ATOMIC
              smf_sub(guess) = smf_sub(guess) + 1
            endif
          else
            !$OMP ATOMIC
            smf_orph(guess) = smf_orph(guess) + 1
          endif

        enddo
      !$OMP ENDDO

      ! !$OMP DO
      !   do i=1, nsamples
      !     f = (logmass(i)-logmass(i-1))*volume
      !     smf_all(i) = smf_all(i)/f
      !     smf_sub(i) = smf_sub(i)/f
      !     smf_main(i) = smf_main(i)/f
      !   enddo
      ! !$OMP ENDDO

    !$OMP END PARALLEL


    fname=TRIM(TRIM(srcdir)//"/smf-all.txt")
    open(unit=666, file=fname, form='formatted')
    write(666, '(4A20)') "Mass bins", "SMF counts main", "SMF counts sub", "SMF counts orph"
    do i=0, nsamples
      write(666, '(F20.6, 3E20.12)') logmass(i), smf_main(i), smf_sub(i), smf_orph(i)
    enddo
    close(666)

    write(*,*) "SMF stuff finished, results written to ", fname
    write(*,*) "Found", sum(smf_main), "central galaxies"
    write(*,*) "     ", sum(smf_sub), "satellite galaxies"
    write(*,*) "     ", sum(smf_orph), "orphan galaxies"


    deallocate(smf_main)
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
  subroutine write_galaxy_images(add_name)
  !========================================
    !--------------------------------------------------------
    ! write images to file that can be
    ! plotted directly with plot_fortran_galaxies.py
    !--------------------------------------------------------

    character(len=3), intent(in) :: add_name
    character(len=80)            :: fname


    fname = TRIM(srcdir//'/density_image-'//add_name//'.dat')
    open(unit=666, file=fname, form='unformatted')
    write(666) nc
    write(666) aexp
    write(666) density_field
    close(666)

    write(*,'(A14,A)') " Written file ", TRIM(fname)
  end subroutine write_galaxy_images 

end program
