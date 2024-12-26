!-------------------------------------------------------------
! Module containing common routines to read in ramses data
!-------------------------------------------------------------


module io_module

  use constants_and_parameters

  implicit none

  ! Global parameters
  !----------------------------

  character(len=12)   :: srcdir                        ! output_XXXXX


  ! Infofile parameters
  !-------------------------

  integer   :: ncpu
  integer   :: ndim
  integer   :: levelmin
  integer   :: levelmax
  integer   :: ngridmax
  integer   :: nstep_coarse

  real(dp)            :: boxlen     ! box size in internal units
  real(dp)            :: time       ! Current time
  real(dp)            :: aexp       ! expansion factor
  real(dp)            :: H0         ! Hubble parameter
  real(dp)            :: unit_l     ! lenght unit
  real(dp)            :: unit_rho   ! density unit
  real(dp)            :: unit_t     ! time unit
  real(dp)            :: Omega_m, Omega_l, Omega_k, Omega_b ! Omega parameters

  ! derived quantities
  real(dp)            :: unit_l_Mpc ! length unit in Mpc
  real(dp)            :: unit_l_cMpc! length unit in COMOVING Mpc
  real(dp)            :: unit_m     ! mass unit
  real(dp)            :: unit_m_MSol! mass unit in M_Sol
  real(dp)            :: volume     ! simulation box volume
  real(dp)            :: volume_Mpc ! simulation box volume in Mpc^3
  real(dp)            :: volume_cMpc! simulation box volume in COMOVING Mpc^3



  ! Global Arrays and Related Integers
  !-------------------------------------

  integer  :: nclumps       ! number of clumps in snapshot
  integer,  dimension(:), allocatable   :: clmp_id ! ID's of halos in snapshot and associated clump IDs
  integer,  dimension(:), allocatable   :: parent_id ! ID's of a clump's parent

  integer  :: ngalaxies     ! number of galaxies in a snapshot
  integer  :: norphans      ! number of orphans in a snapshot
  integer  :: ngalaxies_tot ! total galaxies in snapshot
  real(dp), dimension(:,:), allocatable :: xg     ! positions of galaxies [code units]
  real(dp), dimension(:,:), allocatable :: xo     ! positions of orphans [code units]
  real(dp), dimension(:), allocatable   :: mg     ! masses of galaxies [M_Sol]
  real(dp), dimension(:), allocatable   :: mo     ! masses of orphans [M_Sol]
  integer(i8),dimension(:), allocatable :: ido    ! paritlce id of ophans
  integer , dimension(:), allocatable   :: gal_id ! galaxy associated clump IDs
  logical , dimension(:), allocatable   :: is_central ! is galaxy a central galaxy?



  ! Other Read-In Related Data
  !------------------------------------

  real(dp), dimension(1:3) :: CoM      ! Center of Mass [code units]
  real(dp), dimension(1:3) :: peakpos  ! Peak Position [code units]
  real(dp) :: halomass                 ! (exclusive) clump mass for given halo 



  ! Keeping track of what's been done
  !----------------------------------
  logical :: read_clumps = .false.
  logical :: read_galaxies = .false.
  logical :: infofile_read = .false.


contains


  subroutine get_srcdir_and_halo_from_cmdlineargs(halo)

    !-----------------------------------------------------------------  
    ! reads halo ID and output_XXXXX directory from cmdline args.
    !-----------------------------------------------------------------

    implicit none

    integer, intent(out) :: halo

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
        read(arg, *) halo_temp
        if (halo_temp>=0) halo = halo_temp
      else
        exit
      endif
    end do
  end subroutine get_srcdir_and_halo_from_cmdlineargs


  subroutine get_srcdir_from_cmdlineargs()

    !-----------------------------------------------------------------  
    ! reads output_XXXXX directory from cmdline args.
    !-----------------------------------------------------------------

    implicit none

    character(len=12) :: arg
    character(len=5)  :: dirnr_str
    integer           :: dirnr
    logical           :: exists

    call getarg(1, arg)
    read(arg, *) dirnr
    call title(dirnr, dirnr_str)
    srcdir = TRIM('output_'//dirnr_str)
    inquire(file=srcdir, exist=exists)
    if (.not. exists) then
      write(*,*) "Couldn't find directory ", srcdir, &
                  "pass me only the integer number of the directory"
      stop
    endif

  end subroutine get_srcdir_from_cmdlineargs





  subroutine read_info()

    !-------------------------------------------------------
    ! Read in info_*txt file.
    !-------------------------------------------------------

    implicit none

    character(len=5)   :: outnr_str
    character(len=13)  :: junk
    character(len=80)  :: fname
    logical            :: file_exists

    outnr_str = srcdir(8:12)
    fname=TRIM(srcdir//'/info_'//outnr_str//'.txt')
    
    inquire(file=fname, exist=file_exists)
    if (file_exists) then

      open(unit=666, file=fname)
      read(666, '(A13,I11)') junk, ncpu
      read(666, '(A13,I11)') junk, ndim
      read(666, '(A13,I11)') junk, levelmin
      read(666, '(A13,I11)') junk, levelmax
      read(666, '(A13,I11)') junk, ngridmax
      read(666, '(A13,I11)') junk, nstep_coarse
      read(666,*)! skip empty line
      read(666, '(A13,E23.15)') junk, boxlen
      read(666, '(A13,E23.15)') junk, time
      read(666, '(A13,E23.15)') junk, aexp
      read(666, '(A13,E23.15)') junk, H0
      read(666, '(A13,E23.15)') junk, Omega_m
      read(666, '(A13,E23.15)') junk, Omega_l
      read(666, '(A13,E23.15)') junk, Omega_k
      read(666, '(A13,E23.15)') junk, Omega_b
      read(666, '(A13,E23.15)') junk, unit_l
      read(666, '(A13,E23.15)') junk, unit_rho
      read(666, '(A13,E23.15)') junk, unit_t

      close(666)

    else
      write(*,*) "Didn't find file ", fname
      stop
    endif

    unit_l_Mpc = unit_l / Mpc
    unit_l_cMpc = unit_l / Mpc / aexp
    unit_m = unit_rho * unit_l**3
    unit_m_MSol = unit_m / M_Sol
    volume = (boxlen * unit_l) ** 3
    volume_Mpc = (boxlen * unit_l_Mpc) ** 3
    volume_cMpc = (boxlen * unit_l_cMpc) ** 3

    infofile_read = .true.

  end subroutine read_info




  subroutine read_clump_data(halo)
    !----------------------------------------------------
    ! Read in clump data: clump IDs and their parent IDs.
    ! If halo is positive integer, it will also find
    ! the peak position of the given halo and the clump
    ! (exclusive) mass of the halo, if it exists.
    !----------------------------------------------------

    use omp_lib
    implicit none
    integer, intent(in)         :: halo

    character(len=5)            :: outputnr, cpunr
    character(len=80)           :: fname
    character(len=18)           :: junk_char
    integer, dimension(0:ncpu)  :: counts_c
    integer                     :: i, ind, id, io, ifile, s1, junk1, junk2
    logical                     :: file_exists, skip_this
    real(dp)                    :: xpeak, ypeak, zpeak, clumpmass

    counts_c = 0; counts_c(0)=1;
    outputnr = srcdir(8:12)

    s1=0
    skip_this = .false.

    !$OMP PARALLEL PRIVATE(i, id, io, ifile, file_exists, fname, cpunr, ind) &
    !$OMP DEFAULT(shared)
      id = omp_get_thread_num()

      !----------------------------------
      ! Read file sizes
      !----------------------------------

      !$OMP DO
        do ifile=1, ncpu
          call title(ifile, cpunr)
          fname=TRIM(srcdir//'/clump_'//outputnr//'.txt'//cpunr)
          inquire(file=fname, exist=file_exists)
          if (file_exists) then
            open(unit=600+id, file=fname)
            read(600+id, *) ! skip header
            do 
              read(600+id,*, iostat=io)
              if (io/=0) exit
              counts_c(ifile) = counts_c(ifile)+1
            enddo
            close(600+id)
          else
            write(*,*) "Didn't find file ", fname
          endif
        enddo
      !$OMP ENDDO

      !$OMP DO REDUCTION(+:s1)
        do ifile=0, ncpu
          s1 = s1 + counts_c(ifile)
        enddo
      !$OMP ENDDO

      if (s1==1) then
        !$OMP MASTER
          write(*,*) "No clump files found."
          skip_this = .true.
        !$OMP END MASTER
      endif
    !$OMP END PARALLEL

    if (skip_this) stop


    do i=0, ncpu-1
      counts_c(i+1) = counts_c(i+1)+counts_c(i)
    enddo



    !----------------------------------
    ! Allocate Arrays
    !----------------------------------

    nclumps = counts_c(ncpu)-1
    allocate(clmp_id(1:nclumps))
    allocate(parent_id(1:nclumps))


    !----------------------------------
    ! Read clump data
    !----------------------------------
    ! write(*,*) "Halo is:",halo

    !$OMP PARALLEL PRIVATE(id, ind, io, ifile, file_exists, &
    !$OMP fname, cpunr, xpeak, ypeak, zpeak, junk1, junk2, &
    !$OMP junk_char, clumpmass)
      id = omp_get_thread_num()
      !$OMP DO
        do ifile=1, ncpu
          call title(ifile, cpunr)
          fname=TRIM(srcdir//'/clump_'//outputnr//'.txt'//cpunr)
          open(unit=600+id, file=fname)
          read(600+id,*) ! skip header
          do ind=counts_c(ifile-1), counts_c(ifile)-1
            read(600+id,'(I8,x,I2,x,I10,x,I10,3(X,E18.9E2),3(X, A18),x,E18.9E2,x,A18)') &
                clmp_id(ind), junk1, parent_id(ind), junk2, xpeak, ypeak, zpeak, &
                junk_char, junk_char, junk_char, clumpmass, junk_char
            ! store halo position
            if (clmp_id(ind)==halo) then
              peakpos(1) = xpeak
              peakpos(2) = ypeak
              peakpos(3) = zpeak
              halomass = clumpmass
            endif
          enddo
          close(600+id)
        enddo
      !$OMP END DO

    !$OMP END PARALLEL


    read_clumps = .true.

  end subroutine read_clump_data



  subroutine read_galaxy_data(ignore_centrals_satellite_distinction)
    !------------------------------------------
    ! reads in galaxy data: associated clump,
    ! x, y, z, and associated particle.
    !------------------------------------------

    use omp_lib
    implicit none
    ! ignore distinction into centrals and satellites? Default: No.
    logical, intent(in), optional :: ignore_centrals_satellite_distinction

    logical :: ignoreCSD ! shorthand for above variable
    character(len=5)                      :: outputnr, cpunr
    character(len=80)                     :: fname
    integer                               :: i, ind, id, io, ifile
    integer                               :: igal, iorph, iclump
    integer                               :: cid
    integer(i8)                           :: pid
    real(dp)                              :: sm, xx, yy, zz
    logical                               :: file_exists, skip_this
    integer,  allocatable, dimension(:)   :: counts_g, counts_o
    integer(i8),allocatable,dimension(:)   :: idol
    integer,  allocatable, dimension(:)   :: gal_idl
    logical,  allocatable, dimension(:)   :: is_centrall
    real(dp), allocatable, dimension(:,:) :: xgl, xol
    real(dp), allocatable, dimension(:)   :: mgl, mol

    if (.not.present(ignore_centrals_satellite_distinction)) then
      ignoreCSD = .false.
    else
      ignoreCSD = ignore_centrals_satellite_distinction
    endif

    outputnr = srcdir(8:12)

    skip_this = .false.
    ngalaxies = 0
    norphans = 0

    if (.not.allocated(clmp_id).and..not.ignoreCSD) then
      if (read_clumps) write(*, *) "Clump data should be read...?"
      write(*, *) "Can't read in galaxy data without reading in clump data"
      stop
    endif


    !$OMP PARALLEL PRIVATE(i, id, io, ifile, file_exists, fname, cpunr, ind, &
    !$OMP                   cid, sm, xx, yy, zz, pid) &
    !$OMP DEFAULT(shared)
      id = omp_get_thread_num()

      !----------------------------------
      ! Read file sizes
      !----------------------------------

      !$OMP DO REDUCTION(+:ngalaxies, norphans)
        do ifile=1, ncpu
          call title(ifile, cpunr)
          fname=TRIM(srcdir//'/galaxies_'//outputnr//'.txt'//cpunr)
          inquire(file=fname, exist=file_exists)
          if (file_exists) then
            open(unit=600+id, file=fname)

            read(600+id, *) ! skip header
            do 
              read(600+id,'(I20,x,4(E20.12,x),I20)', iostat=io) &
                  cid, sm, xx, yy, zz, pid
              if (io/=0) exit
              if (cid > 0) then
                ngalaxies = ngalaxies + 1
              else
                norphans = norphans + 1
              endif
            enddo
            close(600+id)
          else
            write(*,*) "Didn't find file ", fname
          endif
        enddo
      !$OMP ENDDO

      if (ngalaxies+norphans==0) then
        !$OMP MASTER
          write(*,*) "No galaxy data found. Stopping."
          skip_this = .true.
        !$OMP END MASTER
      endif
    !$OMP END PARALLEL

    if (skip_this) stop



    !----------------------------------
    ! Allocate Arrays
    !----------------------------------

    allocate(xg(1:ngalaxies,1:3))
    allocate(gal_id(1:ngalaxies))
    allocate(mg(1:ngalaxies))
    if (.not.ignoreCSD) then
      allocate(is_central(1:ngalaxies))
    else
      write(*, *) "WARNING: ALLOCATING DUMMY IS_CENTRAL ARRAY"
      allocate(is_central(1:1))
    endif
    allocate(xo(1:norphans,1:3))
    allocate(ido(1:norphans))
    allocate(mo(1:norphans))
    allocate(counts_g(0:omp_get_max_threads()))
    allocate(counts_o(0:omp_get_max_threads()))
    counts_g = 0; counts_g(0)=1;
    counts_o = 0; counts_o(0)=1;
    ngalaxies_tot = ngalaxies + norphans

    xg = 0.
    gal_id = 0.
    mg = 0.
    if (.not.ignoreCSD) is_central = .false.
    xo = 0.
    ido = 0.
    mo = 0.


    !----------------------------------
    ! Read galaxy data
    !----------------------------------

    !$OMP PARALLEL PRIVATE(id, ind, io, ifile, file_exists, fname, cpunr, &
    !$OMP                  cid, sm, xx, yy, zz, pid, xgl, xol, mgl, mol, &
    !$OMP                  gal_idl, idol, igal, iorph, is_centrall, iclump)
      id = omp_get_thread_num()

      allocate(xgl(1:ngalaxies,1:3))
      allocate(gal_idl(1:ngalaxies))
      allocate(mgl(1:ngalaxies))
      if (.not.ignoreCSD) allocate(is_centrall(1:ngalaxies))
      allocate(xol(1:norphans,1:3))
      allocate(idol(1:norphans))
      allocate(mol(1:norphans))
      igal = 0
      iorph = 0

      !$OMP DO
        do ifile=1, ncpu
          call title(ifile, cpunr)
          fname=TRIM(srcdir//'/galaxies_'//outputnr//'.txt'//cpunr)
          inquire(file=fname, exist=file_exists)
          if (file_exists) then
            open(unit=600+id, file=fname)
            read(600+id,*) ! skip header

            ! loop over lines
            do 
              read(600+id,'(I20,x,4(E20.12,x),I20)', iostat=io) cid, sm, xx, yy, zz, pid
              if (io/=0) exit
              if (cid>0) then
                igal = igal + 1
                xgl(igal,1) = xx
                xgl(igal,2) = yy
                xgl(igal,3) = zz
                mgl(igal) = sm
                gal_idl(igal) = cid
                if (.not.ignoreCSD) then
                  do iclump=1, nclumps
                    if (cid == clmp_id(iclump)) then
                      if (parent_id(iclump) == cid) then
                        is_centrall(igal) = .true.
                      else
                        is_centrall(igal) = .false.
                      endif
                      exit
                    endif
                  enddo
                endif
              else
                iorph = iorph + 1
                xol(iorph,1) = xx
                xol(iorph,2) = yy
                xol(iorph,3) = zz
                mol(iorph) = sm
                idol(iorph) = pid
              endif
            enddo !loop over lines
            close(600+id)
          endif
        enddo 
      !$OMP END DO

      counts_o(id+1) = iorph
      counts_g(id+1) = igal

      !$OMP BARRIER
      !$OMP SINGLE
        do i=1, omp_get_max_threads()
          counts_o(i) = counts_o(i) + counts_o(i-1)
          counts_g(i) = counts_g(i) + counts_g(i-1)
        enddo
      !$OMP END SINGLE
      !$OMP BARRIER

      if (iorph>0) then
        xo(counts_o(id):counts_o(id+1)-1,1:3) = xol(1:iorph,1:3)
        ido(counts_o(id):counts_o(id+1)-1) = idol(1:iorph)
        mo(counts_o(id):counts_o(id+1)-1) = mol(1:iorph)
      endif
      if (igal>0) then
        xg(counts_g(id):counts_g(id+1)-1,1:3) = xgl(1:igal,1:3)
        gal_id(counts_g(id):counts_g(id+1)-1) = gal_idl(1:igal)
        mg(counts_g(id):counts_g(id+1)-1) = mgl(1:igal)
        if (.not.ignoreCSD) then
          is_central(counts_g(id):counts_g(id+1)-1) = is_centrall(1:igal)
        endif
      endif

      deallocate(xol,xgl,gal_idl,idol,mgl,mol)
    !$OMP END PARALLEL

    write(*,*) "Read in ", ngalaxies, "galaxies with clumps and ", norphans, "orphan galaxies."
    write(*,*) "Total galaxies in snapshot:", ngalaxies_tot

    read_galaxies = .true.

  end subroutine read_galaxy_data





  subroutine title(n,nchar)

    !----------------------------------------------------
    ! Get integer N in to zero padded string of length 5
    !----------------------------------------------------

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




  subroutine get_halo_string(halo_str, halo)

    !-----------------------------------------
    ! Writes the halo as a string.
    !-----------------------------------------
    
    implicit none
    integer, intent(in) :: halo
    character(len=20), intent(out) :: halo_str

    character(len=20) :: fmt_string


    fmt_string = '(I1)'
    if (halo/10          > 0) fmt_string = '(I2)'
    if (halo/100         > 0) fmt_string = '(I3)'
    if (halo/1000        > 0) fmt_string = '(I4)'
    if (halo/10000       > 0) fmt_string = '(I5)'
    if (halo/100000      > 0) fmt_string = '(I6)'
    if (halo/1000000     > 0) fmt_string = '(I7)'
    if (halo/10000000    > 0) fmt_string = '(I8)'
    if (halo/100000000   > 0) fmt_string = '(I9)'
    ! if (halo/1000000000  > 0) fmt_string = '(I10)'
    if (halo/1000000000  > 0) then
      write(*, *) "Halo ID is too big. Modify get_halo_string() subroutine"
      stop
    endif

    write(halo_str,fmt_string) halo

  end subroutine get_halo_string



  subroutine deallocate_io_module()
    !-------------------------------------------
    ! Deallocate arrays for proper cleanup
    !-------------------------------------------

    implicit none

    if (read_clumps) then
      deallocate(clmp_id, parent_id)
    endif
    if (read_galaxies) then
      deallocate(xg, xo, mg, mo, ido, gal_id, is_central)
    endif
  end subroutine deallocate_io_module


end module
