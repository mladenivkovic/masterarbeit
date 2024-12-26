!=======================================================================
module io
!=======================================================================

  contains

!###############################################
!###############################################
!###############################################

!=======================================================================
    subroutine readargs( output_start, output_end )
!=======================================================================

    use mergertree_commons
    implicit none

    integer, intent (inout) :: output_start, output_end

    !===========================
    ! read command line args
    !===========================

    character(len = 5) :: arg
    integer :: nargs

    

    nargs = iargc()

    if (nargs == 0) then
      output_start = 1
      output_end = find_output_end(output_start)

    else if (nargs == 1) then

      call getarg(1, arg)
      read(arg, '(I5)' ) output_start
      output_end = find_output_end(output_start)
      
    else if (nargs == 2) then

      call getarg(1, arg)
      read(arg, '(I5)' ) output_start

      call getarg(2, arg)
      read(arg, '(I5)' ) output_end

      if (output_end < output_start) then
        stop("output end can't be smaller than start. Aborting")
      end if

    else
      stop("Too many arguments given. Can't handle that. Aborting.")
    end if


    contains

      integer function find_output_end(output_start)

        implicit none

        integer, intent(in) :: output_start

        integer :: output_end

        character (len=5)  :: dir
        character (len=80) :: dirname
        logical :: dir_exists = .true.


        output_end = output_start

        do while (dir_exists)
          output_end = output_end + 1

          call title(output_end+1, dir)
          dirname=TRIM('output_'//TRIM(dir)//'/')

          inquire(file=TRIM(dirname), exist=dir_exists)  

        end do

        if(debug) write(*,*) "DEBUG 0 Found output_start ", output_start, "output_end ", output_end

        find_output_end = output_end

      end function find_output_end
    end subroutine readargs 

!###############################################
!###############################################
!###############################################


!=======================================================================
    subroutine readdata(step)
!=======================================================================
    use mergertree_commons
    use amr_parameters
    implicit none 

    integer, intent (in) :: step

    character (len=5)  :: dir, idnr
    character (len=80) :: fileloc
    integer :: i=0, ipeak = 0, start_ind =0

    ! temporary
    !integer :: x

    if(verbose) write(*,*) "Reading in data for step", step

    !=========================
    ! read descendants data
    !=========================

    !read first for proc 1 which should always exist, 
    !then loop for all processors


    call title(step+1, dir)
    call title(1, idnr)
    fileloc=TRIM('output_'//TRIM(dir)//'/descendant_data.out'//TRIM(idnr))


    open(unit=666,file=fileloc,form='unformatted')
    read(666)  ncpu
    read(666)  npartmax
    read(666)  npeaks_max
    close(666)
  

    if (debug) write(*,'(A38,I4,2(A15,I8))') " DEBUG RD1 Read in initial info. Ncpu:", ncpu, "npartmax", npartmax, "npeaks_max", npeaks_max


    allocate(hfree(1:ncpu))
    allocate(npeaks(1:ncpu))

    allocate(clmpidp(1:ncpu*npartmax))
    allocate(levelp( 1:ncpu*npartmax))
    allocate(idp( 1:ncpu*npartmax))
    allocate(clmp_mass_pb(1:ncpu*npeaks_max))


    clmpidp = 0
    levelp = 0
    idp = 0
    clmp_mass_pb=0



    ! loop for all output

    start_ind = 0

    do i = 1, ncpu
      
      call title(i, idnr)
      fileloc=TRIM('output_'//TRIM(dir)//'/descendant_data.out'//TRIM(idnr))

      if (debug) write(*,*) "DEBUG RD2 Reading in from ", fileloc

      open(unit=666,file=fileloc,form='unformatted')

      read(666)  ncpu
      read(666)  npartmax
      read(666)  npeaks_max
      read(666)  nmost_bound 
      ! WARNING: assuming npartmax doesn't change during the course of a run
      read(666)  npeaks(i) 
      read(666)  hfree(i)
      read(666)  particle_mass
      read(666)  clmpidp((i-1)*npartmax+1 : i*npartmax) 
      read(666)  idp((i-1)*npartmax+1 : i*npartmax) 
      read(666)  levelp((i-1)*npartmax+1 : i*npartmax) 

      !introduce shifts for uniqueness

      read(666)  clmp_mass_pb(start_ind+1:start_ind+npeaks(i))
      start_ind = start_ind + npeaks(i)
      close(666)

    end do





    !=======================
    ! read progenitor data
    !=======================


    allocate(halocount(1:ncpu))

    do i = 1, ncpu
      call title(step, dir)
      call title(i, idnr)
      fileloc=TRIM('output_'//TRIM(dir)//'/progenitor_halocount.out'//TRIM(idnr))
      if (debug) write(*,*) "DEBUG RD3 Reading in from ", fileloc

      open(unit=666,file=fileloc,form='formatted')
      read(666, '(I10)') halocount(i)
      close(666)

    end do 



    if (debug) then 
      write(*,*) "DEBUG RD4 data read in:"
      write(*,*) "npartmax", npartmax, "npeaks_max", npeaks_max, "ncpu", ncpu 
      write(*,*) "halocount:", halocount
    end if



    allocate(progenitor_list(1:sum(halocount)))
    allocate(progenitor_mass(1:sum(halocount)))
    allocate(progenitor_tracers(1:sum(halocount), 1:nmost_bound))

    progenitor_list = 0
    progenitor_tracers = 0


    start_ind = 0
    do i = 1, ncpu

      call title(i, idnr)
      fileloc=TRIM('output_'//TRIM(dir)//'/progenitor_data.out'//TRIM(idnr))
      if (debug) write(*,*) "DEBUG RD5 Reading in from ", fileloc

      open(unit=666,file=fileloc,form='formatted')
      do ipeak = 1, halocount(i)
        !write(*,*) "Reading progenitors cpu", i, "ipeak", ipeak
        read(666,'(I10,E15.6,250I10)') progenitor_list(ipeak+start_ind), progenitor_mass(ipeak+start_ind), progenitor_tracers(ipeak+start_ind,:)
      end do 

      start_ind = start_ind + halocount(i)

      close(666)

    end do



    !temp for testing
    !do ipeak = 1, halocount_tot 
      !write(*,*) "Peak", progenitor_list(ipeak) 
      !do i = 0, 24
      !  write(*,*) (progenitor_tracers(ipeak, x), x = 10*i+1, 10*(i+1) )
      !end do
      !write(*,*) "==============================="
      !
    !end do


    end subroutine readdata 

!###############################################
!###############################################
!###############################################




!============================================
    subroutine write_trees(step)
!============================================

    use mergertree_commons

    implicit none

    integer, intent(in) :: step

    character (len=5)  :: dir, idnr
    character (len=80) :: fileloc
    integer:: ipeak, iprog

    logical, dimension(1:nprogs) :: printed

    printed = .false.

    call title(step+1, dir)
    call title(1, idnr) !TODO for now
    fileloc=TRIM('output_'//TRIM(dir)//'/mergertree.txt'//TRIM(idnr))

    if(debug) write(*,*) "DEBUG-OUTPUT Writing output to file ", fileloc

    open(unit=666,file=fileloc,form='formatted')
    write(666, '(2(A15))') "clump", "progenitor"
    do ipeak = 1, ncpu*npeaks_max
      if (clmp_mass_pb(ipeak) > 0) then
        write(666,'(2(I15))') ipeak, main_prog(ipeak)
        if (main_prog(ipeak)>0) printed(prog_index(main_prog(ipeak))) = .true.
      end if
    end do

    do iprog = 1, nprogs
      if (.not. printed(iprog) ) then
        write(666, '(2(I15))') main_desc(iprog), sorted_progs(iprog)
      end if
    end do



    close(666)


end subroutine write_trees

!###############################################
!###############################################
!###############################################


end module io

