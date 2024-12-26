!=================================================
! Read maximal values used for npart/ngrids
! Usage:
! read_maxarrayvals.o <outputnr>
!=================================================

program findmax

  implicit none

  integer           :: ncpu
  character(len=12) :: srcdir


  call read_cmdlineargs()
  call read_info()
  call find_maxval('amr_ ', 6)
  call find_maxval('part_', 2)






contains


  !========================================
  subroutine find_maxval(prefix, skip)
  !========================================
    implicit none
    character(len=5), intent(in) :: prefix
    integer, intent(in)          :: skip

    character(len=80):: fname
    character(len=5) :: dirnr, cpunr
    integer :: ifile, current, val, i

    current = 0
    do ifile=1, ncpu
      dirnr = srcdir(8:12)
      call title(ifile, cpunr)
      fname=TRIM(srcdir//'/'//TRIM(prefix)//dirnr//'.out'//cpunr)
      open(unit=666,form='unformatted', file=fname)
      do i=1, skip
        read(666)
      enddo
      read(666) val
      if (current < val) current = val
      close(666)
    enddo

    write(*,*) "Max for ", prefix, " : ", current 

  end subroutine find_maxval




  !========================================
  subroutine read_cmdlineargs()
  !========================================
    !-----------------------
    ! reads cmd line args
    !-----------------------
    character(len=12) :: arg
    character(len=5)  :: dirnr_str
    integer           :: dirnr, i
    logical           :: exists

    do i = 1, iargc()
      call getarg(i, arg)
      read(arg, *) dirnr
      call title(dirnr, dirnr_str)
      srcdir = TRIM('output_'//dirnr_str)
      write(*,*) "Working for srcdir ", srcdir
    end do
  end subroutine read_cmdlineargs





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

    outnr_str = srcdir(8:12)
    fname=TRIM(srcdir//'/info_'//outnr_str//'.txt')
    
    inquire(file=fname, exist=file_exists)
    if (file_exists) then

      open(unit=666, file=fname)
      read(666, '(A13,I11)') junk, ncpu
      close(666)

    else
      write(*,*) "Didn't find file", fname
      stop
    endif
  end subroutine read_info





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





end program findmax
