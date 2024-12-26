!========================================================
! Count how many orphans are duplicate - have same partice IDs
!========================================================


!===================================
program count_cuplicate_orphans
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

  integer , dimension(:), allocatable   :: order    ! sorted orphan order by ID
  integer :: i, doubles
  real(dp) :: mass_av, mass_max

  !---------------------------
  ! Set-up
  !---------------------------

  call get_srcdir_from_cmdlineargs() ! from io_module
  call read_info() ! from io_module
  call read_galaxy_data(ignore_centrals_satellite_distinction=.true.)

  write(*, *) "Working for ", srcdir

  ! sort orphans by their ID
  allocate(order(1:norphans))
  do i = 1, norphans
    order(i) = i
  enddo
  call quick_sort_int(ido, order, norphans)

  doubles = 0
  mass_av = 0.d0
  mass_max = 0.d0

  do i = 1, norphans - 1
    if (ido(i) == ido(i+1)) then
      write(*, "(3I12)") i, i+1, ido(i)
      doubles = doubles + 1
      mass_av = mass_av + mo(order(i))
      mass_max = max(mass_max, mo(order(i)))
    endif
  enddo
  write(*, '(A,x,I12,x,A,x,E12.3,x,A,x,E12.3)') "Total duplicates:", doubles, &
          "Mass max", mass_max, "mass av", mass_av/dble(doubles)

  call deallocate_io_module()
  deallocate(order)


end program
