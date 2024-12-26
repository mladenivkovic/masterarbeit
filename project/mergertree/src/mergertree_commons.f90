module mergertree_commons

  use amr_parameters

  implicit none


  integer :: ncpu = 0
  integer :: nmost_bound = 0


  ! new for mergertree
  integer, allocatable, dimension(:)   :: halocount   !count of progenitor halos per proc
  integer, allocatable, dimension(:)   :: progenitor_list !contains list of progenitors
  real(dp), allocatable, dimension(:)   :: progenitor_mass !contains list of progenitors masses
  integer, allocatable, dimension(:,:) :: progenitor_tracers
  integer, allocatable, dimension(:,:) :: prog_desc_links ! counts particles of each descendant for each progenitor
  integer, allocatable, dimension(:)   :: prog_index    ! gets index for each progenitor for prog_desc_links
  integer, allocatable, dimension(:)   :: sorted_progs  ! list of sorted progenitors global peak id
  integer :: nprogs !number of progenitors after removing doubles

  integer, allocatable, dimension(:) :: main_prog, main_desc
  real(dp), allocatable, dimension(:) :: main_prog_merit, main_desc_merit


  !simulating ramses
  integer :: npeaks_max, npartmax 
  integer, allocatable, dimension(:)  :: npeaks, hfree 
  integer, allocatable, dimension(:)  :: clmpidp, idp, levelp
  real(dp), allocatable, dimension(:) :: clmp_mass_pb 
  integer, allocatable, dimension(:)  :: pind_by_id !"local" particle index by unique particle ID
  real(dp) :: particle_mass



end module mergertree_commons 
