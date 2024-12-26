module amr_parameters

  ! Define real types
  integer,parameter::sp=kind(1.0E0)
  integer,parameter::dp=kind(1.0D0) ! default
  integer,parameter::qdp=kind(1.0_8) ! real*8
  integer,parameter::MAXOUT=1000
  integer,parameter::MAXLEVEL=100

  ! Define integer types (for particle IDs mostly)
  ! Warning: compiler needs to accept fortran standard 2003.
  ! Specific numbers for fortran kinds are, in principle, implementation
  ! dependent, so "i8b=8" with "integer(i8b)" is implementation-dependent.
  ! See portability note for non-gcc compilers: (boud 2016-11-29) -
  ! https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
  ! The selected_int_kind approach below is more likely to be portable:
  integer,parameter::i4b=selected_int_kind(9) ! since log(2*10^9)/log(2)=30.9
  !integer,parameter::i8b=4  ! default long int are 4-byte int
  integer,parameter::i8b=selected_int_kind(9) ! since log(2*10^9)/log(2)=30.9

  ! Number of dimensions
  integer,parameter::ndim=3
  integer,parameter::twotondim=2**ndim
  integer,parameter::threetondim=3**ndim
  integer,parameter::twondim=2*ndim

  ! Vectorization parameter
  integer,parameter::nvector=500  ! Size of vector sweeps

  integer, parameter :: nstride = 65536

  ! Run control
  logical::verbose =.true.   ! Write everything
  logical::hydro   =.false.   ! Hydro activated
  logical::pic     =.false.   ! Particle In Cell activated
  logical::poisson =.false.   ! Poisson solver activated
  logical::cosmo   =.false.   ! Cosmology activated
  logical::star    =.false.   ! Star formation activated
  logical::sink    =.false.   ! Sink particles activated
  logical::rt      =.false.   ! Radiative transfer activated
  logical::debug   =.false.   ! Debug mode activated
  logical::writeinfo   =.true.   ! Debug mode activated
  logical::static  =.false.   ! Static mode activated
  logical::static_dm=.false.  ! Static mode for dm only activated
  logical::static_gas=.false. ! Static mode for gas only activated
  logical::static_stars=.false.! Static mode for stars only activated
  logical::tracer  =.false.   ! Tracer particles activated
  logical::lightcone=.false.  ! Enable lightcone generation
  logical::clumpfind=.false.  ! Enable clump finder
  logical::aton=.false.       ! Enable ATON coarse grid radiation transfer


end module amr_parameters
