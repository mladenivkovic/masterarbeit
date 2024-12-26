! Module containg constants and parameters

module constants_and_parameters

  !------------------
  ! Parameters
  !------------------

  integer, parameter  :: dp = kind(1.0d0) 
  integer, parameter  :: i8 = selected_int_kind(9)


  !------------------
  ! Constants
  !------------------

  real(dp), parameter :: pi        = 3.14159265359d+0
  real(dp), parameter :: twopicube = (2.d0*pi)**3
  real(dp), parameter :: euler     = 2.7182818284590
  real(dp), parameter :: log10_2   = 0.301029995663981 ! log_10(2)



  real(dp), parameter :: Mpc       = 3.08567758d+24         ! Mpc in cm
  real(dp), parameter :: M_Sol     = 1.99847d+33            ! solar mass in g
  real(dp), parameter :: G         = 4.30091d-9             ! Mpc MSol-1 km^2/s^2


end module
