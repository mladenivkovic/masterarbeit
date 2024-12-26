module SMHM_relations
  
  use constants_and_parameters
  use io_module, only: mg, mo, ngalaxies, norphans, aexp
  implicit none

  logical :: SMHM_seed_set = .false.

contains

  subroutine calculate_sm(model)
    !------------------------------------------
    ! Subroutine to compute the stellar mass
    ! assuming the read-in 'stellar_mass' array
    ! doesn't actually contain galaxy masses,
    ! but halo masses to be used for the SMHM 
    ! relation.
    !
    ! model: which Stellar mass - halo - mass
    ! relation to use
    !------------------------------------------
    
    implicit none
    character(len=4) :: model
    integer :: i

    if (model=='behr') then
      !$OMP PARALLEL 
      !$OMP DO PRIVATE(i)
        do i=1, ngalaxies
          mg(i) = smhm_behroozi(mg(i), aexp)
        enddo
      !$OMP END DO NOWAIT
      !$OMP DO PRIVATE(i)
        do i=1, norphans
          mo(i) = smhm_behroozi(mo(i), aexp)
        enddo
      !$OMP END DO
      !$OMP END PARALLEL

    else if (model=='most') then
      !$OMP PARALLEL 
      !$OMP DO PRIVATE(i)
        do i=1, ngalaxies
          mg(i) = smhm_behroozi(mg(i), aexp)
        enddo
      !$OMP END DO NOWAIT
      !$OMP DO PRIVATE(i)
        do i=1, norphans
          mo(i) = smhm_behroozi(mo(i), aexp)
        enddo
      !$OMP END DO
      !$OMP END PARALLEL
    endif
  end subroutine calculate_sm


  real(dp) function smhm_behroozi(m,a)
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
    use constants_and_parameters

    implicit none
    real(dp), intent(in) :: m,a ! mass

    real(dp), parameter :: M_10    =   11.514 
    real(dp), parameter :: M_1a    = -  1.793
    real(dp), parameter :: M_1z    = -  0.251
    real(dp), parameter :: e_0     = -  1.777
    real(dp), parameter :: e_a     = -  0.006
    real(dp), parameter :: e_z     =    0.000
    real(dp), parameter :: e_a2    = -  0.119
    real(dp), parameter :: alpha_0 = -  1.412
    real(dp), parameter :: alpha_a =    0.731
    real(dp), parameter :: delta_0 =    3.508 
    real(dp), parameter :: delta_a =    2.608
    real(dp), parameter :: delta_z = -  0.043
    real(dp), parameter :: gamma_0 =    0.316 
    real(dp), parameter :: gamma_a =    1.319 
    real(dp), parameter :: gamma_z =    0.279
    real(dp), parameter :: xi_0    =    0.218
    real(dp), parameter :: xi_a    = -  0.023

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
    if (.not.SMHM_seed_set) then
      call random_seed(size=n)
      allocate(seed(n))
      call system_clock(count=clock)
      seed = clock + omp_get_thread_num() + 53*[(n-i, i=1, n)]
      call random_seed(put=seed)
      deallocate(seed)
      SMHM_seed_set = .true.
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
    smhm_behroozi = 10.d0**(loge+logM1 + behr_f(log10(m)-logM1, alpha, delta, gam)- f0 + xi)
  end function smhm_behroozi


  real(dp) function behr_f(x, alpha, delta, gam)
    implicit none
    real(dp), intent(in) :: x, alpha, delta, gam
    behr_f = -log10(10.d0**(alpha*x) + 1) + delta*log10(1.d0+exp(x))**gam/(1.d0 + exp(10.d0**(-x)))
  end function behr_f


  real(dp) function smhm_moster(m,a)
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

    real(dp), parameter :: M_10     =  11.590
    real(dp), parameter :: M_11     =   1.195
    real(dp), parameter :: N_10     =   0.0351
    real(dp), parameter :: N_11     = - 0.0247
    real(dp), parameter :: beta_10  =   1.376
    real(dp), parameter :: beta_11  = - 0.826
    real(dp), parameter :: gamma_10 =   0.608
    real(dp), parameter :: gamma_11 =   0.329

    real(dp) :: logM1, M1
    real(dp) :: N   
    real(dp) :: beta
    real(dp) :: gam 
    real(dp) :: am1

    am1 = a - 1.d0

    logM1 = M_10     + M_11*am1
    M1    = 10**logM1
    N     = N_10     + N_11*am1
    beta  = beta_10  + beta_11*am1
    gam   = gamma_10 + gamma_11*am1

    smhm_moster = m*2.d0*N/((m/M1)**(-beta)+(m/M1)**gam)
  end function smhm_moster


end module
