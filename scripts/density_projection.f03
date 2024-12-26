! Module for creating density projections
!--------------------------------------------



module density_projection

  use constants_and_parameters
  use io_module
  use omp_lib
  implicit none

  character(len=3) :: interpolation = 'ngp'    ! which interpolation to use: cic for cloud in cell, ngp for nearest grid point
  integer :: nc ! number of cells for density fields

  real(dp), dimension(:,:),allocatable:: density_field_xy ! Density image on xy plane
  real(dp), dimension(:,:),allocatable:: density_field_yz ! Density image on yz plane
  real(dp), dimension(:,:),allocatable:: density_field_xz ! Density image on xz plane

  real(dp) :: xmax, xmin, ymax, ymin, zmax, zmin ! min/max values of coordinates
  real(dp) :: image_width    ! width of square that will be plotted
  real(dp) :: dc    ! cell width
  real(dp) :: vol   ! cell "volume" (surface actually)

contains


  subroutine get_density_projection(x, m, n)

    !-----------------------------------------------------------------
    ! Get density field projection along all three fundamental planes.
    !
    ! x: array to use for positions. ALL of the values in the array 
    !     will be plotted. If you want a subset of the array, only 
    !     pass the subset!
    ! y: array to use for masses
    ! n: array size
    ! interpolation: 3 letter string.
    !                 'ngp' for nearest grid point
    !                 'cic' for cloud in cell
    ! nc: how many cells in each direction to use
    !-----------------------------------------------------------------

    implicit none

    real(dp), dimension(1:n, 1:3), intent(in) :: x
    real(dp), dimension(1:n), intent(in)      :: m
    integer, intent(in)                       :: n

    real(dp) :: xc, dx, yc, dy, zc, dz

    ! find center and width of field
    !-----------------------------------

    xmax = maxval(x(:, 1))
    xmin = minval(x(:, 1))
    ymax = maxval(x(:, 2))
    ymin = minval(x(:, 2))
    zmax = maxval(x(:, 3))
    zmin = minval(x(:, 3))

    dx = xmax - xmin
    xc = 0.5*(xmax+xmin)
    dy = ymax - ymin
    yc = 0.5*(ymax+ymin)
    dz = zmax - zmin
    zc = 0.5*(zmax+zmin)

    image_width = max(dx, dy, dz)*1.02 ! add a little border
    xmax = xc+image_width/2
    xmin = xc-image_width/2
    ymax = yc+image_width/2
    ymin = yc-image_width/2
    zmax = zc+image_width/2
    zmin = zc-image_width/2

    dc = image_width / nc
    vol = dc*dc



    !---------------------------------------------------------
    write(*,*) "Computing density fields"
    !---------------------------------------------------------

    if (interpolation=='ngp') then
      call get_ngp_density_field(x, m, n)
    else if (interpolation=='cic') then
      call get_cic_density_field(x, m, n)
    else
      write(*,*) "Didn't recognize interpolation method ", interpolation
      stop
    endif

  end subroutine get_density_projection

  
  subroutine get_ngp_density_field(x, m, n)
    !------------------------------------------------------------
    ! Computes the density field and overdensity field delta
    !------------------------------------------------------------

    implicit none

    real(dp), dimension(1:n, 1:3), intent(in) :: x
    real(dp), dimension(1:n), intent(in)      :: m
    integer, intent(in)                       :: n

    integer  :: i, j, k, p
    real(dp) :: f_image 

    f_image = dble(nc) / image_width

    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP PRIVATE(i,j,k,p)
      !$OMP DO
        do p=1, n
          i = int((x(p,1)-xmin)*f_image)+1
          j = int((x(p,2)-ymin)*f_image)+1
          k = int((x(p,3)-zmin)*f_image)+1
          !$OMP CRITICAL
            density_field_xy(i,j) = density_field_xy(i,j)+m(p)
            density_field_yz(j,k) = density_field_yz(j,k)+m(p)
            density_field_xz(i,k) = density_field_xz(i,k)+m(p)
          !$OMP END CRITICAL
        enddo
      !$OMP END DO
    !$OMP END PARALLEL

    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP PRIVATE(i,j,k,p)
      !$OMP DO COLLAPSE(2)
        do i=1, nc
          do j=1, nc
            density_field_xy(i,j) = density_field_xy(i,j)/vol
            density_field_yz(i,j) = density_field_yz(i,j)/vol
            density_field_xz(i,j) = density_field_xz(i,j)/vol
          enddo
        enddo
      !$OMP ENDDO
    !$OMP END PARALLEL

  end subroutine get_ngp_density_field


  subroutine get_cic_density_field(x, m, n)
    !------------------------------------------------------------
    ! Computes the density field and overdensity field delta
    !------------------------------------------------------------

    implicit none

    real(dp), dimension(1:n, 1:3), intent(in) :: x
    real(dp), dimension(1:n), intent(in)      :: m
    integer, intent(in)                       :: n

    integer  :: p, iup, idown, jup, jdown, kup, kdown
    real(dp) :: rho, xup, yup, zup, hdc, cv

    hdc = dc/2
    cv = dc**2

    !$OMP PARALLEL DEFAULT(shared) &
    !$OMP PRIVATE(p,iup,idown,jup,jdown,kup,kdown,rho,xup,yup,zup)
      !$OMP DO
        do p=1, n
          rho = m(p)/cv
          iup   = int((x(p,1)-xmin+hdc)/image_width*nc)+1
          if (iup>nc) write(*,*) "IUP > nc", x(p,1), xmin, hdc, (x(p,1)-xmin+hdc)/image_width 
          if (iup==1) write(*,*) "IUP = 1", x(p,1), xmin, hdc, (x(p,1)-xmin+hdc)/image_width 
          idown = iup-1
          jup   = int((x(p,2)-ymin+hdc)/image_width*nc)+1
          if (jup>nc) write(*,*) "JUP > nc", x(p,2), ymin, hdc, (x(p,2)-ymin+hdc)/image_width 
          if (jup==1) write(*,*) "JUP = 1", x(p,2), ymin, hdc, (x(p,2)-ymin+hdc)/image_width 
          jdown = jup-1
          kup   = int((x(p,3)-zmin+hdc)/image_width*nc)+1
          if (kup>nc) write(*,*) "KUP > nc", x(p,3), zmin, hdc, (x(p,3)-zmin+hdc)/image_width 
          if (kup==1) write(*,*) "KUP = 1", x(p,3), zmin, hdc, (x(p,3)-zmin+hdc)/image_width 
          kdown = kup-1

          xup = x(p,1) -xmin + hdc - (iup-1)*dc ! fraction in x-direction of surface/volume in upper cell
          yup = x(p,2) -ymin + hdc - (jup-1)*dc ! fraction in y-direction of surface/volume in upper cell
          zup = x(p,3) -zmin + hdc - (kup-1)*dc ! fraction in z-direction of surface/volume in upper cell

          ! No periodicity, present, so not needed
          ! if (iup>nc)  iup   = iup-nc
          ! if (idown<1) idown = nc+idown
          ! if (jup>nc)  jup   = jup-nc
          ! if (jdown<1) jdown = nc+jdown

          !$OMP CRITICAL
            density_field_xy(iup,   jup  ) = density_field_xy(iup,   jup  )+ xup      * yup      * rho
            density_field_xy(idown, jup  ) = density_field_xy(idown, jup  )+ (dc-xup) * yup      * rho
            density_field_xy(iup,   jdown) = density_field_xy(iup,   jdown)+ xup      * (dc-yup) * rho
            density_field_xy(idown, jdown) = density_field_xy(idown, jdown)+ (dc-xup) * (dc-yup) * rho

            density_field_yz(jup,   kup  ) = density_field_yz(jup,   kup  )+ yup      * zup      * rho
            density_field_yz(jdown, kup  ) = density_field_yz(jdown, kup  )+ (dc-yup) * zup      * rho
            density_field_yz(jup,   kdown) = density_field_yz(jup,   kdown)+ yup      * (dc-zup) * rho
            density_field_yz(jdown, kdown) = density_field_yz(jdown, kdown)+ (dc-yup) * (dc-zup) * rho

            density_field_xz(iup,   kup  ) = density_field_xz(iup,   kup  )+ xup      * zup      * rho
            density_field_xz(idown, kup  ) = density_field_xz(idown, kup  )+ (dc-xup) * zup      * rho
            density_field_xz(iup,   kdown) = density_field_xz(iup,   kdown)+ xup      * (dc-zup) * rho
            density_field_xz(idown, kdown) = density_field_xz(idown, kdown)+ (dc-xup) * (dc-zup) * rho
          !$OMP END CRITICAL
        enddo
      !$OMP END DO


      !$OMP DO COLLAPSE(2)
        do iup=1, nc
          do jup=1, nc
            density_field_xy(iup,jup) = density_field_xy(iup,jup)/vol
            density_field_yz(iup,jup) = density_field_yz(iup,jup)/vol
            density_field_xz(iup,jup) = density_field_xz(iup,jup)/vol
          enddo
        enddo
      !$OMP ENDDO
    !$OMP END PARALLEL

  end subroutine get_cic_density_field




  subroutine write_density_projections(fname)
    !--------------------------------------------------------
    ! write image of density fields to file that can be
    ! plotted directly with plot_fortran_halo.py
    ! fname: file name to use
    !--------------------------------------------------------

    character(len=80) :: fname

    if (.not.infofile_read) then
      write(*, *) "Missing data from info file to dump density projections :/"
      stop
    endif

    open(unit=666, file=fname, form='unformatted')
    write(666) nc
    write(666) aexp
    write(666) unit_l
    write(666) xmax
    write(666) xmin
    write(666) ymax
    write(666) ymin
    write(666) zmax
    write(666) zmin
    write(666) density_field_xy
    write(666) density_field_yz
    write(666) density_field_xz
    close(666)

    write(*,'(A14,A)') " Written file ", TRIM(fname)

  end subroutine write_density_projections 




  subroutine allocate_density_fields()
    implicit none

    ! Allocate arrays
    !-----------------------

    allocate(density_field_xy(1:nc, 1:nc))
    density_field_xy = 0
    allocate(density_field_yz(1:nc, 1:nc))
    density_field_yz = 0
    allocate(density_field_xz(1:nc, 1:nc))
    density_field_xz = 0
  end subroutine allocate_density_fields




  subroutine deallocate_density_fields()

    ! Clean up after yourself
    !-----------------------------

    deallocate(density_field_xy, density_field_xz, density_field_yz)

  end subroutine deallocate_density_fields



end module
