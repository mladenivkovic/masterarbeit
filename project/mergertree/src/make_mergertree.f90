!=============================================
program make_mergertree
!=============================================

  use mergertree_commons
  use io
  implicit none

  integer :: output_start=0
  integer :: output_end=0


  integer :: step = 0

  ! loopers
  !integer :: i



  !==================
  ! IMPORTANT NOTES
  !==================
  
  !main progenitor = -4: descendant doesn't have any progenitor,
  !but is a new halo.
  !
  !main progenitor = -5: initial value 
  !
  !main progenitor = -6: descendant doesn't have a progenitor, because its progenitor is already taken
  !
  !main descendant = -3 : initial value
  !
  !
  !main descendant = -7: progenitor doesn't have any descendant, it disappears. 
  !






  !====================
  ! read in data
  !====================

  write(*,*) "=================================="
  write(*,*) "Started making merger trees"
  write(*,*) "=================================="



  call readargs(output_start, output_end)


  !====================
  ! Construction loop
  !====================
  do step = output_start, output_end-1

    write(*,*) 
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
    write(*,'(A24,I5)') " Started for output step ", step
    write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"

    !read in data for this step
    call readdata(step)

    if(debug) write(*,*) "main found halocount", sum(halocount)
    
    call prepare_arrays()
    call create_links() 
    call write_trees(step)


    !finish up step
    call deallocate_mergertree_loop()

    write(*,*) "============== Finished step", step, "==============="

  end do

  
  write(*,*) "================= Finished program. ===================="






end program  make_mergertree
!################################################################
!################################################################
!################################################################
!=============================================
subroutine prepare_arrays()
!=============================================
  
  use mergertree_commons

  implicit none

  integer, allocatable, dimension (:) :: ind, ind2, temp, idp_copy
  integer :: i, j
  integer :: n, h, maxind


  ! set up
  n = maxval(progenitor_list)
  h = sum(halocount)
  maxind = maxval(idp)
  allocate(prog_index(1:n))
  allocate(ind(1:h))
  allocate(ind2(1:ncpu*npartmax))
  allocate(temp(1:h))
  allocate(pind_by_id(1:maxind))
  allocate(idp_copy(1:ncpu*npartmax))


  prog_index = 0






  !========================================================
  ! create hash table to get particle index by its id 
  !========================================================

  idp_copy(:) = idp(:)
  call quick_sort_int(idp_copy, ind2, ncpu*npartmax)
  
 
  pind_by_id(:) = ind2(ncpu*npartmax-maxind+1: ncpu*npartmax)


  do i =  ncpu*npartmax-maxind+1, ncpu*npartmax
    if (ind2(i)<1) then
      write(*,*) "Error 1"
    end if
  end do

  deallocate(ind2, idp_copy)


  !========================================================
  ! create hash table to get progenitor position by its ID 
  !========================================================

  temp(:) = progenitor_list(:)


  call quick_sort_int(temp, ind, h )


  !remove doubles
  i = 1
  do while (i < h) 
    if (temp(i+1) == temp(i)) then
      do j = i, h-1
        if (temp(j+1) == temp(i)) then
          temp(j+1) = HUGE(8)
        else
          i = j
          exit
        end if
      end do
    end if
    i = i + 1
  end do



  call quick_sort_int(temp, ind, h )
  

  !assume maximal value for nprogs
  nprogs = h

  do i = 1, h
    if (temp(i) < HUGE(8)) then 
      prog_index(temp(i)) = i
    else
      nprogs = i-1
      exit
    end if
  end do


  allocate (sorted_progs(1:nprogs))
  sorted_progs(1:nprogs) = temp(1:nprogs)

end subroutine prepare_arrays 
!################################################################
!################################################################
!################################################################
!=============================================
subroutine create_links()
!=============================================

  use mergertree_commons
 
  implicit none
  integer, allocatable, dimension (:) :: tracercount
  integer :: ipeak, iprog, ipart, thisprog, thisdesc, part_id
  integer :: h
  real(dp) :: merit_max, merit_calc, npart_desc, npart_prog
  integer :: merit_max_id, store_id = 0
  logical :: is_first


  integer :: counter, c


  ! set up
  h = sum(halocount)

  allocate(tracercount(1:ncpu*npeaks_max))
  tracercount = 0

  allocate(prog_desc_links(1:nprogs, 1:ncpu*npeaks_max))
  prog_desc_links = 0



  !========================================
  ! Count particles of each prog in desc
  !========================================

  counter = 0

  !TODO: not nprogs?
  do iprog = 1, h ! loop over progenitor particle list

    if (progenitor_list(iprog) > 0) then 

      thisprog = prog_index(progenitor_list(iprog))
      c = 0

      do ipart = 1, nmost_bound
        
        part_id = progenitor_tracers(iprog, ipart)
        
        if (part_id > 0 ) then
          thisdesc = abs(clmpidp(pind_by_id(part_id)))

          ! add particle
          if (thisdesc /= 0) then !tracer particle might not be in a halo anymore

            prog_desc_links(thisprog, thisdesc) = prog_desc_links(thisprog, thisdesc) + 1
          ! add to total sum for descendant
            
            tracercount(thisdesc) = tracercount(thisdesc) + 1
            
            counter = counter + 1
            c = c + 1
          end if
        end if
      end do

      if(debug) write(*,*) "DEBUG 2-1 counted particles for progenitor", progenitor_list(iprog), "iprog", iprog, "count", c
    end if
  end do 
  
  if(debug) write(*,*) "DEBUG 2-2 counter after linking:", counter



!  if (debug) then
  !  write(*,*) "DEBUG 2-2-1"
  !  
  !  do iprog = 1, h
  !
  !    write(*,*) "PROGENITOR", progenitor_list(iprog), "iprog", iprog
  !    do pi = 1, 12
  !      write(*,'(20(I7))') (progenitor_tracers(iprog, ipart), ipart = 20*(pi-1)+1, 20*pi)
  !    end do
  !
  !    write(*,'(20(I7))') (progenitor_tracers(iprog, ipart), ipart = 240, 250)
  !    write(*,*)
  !  end do
  !
  !  write(*, *) "=== CLUMPIDP", clmpidp(pind_by_id(10878))
  !  write(*,*)
  !  write(*,*)
  !
  !end if
!


  if (debug) then
    do ipeak = 1, ncpu*npeaks_max
      if (clmp_mass_pb(ipeak) > 0) write(*,*) "DEBUG 2-4-1 Found descendant mass", ipeak, clmp_mass_pb(ipeak)
    end do

    do iprog = 1, nprogs
      write(*,*) "DEBUG 2-4-2 Found progenitor mass", progenitor_list(iprog), progenitor_mass(iprog)
    end do
  end if




  !temp
  do ipeak = 1, ncpu*npeaks_max
    if (tracercount(ipeak) > 0 ) then
      if (debug) write(*,*) "DEBUG 2-3 ipeak", ipeak, "tracercount", tracercount(ipeak)
      !do iprog = 1, h
        !if (progenitor_list(iprog) < HUGE(8)) then
        !  thisprog = prog_index(progenitor_list(iprog))
        !  if (thisprog > 0) then 
        !    write(*,*) "Progenitor:", progenitor_list(iprog), "tracer particles found:", prog_desc_links(thisprog, ipeak)
        !  end if
        !else
        !  exit
        !end if
      !end do
    
    end if
    
  end do






  !===============================================
  ! Determine main descendant for each progenitor
  !===============================================

  allocate(main_desc(1:nprogs))
  main_desc = -3
  allocate(main_prog(1:ncpu*npeaks_max))
  main_prog = -5
  allocate(main_prog_merit(1:ncpu*npeaks_max))
  main_prog_merit = -10
  allocate(main_desc_merit(1:nprogs))
  main_desc_merit = -11




  do iprog = 1, nprogs

    merit_max = 0
    merit_max_id = 0 

    do ipeak = 1, ncpu*npeaks_max 
     
      if (prog_desc_links(iprog, ipeak) > 0) then
      
        if (debug) write(*,*) "DEBUG 2-5-1 mass", clmp_mass_pb(ipeak), "pm", particle_mass, "particle links", prog_desc_links(iprog, ipeak)
        if (debug) write(*,*) "DEBUG 2-5-2 iprog", iprog, "ipeak", ipeak 

        npart_desc =  clmp_mass_pb(ipeak) / particle_mass
        merit_calc =  DBLE(prog_desc_links(iprog, ipeak)) / npart_desc

        if (merit_calc > merit_max) then
          merit_max = merit_calc
          merit_max_id = ipeak
        end if
        
      end if
    end do

    ! TODO: if necessary, insert second condition here
    if (merit_max_id > 0) then 
      main_desc(iprog) = merit_max_id
      main_desc_merit(iprog) = merit_calc
    else
      main_desc(iprog) = -7
    end if


  end do




  !===============================================
  ! Determine main progenitor for each descendant
  !===============================================


  do ipeak = 1, ncpu*npeaks_max

    if (clmp_mass_pb(ipeak) > 0) then

      merit_max = 0
      merit_max_id = 0 
      if (debug) write(*,*) "DEBUG 2-6 checking ipeak", ipeak

      do iprog = 1, nprogs
        if (prog_desc_links(iprog, ipeak) > 0) then
         
          npart_prog =  clmp_mass_pb(ipeak) / particle_mass
          merit_calc =  DBLE(prog_desc_links(iprog, ipeak)) / npart_prog

          if (merit_calc > merit_max) then
            merit_max = merit_calc
            merit_max_id = sorted_progs(iprog)
          end if
         
        end if
      end do
      

      ! TODO: if necessary, insert second condition here
      if (merit_max_id > 0) then 
        main_prog(ipeak) = merit_max_id
        main_prog_merit(ipeak) = merit_calc
      else
        main_prog(ipeak) = -4
      end if
    end if


  end do







  !===============================================
  ! Create tree:
  ! Check progenitors for merging, identify main
  ! progenitor
  !===============================================


  is_first = .true.
  iprog = 1

  do while (iprog <= nprogs)
   
    if (is_first) store_id = main_desc(iprog)
    is_first = .true.

    ! if this progenitor isnt' the main progenitor of 
    ! its own main descendant
    if (sorted_progs(iprog) /= main_prog(main_desc(iprog))) then

      prog_desc_links(iprog, main_desc(iprog)) = -1
      merit_max = 0
      merit_max_id = 0 
      
      !try finding another
      do ipeak = 1, ncpu*npeaks_max
       
        if (prog_desc_links(iprog, ipeak) > 0) then
        
          if (debug) write(*,*) "DEBUG 2-11-1 mass", clmp_mass_pb(ipeak), "pm", particle_mass, "particle links", prog_desc_links(iprog, ipeak)
          if (debug) write(*,*) "DEBUG 2-11-2 iprog", iprog, "ipeak", ipeak 
      
          !npart_prog =  progenitor_mass(iprog) / particle_mass
          npart_prog =  progenitor_mass(iprog) / particle_mass
          npart_desc =  clmp_mass_pb(ipeak) / particle_mass
          merit_calc =  DBLE(prog_desc_links(iprog, ipeak))/ npart_desc
          !merit_calc =  DBLE(prog_desc_links(iprog, ipeak))  * npart_prog/ npart_desc ** 2

          if (merit_calc > merit_max) then
            merit_max = merit_calc
            merit_max_id = ipeak
            main_desc_merit(iprog) = merit_calc
          end if
          
        end if
      end do


      ! TODO: if necessary, insert second condition here
      if (merit_max_id > 0) then 
        main_desc(iprog) = merit_max_id
        !check for next peak only if no other progenitor found;
        !otherwise, repeat
        iprog = iprog - 1 
        is_first = .false.
      else
        if (debug) write(*,*) "DEBUG 2-10-2 store id in use", store_id
        main_desc(iprog) = -store_id
      end if

    end if

    iprog = iprog + 1

  end do











  !===============================================
  ! Create tree:
  ! Check descendants for fracturing, identify main
  ! descendant
  !===============================================




  ipeak = 1

  do while (ipeak <= ncpu*npeaks_max)
   
    if (main_prog(ipeak) > 0) then
      ! if this descendant is not main descendant
      ! of its own main progenitor
      thisprog = prog_index(main_prog(ipeak))
      if (ipeak /= main_desc(thisprog)) then

        prog_desc_links(thisprog, ipeak) = -1
        merit_max = 0
        merit_max_id = 0 
        
        !try finding another
        do iprog = 1, nprogs
         
          if (prog_desc_links(iprog, ipeak) > 0) then
          
            if (debug) write(*,*) "DEBUG 2-12-1 mass", clmp_mass_pb(ipeak), "pm", particle_mass, "particle links", prog_desc_links(iprog, ipeak)
            if (debug) write(*,*) "DEBUG 2-12-2 iprog", iprog, "ipeak", ipeak 
        
            npart_desc =  clmp_mass_pb(ipeak) / particle_mass
            merit_calc =  DBLE(prog_desc_links(iprog, ipeak)) / npart_desc
        
            if (merit_calc > merit_max) then
              merit_max = merit_calc
              merit_max_id = sorted_progs(iprog)
              main_desc_merit(iprog) = merit_calc
            end if
            
          end if
        end do


        ! TODO: if necessary, insert second condition here
        if (merit_max_id > 0) then 
          main_prog(ipeak) = merit_max_id
          !check for next peak only if no other progenitor found;
          !otherwise, repeat
          ipeak = ipeak - 1 
        else
          if (writeinfo) write(*,*) "DEBUG 2-13 found new halo ", ipeak
          main_prog(ipeak) = -6
        end if

      end if
    end if

    ipeak = ipeak + 1

  end do





  !---------------------------------
  ! write found pairs to screen
  !---------------------------------



  if (writeinfo) then
    write(*,*)
    write(*,*) "LOG 1-1 Identified descendants:"
    do iprog = 1, nprogs
      write(*,*) "Progenitor", sorted_progs(iprog), "has main descendant", main_desc(iprog), "merit", main_desc_merit(iprog)
    end do

    write(*,*)
    write(*,*) "LOG 1-2 Identified progenitors:"
    do ipeak = 1, ncpu*npeaks_max
      if (clmp_mass_pb(ipeak) > 0 ) then
        write(*,*) "Descendant", ipeak, "has main progenitor", main_prog(ipeak), "merit", main_prog_merit(ipeak)
      end if
    end do

  end if

  deallocate(tracercount)
  deallocate(main_prog_merit, main_desc_merit)

end subroutine create_links 
!################################################################
!################################################################
!################################################################
!=============================================
subroutine deallocate_mergertree_loop
!=============================================
  use mergertree_commons
  implicit none
  deallocate(halocount)
  deallocate(progenitor_list, progenitor_tracers)
  deallocate(prog_desc_links)
  deallocate(prog_index)
  deallocate(progenitor_mass)
  deallocate(sorted_progs)
  deallocate(main_prog, main_desc)


  deallocate(npeaks, hfree, clmp_mass_pb)
  deallocate(clmpidp)
  deallocate(levelp)
  deallocate(idp)
  deallocate(pind_by_id)

end subroutine deallocate_mergertree_loop
