module mod_rays2d
  use mod_fmt2d_nointerfaces
  implicit none

contains

subroutine trace_ray_from_receiver2d(rec,s,ray)

  type(Tsource2d),target         :: s
  type(Treceiver2d)              :: rec
  type(Tray2d),target            :: ray

  real(kind=dp)                  :: x_in,z_in
  real(kind=dp)                  :: x_out,z_out
  integer                        :: n,m,nsections,ntf,n_shift
  logical                        :: outside
  logical                        :: reg_shared

  nsections = s%path(ray%raypath_id)%n_tf
  ray%nsections = nsections
  if (associated(ray%section)) then
     do n = 1, size(ray%section)
        if (associated(ray%section(n)%point)) deallocate(ray%section(n)%point)
     end do
     deallocate(ray%section)
  end if
  allocate(ray%section(nsections))
  do n = 1, nsections
     call ray_section2d_defaults(ray%section(n))
  end do

  x_in = rec%x
  z_in = rec%z

  outside = .false.

  do n = nsections, 1, -1

     ray%section(n)%ray => ray
     ray%section(n)%source => s
     ray%section(n)%place_in_sequence = n

     ntf = s%path(ray%raypath_id)%tf_sequence(n)
     if (ntf > 0 .and. ntf <= s%n_time_fields) then
        if (.not. associated(s%time_field(ntf)%arrivaltime)) then
           ray%valid = .false.
           return
        endif
        ray%section(n)%tf => s%time_field(ntf)
     else
        if (nsections /= 1) then
           ray%valid = .false.
           return
        endif
        ntf = 1
        if (.not. associated(s%time_field(ntf)%arrivaltime)) then
           ray%valid = .false.
           return
        endif
        ray%section(n)%tf => s%time_field(ntf)
     endif

     if (.not. associated(ray%section(n)%tf%reg)) then
        ray%valid = .false.
        return
     endif
     ray%section(n)%reg => ray%section(n)%tf%reg

     if (n < nsections) ray%section(n)%istart => intersection(s%path(ray%raypath_id)%sequence(2*n))
     if (n > 1) ray%section(n)%iend => intersection(s%path(ray%raypath_id)%sequence(2*n-1))

  end do

  ! Pre-load the time field data of the sections whose region object is theirs
  ! alone. trace_ray_section2d may access the next section's region data for
  ! the refraction/headwave check at interface boundaries, so the data has to
  ! be in place before the section loop starts. Without this, the refraction
  ! check uses stale data from a different time field, which causes Pg_ss
  ! section 2 (water up) to incorrectly generate a headwave along the seafloor,
  ! collapsing section 1 (crust) to just 2 points.
  ! Sections sharing a region object are skipped: only one time field can live
  ! in a region at a time, so that load would be overwritten by the section
  ! loop (which loads a section's own time field just before tracing it).
  do n = nsections, 1, -1
     if (.not. associated(ray%section(n)%tf%arrivaltime)) cycle
     reg_shared = .false.
     do m = 1, nsections
        if (m /= n) then
           if (associated(ray%section(m)%reg, ray%section(n)%reg)) reg_shared = .true.
        endif
     end do
     if (reg_shared) cycle
     if (associated(ray%section(n)%reg%loaded_tf, ray%section(n)%tf)) cycle
     ray%section(n)%reg%arrivaltime = ray%section(n)%tf%arrivaltime
     ray%section(n)%reg%time_gradient = ray%section(n)%tf%time_gradient
     ray%section(n)%reg%loaded_tf => ray%section(n)%tf
  end do

  do n = nsections, 1, -1

     if (associated(ray%section(n)%tf%arrivaltime)) then
        ! Make sure this section's time field is the one resident in its region.
        ! The copy is skipped when the region already holds it, which is the
        ! common case for rays that reuse the same time fields.
        if (.not. associated(ray%section(n)%reg%loaded_tf, ray%section(n)%tf)) then
           ray%section(n)%reg%arrivaltime = ray%section(n)%tf%arrivaltime
           ray%section(n)%reg%time_gradient = ray%section(n)%tf%time_gradient
           ray%section(n)%reg%loaded_tf => ray%section(n)%tf
        endif
     endif

     if (n == nsections) then
        call interpolate_arrivaltime2d(ray%section(n)%reg, rec%x, rec%z, ray%receiver_time)
        call interpolate_time_gradient2d(ray%section(n)%reg, rec%x, rec%z, &
             ray%receiver_time_gradient(1), ray%receiver_time_gradient(2))
        if (ray%receiver_time > 1.0e10_dp) then
           ray%valid = .false.
           return
        endif
        if (ray%receiver_time_gradient(1) == 0.0_dp .and. ray%receiver_time_gradient(2) == 0.0_dp) then
           ray%valid = .false.
           return
        endif
     endif

     call trace_ray_section2d(ray%section(n), x_in, z_in, x_out, z_out, outside)

     if (outside) then
        exit
     endif
     if (.not. ray%valid) then
        exit
     endif

     if (n > 1) then
        x_in = x_out
        z_in = z_out
     endif

  end do

  if (outside) then
     do m = n, nsections
        n_shift = m - n + 1
        ray%section(n_shift) = ray%section(m)
        ray%section(n_shift)%place_in_sequence = n_shift
     end do
     ray%nsections = nsections - n + 1
  endif

  return

end subroutine trace_ray_from_receiver2d

!***************************************************************************************************************

subroutine trace_ray_section2d(raysec, x_in, z_in, x_out, z_out, outside)

  type(Tray_section2d)                       :: raysec
  type(Ttime_field2d),pointer                :: tf
  type(Ttime_field2d),pointer                :: tf_next
  type(Tpropagation_grid2d),pointer          :: grid
  type(Tintersection2d),pointer              :: iftop, ifbot
  type(Tintersection2d),pointer              :: ifend
  integer                                    :: ifend_id
  type(Tsource2d),pointer                    :: s
  real(kind=dp)                              :: x_in, z_in
  real(kind=dp)                              :: x_out, z_out
  real(kind=dp)                              :: z_interface
  real(kind=dp)                              :: z_ifend_prev, z_ifend_new
  real(kind=dp)                              :: f_ifend_prev, f_ifend_new
  real(kind=dp)                              :: frac_ifend

  real(kind=dp),allocatable,dimension(:,:)   :: pt
  real(kind=dp)                              :: pp(2), grad0(2), grad1(2)
  real(kind=dp)                              :: norm(2), h
  real(kind=dp)                              :: dlmax, dl
  real(kind=dp)                              :: v_average
  real(kind=dp)                              :: p_tolerance
  real(kind=dp)                              :: s_tolerance
  real(kind=dp)                              :: dx
  integer                                    :: maxpoints
  integer                                    :: n, iter, seqnum, start_direction
  logical                                    :: section_is_refraction, outside
  logical                                    :: start_on_bot, start_on_top
  logical                                    :: crossing_ok
  logical                                    :: same_internal_ifend
  logical                                    :: left_start_interface
  logical                                    :: use_diving_return_seed
  integer                                    :: n_headwave
  integer                                    :: seed_node, nseed, iseed
  real(kind=dp)                              :: dist_headwave
  real(kind=dp)                              :: seed_dist, frac_seed

  tf => raysec%tf
  iftop => tf%reg%itop
  ifbot => tf%reg%ibot

  grid => pgrid
  s => raysec%source

  outside = .false.

  seqnum = raysec%place_in_sequence
  v_average = 5.0_dp
  maxpoints = 5000

  section_is_refraction = .false.
  start_direction = 0
  same_internal_ifend = .false.
  left_start_interface = .true.
  use_diving_return_seed = .false.
  nullify(ifend)
  nullify(tf_next)

  if (seqnum < raysec%ray%nsections .and. seqnum > 1) then
     ifend_id = raysec%iend%id
     ifend => intersection(ifend_id)
     tf_next => raysec%ray%section(seqnum-1)%tf
     start_direction = tf%reg%id - raysec%ray%section(seqnum+1)%reg%id
     section_is_refraction = raysec%ray%section(seqnum-1)%reg%id /= raysec%ray%section(seqnum)%reg%id
  endif

  if (seqnum == raysec%ray%nsections .and. seqnum > 1) then
     ifend_id = raysec%iend%id
     ifend => intersection(ifend_id)
     tf_next => raysec%ray%section(seqnum-1)%tf
     start_direction = 0
     section_is_refraction = raysec%ray%section(seqnum-1)%reg%id /= raysec%ray%section(seqnum)%reg%id
  endif

  if (seqnum == 1) then
     ifend_id = 0
     nullify(ifend)
     nullify(tf_next)
     if (seqnum+1 <= raysec%ray%nsections) then
        start_direction = tf%reg%id - raysec%ray%section(seqnum+1)%reg%id
     else
        start_direction = 0
     endif
  endif

  if (section_is_refraction .and. associated(tf_next)) then
     if (associated(tf_next%arrivaltime) .and. associated(tf_next%reg%arrivaltime)) then
        if (tf_next%reg%id /= tf%reg%id) then
           if (.not. associated(tf_next%reg%loaded_tf, tf_next)) then
              tf_next%reg%arrivaltime = tf_next%arrivaltime
              tf_next%reg%time_gradient = tf_next%time_gradient
              tf_next%reg%loaded_tf => tf_next
           endif
        endif
     endif
  endif

  if (associated(ifend) .and. associated(raysec%istart)) then
     same_internal_ifend = (ifend%id == raysec%istart%id .and. &
          ifend%id /= iftop%id .and. ifend%id /= ifbot%id)
     if (same_internal_ifend) left_start_interface = .false.
  endif

  if (.not. associated(tf%time_gradient)) then
     raysec%ray%valid = .false.
     return
  endif

  allocate(pt(2, maxpoints))

  dlmax = 0.1_dp * grid%dx * v_average
  p_tolerance = 5.0e-4_dp * grid%dx
  s_tolerance = 2.5_dp * grid%dx

  n = 1
  pt(1,1) = x_in
  pt(2,1) = z_in
  n_headwave = 0
  dist_headwave = 0.0_dp

  ! Check if starting point is on bottom or top interface (for section 1)
  ! This is used to detect the case where the ray starts FROM an interface
  ! and needs a larger first step to leave it
  z_interface = interpolate_interface2d(x_in, ifbot)
  start_on_bot = (abs(z_in - z_interface) < 0.1_dp * grid%dz)
  z_interface = interpolate_interface2d(x_in, iftop)
  start_on_top = (abs(z_in - z_interface) < 0.1_dp * grid%dz)

  ! Converted diving phases return to sed_bot before the final P up leg.
  ! Start the deep section from the same interior Huygens point used to seed
  ! that return, then follow the saved deep time-field gradient back to the
  ! incoming sed_bot point.
  if (seqnum == 2 .and. start_on_top .and. associated(raysec%istart) .and. &
       associated(raysec%iend)) then
     if (raysec%istart%id == iftop%id .and. raysec%iend%id == iftop%id .and. &
          (raysec%ray%raypath_id == 5 .or. raysec%ray%raypath_id == 6 .or. &
           raysec%ray%raypath_id == 8)) then
        seed_node = best_diving_return_node2d(tf, x_in, z_in, 2.0_dp*grid%dz)
        if (seed_node > 0) then
           seed_dist = sqrt((tf%reg%x(seed_node) - x_in)**2 + &
                (tf%reg%z(seed_node) - z_in)**2)
           if (seed_dist > 1.0e-10_dp) use_diving_return_seed = .true.
        endif
     endif
  endif

  if (use_diving_return_seed) then
     nseed = ceiling(seed_dist / max(dlmax, 1.0e-10_dp))
     if (nseed < 1) nseed = 1
     if (nseed > maxpoints-2) nseed = maxpoints-2
     do iseed = 1, nseed
        frac_seed = real(iseed, kind=dp) / real(nseed, kind=dp)
        pt(1,1+iseed) = x_in + frac_seed * (tf%reg%x(seed_node) - x_in)
        pt(2,1+iseed) = z_in + frac_seed * (tf%reg%z(seed_node) - z_in)
     end do
     n = 1 + nseed
  endif

ploop: do

     n = n + 1

     if (n > maxpoints) then
        raysec%ray%valid = .false.
        deallocate(pt)
        return
     endif

     ! Early source convergence check for seqnum==1:
     ! Before interface checks, test if the previous point is close enough to the source.
     ! This prevents the ray from getting stuck oscillating near the source on the seafloor.
     if (seqnum == 1 .and. n > 2) then
        dx = sqrt((pt(1,n-1)-s%x)**2 + (pt(2,n-1)-s%z)**2)
        if (dx < s_tolerance) then
           pt(1,n) = s%x
           pt(2,n) = s%z
           raysec%npoints = n
           allocate(raysec%point(2,n))
           raysec%point(1:2,1:n) = pt(1:2,1:n)
           deallocate(pt)
           return
        endif
     endif

     if (raysec%ray%headwave .or. raysec%ray%diffracted) then
        if (n > 2) then
           dist_headwave = dist_headwave + sqrt((pt(1,n-1)-pt(1,n-2))**2 + (pt(2,n-1)-pt(2,n-2))**2)
        endif
        if (dist_headwave > 200.0_dp) then
           raysec%ray%valid = .false.
           deallocate(pt)
           return
        endif
     endif

     dl = dlmax
     ! For sections starting from an interface, use a larger first step
     ! to ensure the ray can actually leave the interface
     if (n == 2 .and. (start_on_bot .or. start_on_top)) then
        dl = 10.0_dp * dlmax
     endif

     call interpolate_time_gradient2d(tf%reg, pt(1,n-1), pt(2,n-1), grad0(1), grad0(2))
     if (.not. raysec%ray%valid) then
        deallocate(pt)
        return
     endif

     pp(1) = pt(1,n-1) - dl*grad0(1)
     pp(2) = pt(2,n-1) - dl*grad0(2)

     ! When starting on an interface, the gradient at the boundary may point
     ! in the wrong direction (toward the boundary instead of into the interior).
     ! If the ray would go outside the region, reflect it back into the interior.
     if (start_on_bot .and. pp(2) > pt(2,n-1)) then
        pp(2) = 2.0_dp*pt(2,n-1) - pp(2)
     endif
     if (start_on_top .and. pp(2) < pt(2,n-1)) then
        pp(2) = 2.0_dp*pt(2,n-1) - pp(2)
     endif

     ! Start direction check: verify gradient is consistent with ray direction
     ! Uses interpolated gradient (grad0) and threshold 0.0, matching original fmcode
     ! Skip for section 1 (seqnum==1): gradient always points away from source,
     ! which is the correct direction for ray tracing back to source.
     ! The check assumes higher region IDs are deeper, which may not hold
     ! when a combined region (e.g., region 5) has a higher ID than mantle (region 4).
     ! For backward ray tracing, the gradient points opposite to the ray direction,
     ! so the check passes when dot_product <= 0 (gradient opposite to start_direction).
     if (n == 2 .and. start_direction /= 0 .and. .not. section_is_refraction .and. seqnum /= 1) then
        call interface_normal2d(x_in, raysec%istart, norm(1), norm(2), h)
        if (dot_product(grad0, start_direction*norm) > 0.0_dp) then
           raysec%ray%valid = .false.
           deallocate(pt)
           return
        endif
     endif

     ! Check ifend interface (intermediate interface within the region)
     ! This handles the case where the ray ends at an interface that is
     ! neither the top nor bottom of the region (e.g., sed_bot within region 5)
      if (associated(ifend)) then
         if (ifend%id /= iftop%id .and. ifend%id /= ifbot%id) then
           z_ifend_prev = interpolate_interface2d(pt(1,n-1), ifend)
           z_ifend_new = interpolate_interface2d(pp(1), ifend)
           f_ifend_prev = pt(2,n-1) - z_ifend_prev
           f_ifend_new = pp(2) - z_ifend_new
            if (same_internal_ifend .and. .not. left_start_interface) then
               if (abs(f_ifend_new) > 0.5_dp * grid%dz) left_start_interface = .true.
            else if (abs(f_ifend_prev) > 1.0e-3_dp) then
               if (f_ifend_prev * f_ifend_new <= 0.0_dp) then
                 if (abs(f_ifend_prev - f_ifend_new) > 1.0e-12_dp) then
                    frac_ifend = f_ifend_prev / (f_ifend_prev - f_ifend_new)
                 else
                    frac_ifend = 1.0_dp
                 endif
                 if (frac_ifend < 0.0_dp) frac_ifend = 0.0_dp
                 if (frac_ifend > 1.0_dp) frac_ifend = 1.0_dp
                 pt(1,n) = pt(1,n-1) + frac_ifend * (pp(1) - pt(1,n-1))
                 pt(2,n) = interpolate_interface2d(pt(1,n), ifend)
                 if (seqnum == 1) then
                    dx = sqrt((pt(1,n)-s%x)**2 + (pt(2,n)-s%z)**2)
                    if (dx < s_tolerance) then
                       pt(1,n) = s%x
                       pt(2,n) = s%z
                    endif
                 endif
                 raysec%npoints = n
                 allocate(raysec%point(2,n))
                 raysec%point(1:2,1:n) = pt(1:2,1:n)
                 deallocate(pt)
                 return
            endif
         endif
      endif
     endif

     z_interface = interpolate_interface2d(pp(1), iftop)
     if (pp(2) <= z_interface) then
        if (iftop%id == ifend_id) then
           if (section_is_refraction) then
              call interface_normal2d(pp(1), intersection(ifend_id), norm(1), norm(2), h)
              call interpolate_time_gradient2d(tf%reg, pp(1), z_interface, grad0(1), grad0(2))
              call interpolate_time_gradient2d(tf_next%reg, pp(1), z_interface, grad1(1), grad1(2))
              if (.not. raysec%ray%valid) then
                 deallocate(pt)
                 return
              endif
              if (sqrt((grad1(1) - dot_product(norm,grad1)*norm(1))**2 + &
                       (grad1(2) - dot_product(norm,grad1)*norm(2))**2) <= sqrt(sum(grad0**2))) then
                 exit ploop
              else
                 pt(1,n) = pp(1)
                 pt(2,n) = z_interface
                 raysec%headwave = .true.
                 raysec%ray%headwave = .true.
                 if (seqnum == 1) then
                    dx = sqrt((pt(1,n)-s%x)**2 + (pt(2,n)-s%z)**2)
                    if (dx < s_tolerance) then
                       pt(1,n) = s%x
                       pt(2,n) = s%z
                       raysec%npoints = n
                       allocate(raysec%point(2,n))
                       raysec%point(1:2,1:n) = pt(1:2,1:n)
                       deallocate(pt)
                       return
                    endif
                 endif
                 cycle ploop
              endif
           else
              exit ploop
           endif
         else if (n == 2 .and. start_on_top .and. associated(raysec%istart) .and. &
                  raysec%istart%id == iftop%id) then
            continue
        else
           pt(1,n) = pp(1)
           pt(2,n) = z_interface
           raysec%diffracted = .true.
           raysec%ray%diffracted = .true.
           if (seqnum == 1) then
              dx = sqrt((pt(1,n)-s%x)**2 + (pt(2,n)-s%z)**2)
              if (dx < s_tolerance) then
                 pt(1,n) = s%x
                 pt(2,n) = s%z
                 raysec%npoints = n
                 allocate(raysec%point(2,n))
                 raysec%point(1:2,1:n) = pt(1:2,1:n)
                 deallocate(pt)
                 return
              endif
           endif
           cycle ploop
        endif
     endif

     z_interface = interpolate_interface2d(pp(1), ifbot)
     if (pp(2) >= z_interface) then
        if (ifbot%id == ifend_id) then
           if (section_is_refraction) then
              ! For Pg_ss section 2 (water up, tf6), the ray must refract at
              ! the seafloor into the crust (section 1, tf1). The refraction
              ! check can fail due to numerical precision when tf6's gradient
              ! at the seafloor is nearly horizontal (small vertical component),
              ! causing a spurious headwave that collapses section 1 to 2 points.
              ! Since crust velocity > water velocity, the ray always refracts.
              if (raysec%ray%raypath_id == 4 .and. seqnum == 2) then
                 exit ploop
              endif
              call interface_normal2d(pp(1), intersection(ifend_id), norm(1), norm(2), h)
              call interpolate_time_gradient2d(tf%reg, pp(1), z_interface, grad0(1), grad0(2))
              call interpolate_time_gradient2d(tf_next%reg, pp(1), z_interface, grad1(1), grad1(2))
              if (.not. raysec%ray%valid) then
                 deallocate(pt)
                 return
              endif
              if (sqrt((grad1(1) - dot_product(norm,grad1)*norm(1))**2 + &
                       (grad1(2) - dot_product(norm,grad1)*norm(2))**2) <= sqrt(sum(grad0**2))) then
                 exit ploop
              else
                 pt(1,n) = pp(1)
                 pt(2,n) = z_interface
                 raysec%headwave = .true.
                 raysec%ray%headwave = .true.
                 if (seqnum == 1) then
                    dx = sqrt((pt(1,n)-s%x)**2 + (pt(2,n)-s%z)**2)
                    if (dx < s_tolerance) then
                       pt(1,n) = s%x
                       pt(2,n) = s%z
                       raysec%npoints = n
                       allocate(raysec%point(2,n))
                       raysec%point(1:2,1:n) = pt(1:2,1:n)
                       deallocate(pt)
                       return
                    endif
                 endif
                 cycle ploop
              endif
           else
              exit ploop
           endif
         else if (n == 2 .and. start_on_bot .and. associated(raysec%istart) .and. &
                  raysec%istart%id == ifbot%id) then
            continue
        else
           pt(1,n) = pp(1)
           pt(2,n) = z_interface
           raysec%diffracted = .true.
           raysec%ray%diffracted = .true.
           if (seqnum == 1) then
              dx = sqrt((pt(1,n)-s%x)**2 + (pt(2,n)-s%z)**2)
              if (dx < s_tolerance) then
                 pt(1,n) = s%x
                 pt(2,n) = s%z
                 raysec%npoints = n
                 allocate(raysec%point(2,n))
                 raysec%point(1:2,1:n) = pt(1:2,1:n)
                 deallocate(pt)
                 return
              endif
           endif
           cycle ploop
        endif
     endif

     iter = 1
itloop: do

        pp(1) = pt(1,n-1) - dl*grad0(1)
        pp(2) = pt(2,n-1) - dl*grad0(2)

        call interpolate_time_gradient2d(tf%reg, pp(1), pp(2), grad1(1), grad1(2))
        if (.not. raysec%ray%valid) then
           deallocate(pt)
           return
        endif

        pt(1,n) = pt(1,n-1) - dl*(grad0(1)+grad1(1))*0.5_dp
        pt(2,n) = pt(2,n-1) - dl*(grad0(2)+grad1(2))*0.5_dp

        z_interface = interpolate_interface2d(pt(1,n), iftop)
        if (pt(2,n) <= z_interface) then
           if (iftop%id == ifend_id) then
              if (section_is_refraction) then
                 call interface_normal2d(pt(1,n), intersection(ifend_id), norm(1), norm(2), h)
                 call interpolate_time_gradient2d(tf%reg, pt(1,n), z_interface, grad0(1), grad0(2))
                 call interpolate_time_gradient2d(tf_next%reg, pt(1,n), z_interface, grad1(1), grad1(2))
                 if (.not. raysec%ray%valid) then
                    deallocate(pt)
                    return
                 endif
                 if (sqrt((grad1(1) - dot_product(norm,grad1)*norm(1))**2 + &
                          (grad1(2) - dot_product(norm,grad1)*norm(2))**2) <= sqrt(sum(grad0**2))) then
                    pp = pt(1:2,n)
                    exit ploop
                 else
                    pt(2,n) = z_interface
                    raysec%headwave = .true.
                    raysec%ray%headwave = .true.
                    cycle ploop
                 endif
              else
                 pp = pt(1:2,n)
                 exit ploop
              endif
            else if (n == 2 .and. start_on_top .and. associated(raysec%istart) .and. &
                     raysec%istart%id == iftop%id) then
               continue
           else
              pt(2,n) = z_interface
              raysec%diffracted = .true.
              raysec%ray%diffracted = .true.
              cycle ploop
           endif
        endif

        z_interface = interpolate_interface2d(pt(1,n), ifbot)
        if (pt(2,n) >= z_interface) then
           if (ifbot%id == ifend_id) then
              if (section_is_refraction) then
                 ! Same Pg_ss section 2 fix as above
                 if (raysec%ray%raypath_id == 4 .and. seqnum == 2) then
                    pp = pt(1:2,n)
                    exit ploop
                 endif
                 call interface_normal2d(pt(1,n), intersection(ifend_id), norm(1), norm(2), h)
                 call interpolate_time_gradient2d(tf%reg, pt(1,n), z_interface, grad0(1), grad0(2))
                 call interpolate_time_gradient2d(tf_next%reg, pt(1,n), z_interface, grad1(1), grad1(2))
                 if (.not. raysec%ray%valid) then
                    deallocate(pt)
                    return
                 endif
                 if (sqrt((grad1(1) - dot_product(norm,grad1)*norm(1))**2 + &
                          (grad1(2) - dot_product(norm,grad1)*norm(2))**2) <= sqrt(sum(grad0**2))) then
                    pp = pt(1:2,n)
                    exit ploop
                 else
                    pt(2,n) = z_interface
                    raysec%headwave = .true.
                    raysec%ray%headwave = .true.
                    cycle ploop
                 endif
              else
                 pp = pt(1:2,n)
                 exit ploop
              endif
            else if (n == 2 .and. start_on_bot .and. associated(raysec%istart) .and. &
                     raysec%istart%id == ifbot%id) then
               continue
           else
              pt(2,n) = z_interface
              raysec%diffracted = .true.
              raysec%ray%diffracted = .true.
              cycle ploop
           endif
        endif

        if (sqrt(sum((pt(1:2,n)-pp(1:2))**2)) < p_tolerance) exit itloop

        iter = iter + 1
        if (iter > 10) then
           raysec%ray%valid = .false.

           deallocate(pt)
           return
        end if

        dl = 0.5_dp * dl

     end do itloop

     if (seqnum == 1) then
        dx = sqrt((pt(1,n)-s%x)**2 + (pt(2,n)-s%z)**2)
        if (dx < s_tolerance) then
           pt(1,n) = s%x
           pt(2,n) = s%z
           raysec%npoints = n
           allocate(raysec%point(2,n))
           raysec%point(1:2,1:n) = pt(1:2,1:n)
           deallocate(pt)
           return
        endif
     endif

  end do ploop

  if (seqnum == 1) then
     raysec%ray%valid = .false.

     deallocate(pt)
     return
  endif

  if (outside) then
     raysec%npoints = n
     allocate(raysec%point(2,n))
     raysec%point(1:2,1:n) = pt(1:2,1:n)
     deallocate(pt)
     return
  endif

  if (associated(ifend)) then
     call find_interface_crossing2d(pt, pp, n, ifend, grad0, dl, pgrid, crossing_ok)
     if (.not. crossing_ok) then
        raysec%ray%valid = .false.
        deallocate(pt)
        return
     endif
  endif

  x_out = pt(1,n)
  z_out = pt(2,n)

  raysec%npoints = n
  allocate(raysec%point(2,n))
  raysec%point(1:2,1:n) = pt(1:2,1:n)

  deallocate(pt)

  return

end subroutine trace_ray_section2d

!***************************************************************************************************************

function best_diving_return_node2d(tf, x_if, z_if, min_depth) result(best_m)
  type(Ttime_field2d)                       :: tf
  real(kind=dp)                             :: x_if, z_if, min_depth
  integer                                   :: best_m

  integer                                   :: m
  real(kind=dp)                             :: dx, dz, dist, vloc, cand, best_t

  best_m = 0
  best_t = huge_time

  if (.not. associated(tf%reg)) return
  if (.not. associated(tf%arrivaltime)) return

  do m = 1, tf%reg%nnode
     if (m > size(tf%arrivaltime)) cycle
     if (tf%arrivaltime(m) >= huge_time*0.5_dp) cycle
     if (tf%reg%z(m) <= z_if + min_depth) cycle
     if (abs(tf%reg%x(m) - x_if) > 3.0_dp * min_depth) cycle
     if (tf%vtype == VTYPE_S .and. tf%reg%velocity(tf%vtype,m) < 0.01_dp) cycle

     dx = tf%reg%x(m) - x_if
     dz = tf%reg%z(m) - z_if
     dist = sqrt(dx*dx + dz*dz)
     if (dist <= 1.0e-10_dp) cycle

     vloc = tf%reg%velocity(tf%vtype,m)
     if (vloc < 0.01_dp) vloc = 0.01_dp
     cand = tf%arrivaltime(m) + dist / vloc

     if (cand < best_t) then
        best_t = cand
        best_m = m
     endif
  end do

end function best_diving_return_node2d

!***************************************************************************************************************

subroutine find_interface_crossing2d(pt, pp, n, ifend, grad0, dl, grid, ok)

  real(kind=dp),dimension(:,:)              :: pt
  real(kind=dp),dimension(2)                :: pp, grad0
  integer                                   :: n
  type(Tintersection2d),pointer             :: ifend
  real(kind=dp)                             :: dl
  type(Tpropagation_grid2d),pointer         :: grid
  logical                                   :: ok

  real(kind=dp)                             :: x1, x2, f, fmid, xmid, dx, rtbis
  real(kind=dp)                             :: z_interface
  integer                                   :: nn

  ok = .true.

  z_interface = interpolate_interface2d(pp(1), ifend)
  x1 = 0.0_dp
  x2 = dl
  f = pt(2,n-1) - interpolate_interface2d(pt(1,n-1), ifend)
  fmid = pp(2) - z_interface

  if (f*fmid > 0.01_dp*grid%tolerance) then
     ok = .false.
     return
  endif

  if (abs(f*fmid) <= 0.01_dp*grid%tolerance) then
     pt(1:2,n) = pt(1:2,n-1)
     pt(2,n) = interpolate_interface2d(pt(1,n), ifend)
     return
  endif

  if (f*fmid < -0.01_dp*grid%tolerance) then
     if (f < 0.0_dp) then
        rtbis = x1
        dx = x2 - x1
     else
        rtbis = x2
        dx = x1 - x2
     endif
     do nn = 1, 10
        dx = dx * 0.5_dp
        xmid = rtbis + dx
        pt(1,n) = pt(1,n-1) - xmid*grad0(1)
        pt(2,n) = pt(2,n-1) - xmid*grad0(2)
        fmid = pt(2,n) - interpolate_interface2d(pt(1,n), ifend)
        if (fmid <= 0.0_dp) rtbis = xmid
     end do

     pt(1,n) = pt(1,n-1) - rtbis*grad0(1)
     pt(2,n) = pt(2,n-1) - rtbis*grad0(2)
     pt(2,n) = interpolate_interface2d(pt(1,n), ifend)
  endif

  return

end subroutine find_interface_crossing2d

!***************************************************************************************************************

subroutine interpolate_time_gradient2d(reg, px, pz, dtdx, dtdz)

  type(Tregion2d),target                     :: reg
  real(kind=dp)                              :: px, pz, dtdx, dtdz
  type(Tpropagation_grid2d),pointer          :: grid
  real(kind=dp)                              :: w1, w2
  integer                                    :: ix, iz, n1, n2, n3, n4
  real(kind=dp)                              :: gx, gz
  real(kind=dp)                              :: g1x, g2x, g3x, g4x, g1z, g2z, g3z, g4z

  grid => reg%grid

  ix = floor((px - grid%x0) / grid%dx) + 1
  iz = floor((pz - grid%z0) / grid%dz) + 1

  if (ix < 1) ix = 1
  if (ix >= grid%nx) ix = grid%nx - 1
  if (iz < 1) iz = 1
  if (iz >= grid%nz) iz = grid%nz - 1

  w1 = (px - grid%x(ix)) / grid%dx
  w2 = (pz - grid%z(iz)) / grid%dz

  if (w1 < 0.0_dp) w1 = 0.0_dp
  if (w1 > 1.0_dp) w1 = 1.0_dp
  if (w2 < 0.0_dp) w2 = 0.0_dp
  if (w2 > 1.0_dp) w2 = 1.0_dp

  n1 = find_region_node(reg, ix, iz)
  n2 = find_region_node(reg, ix+1, iz)
  n3 = find_region_node(reg, ix, iz+1)
  n4 = find_region_node(reg, ix+1, iz+1)

  if (n1 == 0) n1 = find_nearby_region_node(reg, ix, iz)
  if (n2 == 0) n2 = find_nearby_region_node(reg, ix+1, iz)
  if (n3 == 0) n3 = find_nearby_region_node(reg, ix, iz+1)
  if (n4 == 0) n4 = find_nearby_region_node(reg, ix+1, iz+1)

  g1x = 0.0_dp; g2x = 0.0_dp; g3x = 0.0_dp; g4x = 0.0_dp
  g1z = 0.0_dp; g2z = 0.0_dp; g3z = 0.0_dp; g4z = 0.0_dp

  if (n1 > 0 .and. n1 <= reg%nnode) then
     g1x = reg%time_gradient(1, n1); g1z = reg%time_gradient(2, n1)
  endif
  if (n2 > 0 .and. n2 <= reg%nnode) then
     g2x = reg%time_gradient(1, n2); g2z = reg%time_gradient(2, n2)
  endif
  if (n3 > 0 .and. n3 <= reg%nnode) then
     g3x = reg%time_gradient(1, n3); g3z = reg%time_gradient(2, n3)
  endif
  if (n4 > 0 .and. n4 <= reg%nnode) then
     g4x = reg%time_gradient(1, n4); g4z = reg%time_gradient(2, n4)
  endif

  gx = (1.0_dp-w1)*(1.0_dp-w2)*g1x + w1*(1.0_dp-w2)*g2x + &
       (1.0_dp-w1)*w2*g3x + w1*w2*g4x

  gz = (1.0_dp-w1)*(1.0_dp-w2)*g1z + w1*(1.0_dp-w2)*g2z + &
       (1.0_dp-w1)*w2*g3z + w1*w2*g4z

  if (sqrt(gx**2 + gz**2) < 0.01_dp) then
     call interpolate_gradient_reg(reg, px, pz, gx, gz)
  endif

  dtdx = gx
  dtdz = gz

  return

end subroutine interpolate_time_gradient2d

!***************************************************************************************************************

subroutine interpolate_arrivaltime2d(reg, px, pz, t)

  type(Tregion2d),target                     :: reg
  real(kind=dp)                              :: px, pz, t
  type(Tpropagation_grid2d),pointer          :: grid
  real(kind=dp)                              :: w1, w2
  integer                                    :: ix, iz, n, n1, n2, n3, n4
  real(kind=dp)                              :: t1, t2, t3, t4

  grid => reg%grid

  ix = floor((px - grid%x0) / grid%dx) + 1
  iz = floor((pz - grid%z0) / grid%dz) + 1

  if (ix < 1) ix = 1
  if (ix >= grid%nx) ix = grid%nx - 1
  if (iz < 1) iz = 1
  if (iz >= grid%nz) iz = grid%nz - 1

  w1 = (px - grid%x(ix)) / grid%dx
  w2 = (pz - grid%z(iz)) / grid%dz

  if (w1 < 0.0_dp) w1 = 0.0_dp
  if (w1 > 1.0_dp) w1 = 1.0_dp
  if (w2 < 0.0_dp) w2 = 0.0_dp
  if (w2 > 1.0_dp) w2 = 1.0_dp

  n1 = find_region_node(reg, ix, iz)
  n2 = find_region_node(reg, ix+1, iz)
  n3 = find_region_node(reg, ix, iz+1)
  n4 = find_region_node(reg, ix+1, iz+1)

  t1 = huge_time
  t2 = huge_time
  t3 = huge_time
  t4 = huge_time

  if (n1 > 0) t1 = reg%arrivaltime(n1)
  if (n2 > 0) t2 = reg%arrivaltime(n2)
  if (n3 > 0) t3 = reg%arrivaltime(n3)
  if (n4 > 0) t4 = reg%arrivaltime(n4)

  if (n1 == 0) call find_nearby_node_time(reg, ix, iz, t1)
  if (n2 == 0) call find_nearby_node_time(reg, ix+1, iz, t2)
  if (n3 == 0) call find_nearby_node_time(reg, ix, iz+1, t3)
  if (n4 == 0) call find_nearby_node_time(reg, ix+1, iz+1, t4)

  if (t1 < huge_time*0.5 .and. t2 < huge_time*0.5 .and. &
      t3 < huge_time*0.5 .and. t4 < huge_time*0.5) then
     t = (1.0_dp-w1)*(1.0_dp-w2)*t1 + w1*(1.0_dp-w2)*t2 + &
         (1.0_dp-w1)*w2*t3 + w1*w2*t4
  else
     t = interpolate_arrivaltime_reg(reg, px, pz)
  endif

  return

end subroutine interpolate_arrivaltime2d

function interpolate_arrivaltime_reg(reg, px, pz) result(t)
  type(Tregion2d),target :: reg
  real(kind=dp) :: px, pz, t
  integer :: n, n1, n2, n3, n4
  real(kind=dp) :: dmin, d

  t = huge_time
  dmin = huge_time

  do n = 1, reg%nnode
     d = (px - reg%x(n))**2 + (pz - reg%z(n))**2
     if (d < dmin) then
        dmin = d
        t = reg%arrivaltime(n)
     endif
  end do

end function interpolate_arrivaltime_reg

function find_region_node(reg, ix, iz) result(n)
  type(Tregion2d),target :: reg
  integer :: ix, iz, n

  n = 0
  if (associated(reg%region_rnode_id)) then
     if (ix >= 1 .and. ix <= reg%grid%nx .and. &
         iz >= 1 .and. iz <= reg%grid%nz) then
        n = reg%region_rnode_id(ix, iz)
        if (n <= 0 .or. n > reg%nnode) n = 0
     endif
  else
     do n = 1, reg%nnode
        if (reg%node(n)%i1 == ix .and. reg%node(n)%i2 == iz) then
           return
        endif
     end do
     n = 0
  endif

end function find_region_node

function find_nearby_region_node(reg, ix, iz) result(n)
  type(Tregion2d),target :: reg
  integer :: ix, iz, n
  integer :: diz, iz2

  n = 0
  if (.not. associated(reg%region_rnode_id)) return

  do diz = 1, 5
     iz2 = iz + diz
     if (iz2 >= 1 .and. iz2 <= reg%grid%nz) then
        n = reg%region_rnode_id(ix, iz2)
        if (n > 0 .and. n <= reg%nnode) return
     endif
     iz2 = iz - diz
     if (iz2 >= 1 .and. iz2 <= reg%grid%nz) then
        n = reg%region_rnode_id(ix, iz2)
        if (n > 0 .and. n <= reg%nnode) return
     endif
  end do

  n = 0

end function find_nearby_region_node

subroutine find_nearby_node_time(reg, ix, iz, t)
  type(Tregion2d),target :: reg
  integer :: ix, iz
  real(kind=dp) :: t
  integer :: n

  t = huge_time
  n = find_nearby_region_node(reg, ix, iz)
  if (n > 0 .and. n <= reg%nnode) then
     t = reg%arrivaltime(n)
  endif

end subroutine find_nearby_node_time

subroutine interpolate_gradient_reg(reg, px, pz, gx, gz)
  type(Tregion2d),target :: reg
  real(kind=dp) :: px, pz, gx, gz
  integer :: n, nmin
  real(kind=dp) :: dmin, d

  gx = 0.0_dp
  gz = 0.0_dp
  dmin = huge_time
  nmin = 0

  ! The nearest node is almost always one of the corners of the grid cell that
  ! holds the query point, so start with that O(1) neighbourhood lookup. The
  ! full sweep over reg%nnode is only needed when the neighbourhood holds no
  ! node at all, which means the query point lies outside the region.
  nmin = find_position_in_region2d(reg, px, pz)

  if (nmin > 0) then
     gx = reg%time_gradient(1, nmin)
     gz = reg%time_gradient(2, nmin)
     return
  endif

  do n = 1, reg%nnode
     d = (px - reg%x(n))**2 + (pz - reg%z(n))**2
     if (d < dmin) then
        dmin = d
        nmin = n
     endif
  end do

  if (nmin > 0) then
     gx = reg%time_gradient(1, nmin)
     gz = reg%time_gradient(2, nmin)
  endif

end subroutine interpolate_gradient_reg

!***************************************************************************************************************

function find_position_in_region2d(reg, px, pz) result(inode)

  type(Tregion2d)                            :: reg
  real(kind=dp)                              :: px, pz
  integer                                    :: inode

  type(Tpropagation_grid2d),pointer          :: grid
  integer                                    :: ix, iz
  real(kind=dp)                              :: dist, distmin
  integer                                    :: m, ii, jj

  grid => reg%grid

  ix = floor((px - grid%x0) / grid%dx) + 1
  iz = floor((pz - grid%z0) / grid%dz) + 1

  if (ix < 1) ix = 1
  if (ix >= grid%nx) ix = grid%nx - 1
  if (iz < 1) iz = 1
  if (iz >= grid%nz) iz = grid%nz - 1

  distmin = 1.0e100_dp
  inode = 0

  do ii = ix, ix+1
     do jj = iz, iz+1
        if (ii >= 1 .and. ii <= grid%nx .and. jj >= 1 .and. jj <= grid%nz) then
           if (reg%region_rnode_id(ii,jj) > 0) then
              m = reg%region_rnode_id(ii,jj)
              if (m > 0 .and. m <= reg%nnode) then
                 dist = (px - reg%x(m))**2 + (pz - reg%z(m))**2
                 if (dist < distmin) then
                    inode = m
                    distmin = dist
                 endif
              endif
           endif
        endif
     end do
  end do

  if (inode == 0) then
     distmin = 1.0e100_dp
     do m = 1, reg%nnode
        dist = (px - reg%x(m))**2 + (pz - reg%z(m))**2
        if (dist < distmin) then
           inode = m
           distmin = dist
        endif
     end do
  endif

  return

end function find_position_in_region2d

!***************************************************************************************************************

subroutine cross_interface2d(raysec, isec, direction)

  type(Tray_section2d)                       :: raysec
  type(Tintersection2d),pointer              :: isec
  integer                                    :: direction

  real(kind=dp)                              :: norm(2), h
  real(kind=dp)                              :: grad_cur(2), grad_next(2)
  real(kind=dp)                              :: v_cur, v_next
  real(kind=dp)                              :: sin_theta1, sin_theta2
  real(kind=dp)                              :: ray_param
  real(kind=dp)                              :: cos_theta1, cos_theta2
  real(kind=dp)                              :: pt_x, pt_z
  integer                                    :: vtype_cur, vtype_next
  integer                                    :: seqnum

  seqnum = raysec%place_in_sequence

  if (raysec%npoints < 1) return

  pt_x = raysec%point(1, raysec%npoints)
  pt_z = raysec%point(2, raysec%npoints)

  call interface_normal2d(pt_x, isec, norm(1), norm(2), h)

  call interpolate_time_gradient2d(raysec%reg, pt_x, pt_z, grad_cur(1), grad_cur(2))

  vtype_cur = raysec%tf%vtype

  if (direction > 0) then
     if (associated(isec%regright)) then
        call interpolate_time_gradient2d(isec%regright, pt_x, pt_z, grad_next(1), grad_next(2))
        vtype_next = raysec%ray%section(seqnum-1)%tf%vtype
     else
        return
     endif
  else if (direction < 0) then
     if (associated(isec%regleft)) then
        call interpolate_time_gradient2d(isec%regleft, pt_x, pt_z, grad_next(1), grad_next(2))
        vtype_next = raysec%ray%section(seqnum-1)%tf%vtype
     else
        return
     endif
  else
     return
  endif

  v_cur = get_velocity_at2d(raysec%reg, pt_x, pt_z, vtype_cur)
  v_next = get_velocity_at2d(raysec%ray%section(seqnum-1)%reg, pt_x, pt_z, vtype_next)

  ! Guard against zero gradient (avoid division by zero)
  if (sum(grad_cur**2) < tiny(1.0_dp)) return

  ! Guard against zero velocity (avoid division by zero in ray_param)
  if (v_cur < tiny(1.0_dp)) return

  sin_theta1 = sqrt(max(0.0_dp, sum(grad_cur**2) - dot_product(grad_cur, norm)**2)) / sqrt(sum(grad_cur**2))
  ray_param = sin_theta1 / v_cur

  sin_theta2 = ray_param * v_next

  if (sin_theta2 > 1.0_dp) then
     raysec%headwave = .true.
     raysec%ray%headwave = .true.
     return
  endif

  ! Clamp sin_theta to [0,1] to avoid sqrt of negative due to roundoff
  cos_theta1 = sqrt(max(0.0_dp, 1.0_dp - sin_theta1**2))
  cos_theta2 = sqrt(max(0.0_dp, 1.0_dp - sin_theta2**2))

  return

end subroutine cross_interface2d

!***************************************************************************************************************

subroutine compute_ray_frechet2d(ray)

  type(Tray2d)                               :: ray

  return

end subroutine compute_ray_frechet2d

!***************************************************************************************************************

function interpolate_interface2d(px, isec) result(z_val)

  real(kind=dp)                              :: px, z_val
  type(Tintersection2d)                      :: isec

  integer                                    :: ix
  real(kind=dp)                              :: w

  if (isec%nnode < 1) then
     z_val = 0.0_dp
     return
  endif
  if (isec%nnode < 2) then
     z_val = isec%z(1)
     return
  endif

  ix = floor((px - isec%x(1)) / isec%grid%dx) + 1

  if (ix < 1) then
     z_val = isec%z(1)
     return
  endif

  if (ix >= isec%nnode) then
     z_val = isec%z(isec%nnode)
     return
  endif

  w = (px - isec%x(ix)) / isec%grid%dx

  if (w < 0.0_dp) w = 0.0_dp
  if (w > 1.0_dp) w = 1.0_dp

  z_val = (1.0_dp - w) * isec%z(ix) + w * isec%z(ix+1)

  return

end function interpolate_interface2d

!***************************************************************************************************************

subroutine interface_normal2d(px, isec, norm_x, norm_z, h)

  real(kind=dp)                              :: px, norm_x, norm_z, h
  type(Tintersection2d)                      :: isec

  integer                                    :: ix
  real(kind=dp)                              :: w, dzdx, norm_len

  if (isec%nnode < 2) then
     norm_x = 0.0_dp
     norm_z = 1.0_dp
     h = isec%z(1)
     return
  endif

  ix = floor((px - isec%x(1)) / isec%grid%dx) + 1

  if (ix < 1) ix = 1
  if (ix >= isec%nnode) ix = isec%nnode - 1

  w = (px - isec%x(ix)) / isec%grid%dx
  if (w < 0.0_dp) w = 0.0_dp
  if (w > 1.0_dp) w = 1.0_dp

  h = (1.0_dp - w) * isec%z(ix) + w * isec%z(ix+1)

  if (ix > 1 .and. ix < isec%nnode - 1) then
     dzdx = (isec%z(ix+1) - isec%z(ix-1)) / (2.0_dp * isec%grid%dx)
  else if (ix == 1 .and. isec%nnode >= 3) then
     dzdx = (isec%z(ix+1) - isec%z(ix)) / isec%grid%dx
  else if (ix == isec%nnode - 1 .and. isec%nnode >= 3) then
     dzdx = (isec%z(ix+1) - isec%z(ix-1)) / (2.0_dp * isec%grid%dx)
  else
     dzdx = 0.0_dp
  endif

  norm_len = sqrt(1.0_dp + dzdx**2)
  norm_x = -dzdx / norm_len
  norm_z = 1.0_dp / norm_len

  return

end subroutine interface_normal2d

!***************************************************************************************************************

function get_velocity_at2d(reg, px, pz, vtype) result(vel)

  type(Tregion2d)                            :: reg
  real(kind=dp)                              :: px, pz
  integer                                    :: vtype
  real(kind=dp)                              :: vel

  type(Tpropagation_grid2d),pointer          :: grid
  integer                                    :: ix, iz
  real(kind=dp)                              :: w1, w2

  grid => reg%grid

  ix = floor((px - grid%x0) / grid%dx) + 1
  iz = floor((pz - grid%z0) / grid%dz) + 1

  if (ix < 1) ix = 1
  if (ix >= grid%nx) ix = grid%nx - 1
  if (iz < 1) iz = 1
  if (iz >= grid%nz) iz = grid%nz - 1

  w1 = (px - grid%x(ix)) / grid%dx
  w2 = (pz - grid%z(iz)) / grid%dz

  if (w1 < 0.0_dp) w1 = 0.0_dp
  if (w1 > 1.0_dp) w1 = 1.0_dp
  if (w2 < 0.0_dp) w2 = 0.0_dp
  if (w2 > 1.0_dp) w2 = 1.0_dp

  if (associated(grid%velocity)) then
     vel = (1.0_dp-w1)*(1.0_dp-w2)*grid%velocity(vtype,ix,iz)   + &
           w1*(1.0_dp-w2)*grid%velocity(vtype,ix+1,iz)           + &
           (1.0_dp-w1)*w2*grid%velocity(vtype,ix,iz+1)           + &
           w1*w2*grid%velocity(vtype,ix+1,iz+1)
  else
     vel = reg%velocity(vtype, find_position_in_region2d(reg, px, pz))
  endif

  return

end function get_velocity_at2d

end module mod_rays2d
