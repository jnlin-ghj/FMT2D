module mod_frechet2d

use mod_fmt2d
use mod_rays2d

implicit none

private
public :: initialize_inversion2d, ray_partials2d, tag_active_vgrid_nodes2d, &
          decode_global_index2d

contains

subroutine initialize_inversion2d(do_frechet_derivatives)
  implicit none

  logical :: do_frechet_derivatives

  integer :: ngrids, nifaces, i, igrid, iiface, group, last_index, vtype, iz, ix
  integer :: n_active
  integer :: input_list(1000), type_list(100)
  type(Tvelocity_grid2d), pointer :: vg

  n_inv_parms = 0
  n_inv_active = 0
  group = 0

  open(1, file='frechet.in', status='old', action='read')

  read(1, *, end=100) i
  if (i == 0) then
     do_frechet_derivatives = .false.
     close(1)
     return
  else
     do_frechet_derivatives = .true.
  endif

  read(1, *, end=100) ngrids
  n_inv_vgrid = ngrids

  if (n_inv_vgrid > 0) then

     read(1, *, end=100) input_list(1:ngrids)
     read(1, *, end=100) type_list(1:ngrids)

     allocate(vgrids_to_be_inv(ngrids))

     do i = 1, ngrids
        igrid = input_list(i)
        vtype = type_list(i)
        vgrids_to_be_inv(i)%igrid = igrid
        vgrids_to_be_inv(i)%vtype = vtype
        vgrid(igrid, vtype)%to_be_inverted = .true.
        allocate(vgrid(igrid, vtype)%inv_idx(vgrid(igrid, vtype)%nx, vgrid(igrid, vtype)%nz))
        vgrid(igrid, vtype)%inv_idx = 0
        do iz = 1, vgrid(igrid, vtype)%nz
           do ix = 1, vgrid(igrid, vtype)%nx
              if (vgrid(igrid, vtype)%active(vtype, ix, iz)) then
                 n_inv_parms = n_inv_parms + 1
                 vgrid(igrid, vtype)%inv_idx(ix, iz) = n_inv_parms
              end if
           end do
        end do
        vgrid(igrid, vtype)%start_index = 0
        n_inv_active = n_inv_parms
     end do

     group = n_inv_vgrid
     last_index = n_inv_parms

  endif

  read(1, *, end=100) nifaces
  n_inv_iface = nifaces

  if (n_inv_iface > 0) then

     read(1, *, end=100) input_list(1:nifaces)

     allocate(ifaces_to_be_inv(nifaces))

     do i = 1, nifaces
        iiface = input_list(i)
        ifaces_to_be_inv(i) = iiface
        intrface(iiface)%to_be_inverted = .true.
        group = group + 1
        if (group == 1) then
           intrface(iiface)%start_index = 0
        else
           intrface(iiface)%start_index = last_index
        endif
        last_index = intrface(iiface)%start_index + intrface(iiface)%nnode
        n_inv_parms = n_inv_parms + intrface(iiface)%nnode
        n_inv_active = n_inv_active + intrface(iiface)%nnode
     end do

  endif

  close(1)

  write(*,'(A,I8)') '  Inversion parameters: n_inv_parms = ', n_inv_parms
  write(*,'(A,I8)') '  Inversion active:     n_inv_active = ', n_inv_active
  write(*,'(A,I4)') '  Velocity grids to invert: ', n_inv_vgrid
  write(*,'(A,I4)') '  Interfaces to invert:     ', n_inv_iface

  do i = 1, n_inv_vgrid
    igrid = vgrids_to_be_inv(i)%igrid
    vtype = vgrids_to_be_inv(i)%vtype
    vg => vgrid(igrid, vtype)
    n_active = 0
    do iz = 1, vg%nz
      do ix = 1, vg%nx
        if (vg%active(vtype, ix, iz)) n_active = n_active + 1
      end do
    end do
    write(*,'(A,I2,A,I6,A,I6,A,I8,A,F8.3,A,F8.3)') '  vgrid(', igrid, '): nx=', vg%nx, ' nz=', vg%nz, &
      ' n_active=', n_active, ' z(1)=', vg%z(1), ' z(nz)=', vg%z(vg%nz)
  end do

  return

100 continue
  write(*,'(A)') 'Error reading frechet.in'
  close(1)
  stop

end subroutine initialize_inversion2d


subroutine ray_partials2d(ray)
  implicit none

  type(Tray2d) :: ray

  type(Tray_section2d), pointer :: raysec, raysec_prev
  type(Tvelocity_grid2d), pointer :: gridv
  type(Tintersection2d), pointer :: isec
  type(Tinterface2d), pointer :: iface

  real(kind=dp), allocatable :: dtdpar(:)
  integer :: start_index, glob_inv_index
  integer :: i, n, m, ix, iz, ivx, ivz, vtype
  real(kind=dp) :: px, pz, dl, vel, weight
  real(kind=dp) :: wx, wz
  real(kind=dp) :: gradt_in(2), gradt_out(2), normal(2), h
  real(kind=dp) :: gradt_in_perp, gradt_out_perp, gradt_in_par
  real(kind=dp) :: vel_in, vel_out, geo_factor, det
  real(kind=dp) :: norm
  logical :: rec_region_done, do_interface_partials

  if (n_inv_parms <= 0) return

  allocate(dtdpar(n_inv_parms))
  dtdpar = 0.0_dp

  nullify(raysec_prev)
  rec_region_done = .false.
  do_interface_partials = .false.

  do m = ray%nsections, 1, -1

     raysec => ray%section(m)

     if (raysec%npoints <= 2) then
        cycle
     else
        if (.not. rec_region_done) then
           raysec_prev => raysec
           rec_region_done = .true.
        else
           do_interface_partials = .true.
        endif
     endif

     ! Skip water-layer sections — water velocity is fixed (not inverted)
     if (raysec%reg%id == 1) cycle

     vtype = raysec%tf%vtype
     if (vtype == 0) then
        print *, 'illegal vtype in ray_partials2d'
        stop
     endif

     gridv => vgrid(raysec%reg%ivgrid, vtype)

     if (gridv%to_be_inverted) then

        do n = 2, raysec%npoints

           px = 0.5_dp * (raysec%point(1,n) + raysec%point(1,n-1))
           pz = 0.5_dp * (raysec%point(2,n) + raysec%point(2,n-1))

           dl = sqrt((raysec%point(1,n) - raysec%point(1,n-1))**2 + &
                     (raysec%point(2,n) - raysec%point(2,n-1))**2)

           ivx = floor((px - gridv%x(1)) / gridv%dx) + 1

           if (pz < gridv%z(1)) then
              cycle
           else if (pz >= gridv%z(gridv%nz)) then
              ivz = gridv%nz - 1
              wz = 1.0_dp
           else
              ivz = floor((pz - gridv%z(1)) / gridv%dz) + 1
              if (ivz < 1) ivz = 1
              if (ivz >= gridv%nz) ivz = gridv%nz - 1
              wz = (pz - gridv%z(ivz)) / gridv%dz
           endif

           if (ivx < 1) ivx = 1
           if (ivx >= gridv%nx) ivx = gridv%nx - 1

           wx = (px - gridv%x(ivx)) / gridv%dx
           if (wx < 0.0_dp) wx = 0.0_dp
           if (wx > 1.0_dp) wx = 1.0_dp
           if (wz < 0.0_dp) wz = 0.0_dp
           if (wz > 1.0_dp) wz = 1.0_dp

           vel = (1.0_dp-wx)*(1.0_dp-wz)*gridv%velocity(vtype,ivx,ivz) + &
                 wx*(1.0_dp-wz)*gridv%velocity(vtype,ivx+1,ivz) + &
                 (1.0_dp-wx)*wz*gridv%velocity(vtype,ivx,ivz+1) + &
                 wx*wz*gridv%velocity(vtype,ivx+1,ivz+1)

           if (vel < 1.0e-6_dp) vel = 1.0e-6_dp

           glob_inv_index = gridv%inv_idx(ivx, ivz)
           if (glob_inv_index > 0) then
             dtdpar(glob_inv_index) = dtdpar(glob_inv_index) &
               + dl*(1.0_dp-wx)*(1.0_dp-wz)
           end if

           glob_inv_index = gridv%inv_idx(ivx+1, ivz)
           if (glob_inv_index > 0) then
             dtdpar(glob_inv_index) = dtdpar(glob_inv_index) &
               + dl*wx*(1.0_dp-wz)
           end if

           glob_inv_index = gridv%inv_idx(ivx, ivz+1)
           if (glob_inv_index > 0) then
             dtdpar(glob_inv_index) = dtdpar(glob_inv_index) &
               + dl*(1.0_dp-wx)*wz
           end if

           glob_inv_index = gridv%inv_idx(ivx+1, ivz+1)
           if (glob_inv_index > 0) then
             dtdpar(glob_inv_index) = dtdpar(glob_inv_index) &
               + dl*wx*wz
           end if

        end do

     endif

     if (do_interface_partials) then

        isec => raysec%iend
        iface => intrface(isec%iface_id)

        if (iface%to_be_inverted) then

           px = raysec%point(1, raysec%npoints)
           pz = raysec%point(2, raysec%npoints)

           call interpolate_time_gradient2d(raysec%reg, px, pz, gradt_in(1), gradt_in(2))

           vel_in = get_velocity_at2d(raysec%reg, px, pz, vtype)
           norm = sqrt(sum(gradt_in**2))
           if (norm > 0.0_dp) then
              gradt_in = gradt_in / (vel_in * norm)
           endif

           call interface_normal2d(px, isec, normal(1), normal(2), h)

           gradt_in_perp = dot_product(gradt_in, normal)

           if (raysec%reg%id == raysec_prev%reg%id) then
              gradt_out_perp = -gradt_in_perp
           else
              gradt_in_par = sqrt(sum((gradt_in - gradt_in_perp*normal)**2))
              vel_out = get_velocity_at2d(raysec_prev%reg, px, pz, raysec_prev%tf%vtype)
              det = 1.0_dp/vel_out**2 - gradt_in_par**2
              if (det < 0.0_dp) then
                 det = 0.0_dp
              endif
              gradt_out_perp = sign(sqrt(det), gradt_in_perp)
           endif

           geo_factor = (gradt_in_perp - gradt_out_perp) * normal(2)

           ivx = floor((px - iface%x(1)) / iface%dx) + 1
           if (ivx < 1) ivx = 1
           if (ivx >= iface%nx) ivx = iface%nx - 1

           wx = (px - iface%x(ivx)) / iface%dx
           if (wx < 0.0_dp) wx = 0.0_dp
           if (wx > 1.0_dp) wx = 1.0_dp

           start_index = iface%start_index

           glob_inv_index = start_index + ivx - 1
           dtdpar(glob_inv_index + 1) = dtdpar(glob_inv_index + 1) + (1.0_dp - wx)*geo_factor

           glob_inv_index = start_index + ivx
           dtdpar(glob_inv_index + 1) = dtdpar(glob_inv_index + 1) + wx*geo_factor

        endif

     endif

     raysec_prev => raysec

  end do

  ray%n_pdev = count(dtdpar /= 0.0_dp)
  if (associated(ray%pdev_indx)) deallocate(ray%pdev_indx)
  if (associated(ray%pdev)) deallocate(ray%pdev)
  allocate(ray%pdev_indx(ray%n_pdev), ray%pdev(ray%n_pdev))
  m = 0
  do n = 1, n_inv_parms
     if (dtdpar(n) /= 0.0_dp) then
        m = m + 1
        ray%pdev(m) = dtdpar(n)
        ray%pdev_indx(m) = n
     endif
  end do

  deallocate(dtdpar)

  return

end subroutine ray_partials2d


subroutine tag_active_vgrid_nodes2d(reg)
  implicit none

  type(Tregion2d) :: reg
  type(Tvelocity_grid2d), pointer :: gridv

  integer :: n, vtype, ivx, ivz

  do vtype = 1, NVTYPE

     gridv => vgrid(reg%ivgrid, vtype)

     if (.not. associated(gridv%active)) cycle

     do n = 1, reg%nnode

        if (reg%node(n)%i1 == 0) cycle

        ivx = floor((reg%x(n) - gridv%x(1)) / gridv%dx) + 1

        ivz = 1
        if (reg%z(n) >= gridv%z(1)) then
           do ivz = 1, gridv%nz - 1
              if (reg%z(n) >= gridv%z(ivz) .and. reg%z(n) < gridv%z(ivz+1)) exit
           end do
           if (ivz >= gridv%nz) ivz = gridv%nz - 1
        endif

        if (ivx >= 1 .and. ivx <= gridv%nx .and. ivz >= 1 .and. ivz <= gridv%nz) then
           if (associated(gridv%active)) then
              gridv%active(vtype, max(1,ivx-1):min(gridv%nx,ivx+1), &
                                max(1,ivz-1):min(gridv%nz,ivz+1)) = .true.
           endif
        endif

     end do

  end do

end subroutine tag_active_vgrid_nodes2d


subroutine decode_global_index2d(glob_index, i1, i2, p_type, p_index, p_subindex)
  implicit none

  integer :: glob_index, i1, i2, p_type, p_index, p_subindex
  integer :: i, n, m, group, loc_index
  integer :: start(31), n_groups

  p_subindex = 0

  n_groups = n_inv_vgrid + n_inv_iface
  if (n_groups > 30) stop 'too many groups in decode_global_index2d'

  m = 0
  do i = 1, n_inv_vgrid
     m = m + 1
     start(m) = vgrid(vgrids_to_be_inv(i)%igrid, vgrids_to_be_inv(i)%vtype)%start_index
  end do

  do i = 1, n_inv_iface
     m = m + 1
     start(m) = intrface(ifaces_to_be_inv(i))%start_index
  end do

  start(n_groups+1) = n_inv_parms + 1

  group = 0
  do n = 1, n_groups
     if (glob_index < start(n+1)) then
        group = n
        exit
     endif
  end do

  if (group == 0) then
     p_type = 0
     p_index = 0
     return
  endif

  if (group <= n_inv_vgrid) then

     p_type = 1
     p_index = vgrids_to_be_inv(group)%igrid
     p_subindex = vgrids_to_be_inv(group)%vtype

     loc_index = glob_index - start(group)
     i2 = loc_index / vgrid(p_index, p_subindex)%nx + 1
     i1 = loc_index - (i2-1)*vgrid(p_index, p_subindex)%nx + 1

     return

  endif

  if (group <= n_inv_vgrid + n_inv_iface) then

     p_type = 2
     p_index = ifaces_to_be_inv(group - n_inv_vgrid)

     loc_index = glob_index - start(group)
     i1 = loc_index + 1
     i2 = 0

     return

  endif

end subroutine decode_global_index2d


end module mod_frechet2d