program fmt2d

use type_definitions2d
use global_variables2d
use mod_rays2d
use read_tomo2d
use mod_fmt2d

implicit none

integer :: n, m, i, j, k, ns, ipath, itf, step, iz, ix
integer :: ntf, vtype, prev_tf
integer :: t1, t2, t3, t4, count_rate, count_max
integer :: path_id, nrec_total, iobs, idx
integer :: nobs_geom, nx_smesh, nz_smesh, nx_refl, nx_sed
integer :: ntf_max, n_paths_per_source
integer :: n_valid_pg, n_valid_psedp, n_valid_pmp, n_valid_pgss, n_valid_pgss2
integer :: nn, n_valid_top, n_alive_sweep, n_valid_fmm
integer :: n_updated, mm, icell, ixx, izz, i1, i2, i3, ii, jj, m_top, m_below

real(kind=dp) :: vwater, vair, z_seafloor, z_sedbot, z_moho, t_arrival
real(kind=dp) :: dtdx, dtdz, dtdz_sq, dx_if, v_loc, dtdz_orig
real(kind=dp) :: avg_sed_thick, avg_moho_depth
real(kind=dp) :: psedp_max_offset, pmp_min_offset, pmp_max_offset
real(kind=dp) :: src_x, rec_offset
real(kind=dp) :: dist_direct, t_direct_water, t_pg_fmm, t_pg_total
type(Ttime_field2d), pointer :: tf_ptr

real(kind=dp), allocatable :: xpos(:), topo(:), zpos_rel(:), vel(:,:)
real(kind=dp), allocatable :: refl_x(:), refl_z(:)
real(kind=dp), allocatable :: sed_x(:), sed_z(:)
real(kind=dp), allocatable :: obs_x(:), obs_z(:)
real(kind=dp), allocatable :: rec_x(:), rec_z(:)
integer, allocatable :: obs_id(:), rec_phase(:)
real(kind=dp), allocatable :: geom_obs_time(:), geom_uncertainty(:)

real(kind=dp), allocatable :: water_corr(:)
real(kind=dp), allocatable :: entry_x_arr(:,:), entry_z_arr(:,:)

type(Tsource2d), pointer :: s
type(Treceiver2d), pointer :: rec
type(Tray2d), pointer :: ray
type(Tregion2d), pointer :: reg
type(Tintersection2d), pointer :: isec

character(len=256) :: smesh_vp_file, refl_file, sed_file, geom_file
character(len=256) :: arrivals_file, obs_file, rays_prefix
logical :: write_rays

call system_clock(t1, count_rate, count_max)

write(*,'(A)') '============================================'
write(*,'(A)') '  FMT2D - Fast Marching Tomography 2D'
write(*,'(A)') '  Pure P-wave forward modeling'
write(*,'(A)') '============================================'
write(*,*)

! Read input file (8 lines)
open(1, file='fmt2d.in', status='old', action='read')
read(1, '(A)') smesh_vp_file
read(1, '(A)') refl_file
read(1, '(A)') sed_file
read(1, '(A)') geom_file
read(1, *) write_rays
read(1, '(A)') arrivals_file
read(1, '(A)') obs_file
read(1, '(A)') rays_prefix
close(1)

write(*,'(A,A)') '  smesh_vp:   ', trim(smesh_vp_file)
write(*,'(A,A)') '  refl_file:  ', trim(refl_file)
write(*,'(A,A)') '  sed_file:   ', trim(sed_file)
write(*,'(A,A)') '  geom_file:  ', trim(geom_file)
write(*,'(A,L)') '  write_rays: ', write_rays
write(*,'(A,A)') '  arrivals:   ', trim(arrivals_file)
write(*,'(A,A)') '  obs_file:   ', trim(obs_file)
write(*,'(A,A)') '  rays_prefix:', trim(rays_prefix)
write(*,*)

! Read velocity mesh, reflection interface, sediment interface, geometry
call read_smesh(smesh_vp_file, nx_smesh, nz_smesh, vwater, xpos, topo, zpos_rel, vel, vair)
vwater_global = vwater
vair_global = vair
call read_refl(refl_file, nx_refl, refl_x, refl_z)
call read_refl(sed_file, nx_sed, sed_x, sed_z)
call read_geom(geom_file, nobs_geom, obs_x, obs_z, nrec_total, rec_x, rec_z, obs_id, rec_phase, geom_obs_time, geom_uncertainty)

write(*,*)
write(*,'(A)') '  Building FMT2D model...'

! ============================================================
! 5 interfaces: top, seafloor, sed_bot, Moho, bottom
! ============================================================
n_interfaces = 5
allocate(intrface(n_interfaces))

do n = 1, n_interfaces
   call interface2d_defaults(intrface(n))
   intrface(n)%id = n
end do

! Interface 1: top (z=0)
intrface(1)%nx = nx_smesh
intrface(1)%nnode = nx_smesh
intrface(1)%dx = xpos(2) - xpos(1)
intrface(1)%x0 = xpos(1)
allocate(intrface(1)%x(nx_smesh), intrface(1)%z_nodes(nx_smesh))
intrface(1)%x = xpos
intrface(1)%z_nodes = 0.0_dp

! Interface 2: seafloor
intrface(2)%nx = nx_smesh
intrface(2)%nnode = nx_smesh
intrface(2)%dx = xpos(2) - xpos(1)
intrface(2)%x0 = xpos(1)
allocate(intrface(2)%x(nx_smesh), intrface(2)%z_nodes(nx_smesh))
intrface(2)%x = xpos
intrface(2)%z_nodes = topo

! Interface 3: sediment bottom
intrface(3)%nx = nx_sed
intrface(3)%nnode = nx_sed
intrface(3)%dx = sed_x(2) - sed_x(1)
intrface(3)%x0 = sed_x(1)
allocate(intrface(3)%x(nx_sed), intrface(3)%z_nodes(nx_sed))
intrface(3)%x = sed_x
intrface(3)%z_nodes = sed_z

! Interface 4: Moho
intrface(4)%nx = nx_refl
intrface(4)%nnode = nx_refl
intrface(4)%dx = refl_x(2) - refl_x(1)
intrface(4)%x0 = refl_x(1)
allocate(intrface(4)%x(nx_refl), intrface(4)%z_nodes(nx_refl))
intrface(4)%x = refl_x
intrface(4)%z_nodes = refl_z

! Interface 5: bottom
intrface(5)%nx = nx_smesh
intrface(5)%nnode = nx_smesh
intrface(5)%dx = xpos(2) - xpos(1)
intrface(5)%x0 = xpos(1)
allocate(intrface(5)%x(nx_smesh), intrface(5)%z_nodes(nx_smesh))
intrface(5)%x = xpos
intrface(5)%z_nodes = zpos_rel(nz_smesh) + maxval(topo)

write(*,'(A)') '  Interfaces set up'
write(*,'(A,F8.2,A,F8.2)') '    Iface 1 (top):     z = ', minval(intrface(1)%z_nodes), ' to ', maxval(intrface(1)%z_nodes)
write(*,'(A,F8.2,A,F8.2)') '    Iface 2 (seafloor):z = ', minval(intrface(2)%z_nodes), ' to ', maxval(intrface(2)%z_nodes)
write(*,'(A,F8.2,A,F8.2)') '    Iface 3 (sed_bot): z = ', minval(intrface(3)%z_nodes), ' to ', maxval(intrface(3)%z_nodes)
write(*,'(A,F8.2,A,F8.2)') '    Iface 4 (Moho):    z = ', minval(intrface(4)%z_nodes), ' to ', maxval(intrface(4)%z_nodes)
write(*,'(A,F8.2,A,F8.2)') '    Iface 5 (bottom):  z = ', minval(intrface(5)%z_nodes), ' to ', maxval(intrface(5)%z_nodes)

! ============================================================
! Velocity grid (P-wave only)
! ============================================================
n_vgrids = 1
n_vtypes = 1
allocate(vgrid(n_vgrids, NVTYPE))

do n = 1, n_vgrids
   do m = 1, NVTYPE
      call vgrid2d_defaults(vgrid(n, m))
   end do
end do

vgrid(1, VTYPE_P)%nx = nx_smesh
vgrid(1, VTYPE_P)%nz = nz_smesh
vgrid(1, VTYPE_P)%dx = xpos(2) - xpos(1)
vgrid(1, VTYPE_P)%dz = zpos_rel(2) - zpos_rel(1)
vgrid(1, VTYPE_P)%x0 = xpos(1)
vgrid(1, VTYPE_P)%z0 = zpos_rel(1)
vgrid(1, VTYPE_P)%nnode = nx_smesh * nz_smesh

allocate(vgrid(1, VTYPE_P)%x(nx_smesh))
allocate(vgrid(1, VTYPE_P)%z(nz_smesh))
allocate(vgrid(1, VTYPE_P)%velocity(NVTYPE, nx_smesh, nz_smesh))
allocate(vgrid(1, VTYPE_P)%velocity_init(NVTYPE, nx_smesh, nz_smesh))
allocate(vgrid(1, VTYPE_P)%active(NVTYPE, nx_smesh, nz_smesh))

vgrid(1, VTYPE_P)%x = xpos
vgrid(1, VTYPE_P)%z = zpos_rel
vgrid(1, VTYPE_P)%active = .false.

do n = 1, nx_smesh
   do m = 1, nz_smesh
      vgrid(1, VTYPE_P)%velocity(VTYPE_P, n, m) = vel(n, m)
      vgrid(1, VTYPE_P)%velocity_init(VTYPE_P, n, m) = vel(n, m)
   end do
end do

write(*,'(A)') '  Velocity grid set up'

! ============================================================
! Propagation grid (resample velocity to uniform dz)
! ============================================================
allocate(pgrid)
call pgrid2d_defaults(pgrid)

pgrid%nx = nx_smesh
pgrid%dz = xpos(2) - xpos(1)
pgrid%dx = pgrid%dz
pgrid%z0 = 0.0_dp
pgrid%zmax = zpos_rel(nz_smesh) + maxval(topo) + 5.0_dp
pgrid%nz = nint((pgrid%zmax - pgrid%z0) / pgrid%dz) + 1
pgrid%x0 = xpos(1)
pgrid%xmax = xpos(nx_smesh)
pgrid%tolerance = min(pgrid%dx, pgrid%dz) * 0.01_dp
pgrid%is_main_grid = .true.
pgrid%index_x0 = 1

allocate(pgrid%x(pgrid%nx), pgrid%z(pgrid%nz))
do n = 1, pgrid%nx
   pgrid%x(n) = pgrid%x0 + (n-1)*pgrid%dx
end do
do n = 1, pgrid%nz
   pgrid%z(n) = pgrid%z0 + (n-1)*pgrid%dz
end do

write(*,'(A,I6,A,I6)') '  Propagation grid: nx=', pgrid%nx, ' nz=', pgrid%nz
write(*,'(A,F8.2,A,F8.2)') '  x range: ', pgrid%x(1), ' to ', pgrid%x(pgrid%nx)
write(*,'(A,F8.2,A,F8.2)') '  z range: ', pgrid%z(1), ' to ', pgrid%z(pgrid%nz)

! Resample velocity to propagation grid
block
  integer :: nz_new, iz_new, iz_old, ix_new
  real(kind=dp) :: dz_new, z_val, wz_old
  real(kind=dp), allocatable :: zpos_new(:), vel_new(:,:)
  real(kind=dp) :: pz_rel_node, z_sf_node, wx_z_node
  integer :: ivx_z_node

  dz_new = pgrid%dz
  nz_new = pgrid%nz

  write(*,'(A,I6,A,F8.4)') '  Resampling: nz_new=', nz_new, ' dz_new=', dz_new

  allocate(zpos_new(nz_new), vel_new(nx_smesh, nz_new))
  do iz_new = 1, nz_new
    zpos_new(iz_new) = (iz_new - 1) * dz_new
  end do

  do ix_new = 1, nx_smesh
    ivx_z_node = ix_new
    if (ivx_z_node < 1) ivx_z_node = 1
    if (ivx_z_node >= nx_smesh) ivx_z_node = nx_smesh - 1
    wx_z_node = dble(ix_new - ivx_z_node)
    if (wx_z_node < 0.0_dp) wx_z_node = 0.0_dp
    if (wx_z_node > 1.0_dp) wx_z_node = 1.0_dp
    z_sf_node = (1.0_dp - wx_z_node) * topo(ivx_z_node) + wx_z_node * topo(ivx_z_node + 1)

    do iz_new = 1, nz_new
      z_val = zpos_new(iz_new)
      pz_rel_node = z_val - z_sf_node

      if (pz_rel_node < 0.0_dp) then
        vel_new(ix_new, iz_new) = vwater_global
      else if (pz_rel_node >= zpos_rel(nz_smesh)) then
        vel_new(ix_new, iz_new) = vel(ix_new, nz_smesh)
      else
        do iz_old = 1, nz_smesh - 1
          if (pz_rel_node >= zpos_rel(iz_old) .and. pz_rel_node < zpos_rel(iz_old+1)) exit
        end do
        if (iz_old >= nz_smesh) iz_old = nz_smesh - 1
        wz_old = (pz_rel_node - zpos_rel(iz_old)) / (zpos_rel(iz_old+1) - zpos_rel(iz_old))
        vel_new(ix_new, iz_new) = (1.0_dp - wz_old) * vel(ix_new, iz_old) + &
                                   wz_old * vel(ix_new, iz_old+1)
      end if
    end do
  end do

  ! Replace velocity grid with resampled version
  if (associated(vgrid(1, VTYPE_P)%x)) deallocate(vgrid(1, VTYPE_P)%x)
  if (associated(vgrid(1, VTYPE_P)%z)) deallocate(vgrid(1, VTYPE_P)%z)
  if (associated(vgrid(1, VTYPE_P)%velocity)) deallocate(vgrid(1, VTYPE_P)%velocity)
  if (associated(vgrid(1, VTYPE_P)%velocity_init)) deallocate(vgrid(1, VTYPE_P)%velocity_init)
  if (associated(vgrid(1, VTYPE_P)%active)) deallocate(vgrid(1, VTYPE_P)%active)

  nz_smesh = nz_new
  vgrid(1, VTYPE_P)%nz = nz_new
  vgrid(1, VTYPE_P)%dz = dz_new
  vgrid(1, VTYPE_P)%z0 = zpos_new(1)
  vgrid(1, VTYPE_P)%nnode = nx_smesh * nz_new

  allocate(vgrid(1, VTYPE_P)%x(nx_smesh))
  allocate(vgrid(1, VTYPE_P)%z(nz_new))
  allocate(vgrid(1, VTYPE_P)%velocity(NVTYPE, nx_smesh, nz_new))
  allocate(vgrid(1, VTYPE_P)%velocity_init(NVTYPE, nx_smesh, nz_new))
  allocate(vgrid(1, VTYPE_P)%active(NVTYPE, nx_smesh, nz_new))

  vgrid(1, VTYPE_P)%x = xpos
  vgrid(1, VTYPE_P)%z = zpos_new
  vgrid(1, VTYPE_P)%active = .false.

  do ix_new = 1, nx_smesh
    ivx_z_node = ix_new
    if (ivx_z_node < 1) ivx_z_node = 1
    if (ivx_z_node >= nx_smesh) ivx_z_node = nx_smesh - 1
    wx_z_node = dble(ix_new - ivx_z_node)
    if (wx_z_node < 0.0_dp) wx_z_node = 0.0_dp
    if (wx_z_node > 1.0_dp) wx_z_node = 1.0_dp
    z_sf_node = (1.0_dp - wx_z_node) * topo(ivx_z_node) + wx_z_node * topo(ivx_z_node + 1)

    do iz_new = 1, nz_new
      vgrid(1, VTYPE_P)%velocity(VTYPE_P, ix_new, iz_new) = vel_new(ix_new, iz_new)
      vgrid(1, VTYPE_P)%velocity_init(VTYPE_P, ix_new, iz_new) = vel_new(ix_new, iz_new)
      pz_rel_node = zpos_new(iz_new) - z_sf_node
      if (pz_rel_node >= 0.0_dp) then
        vgrid(1, VTYPE_P)%active(VTYPE_P, ix_new, iz_new) = .true.
      end if
    end do
  end do

  deallocate(zpos_rel, vel)
  allocate(zpos_rel(nz_new), vel(nx_smesh, nz_new))
  zpos_rel = zpos_new
  vel = vel_new

  deallocate(zpos_new, vel_new)

  write(*,'(A,I6,A,I6)') '  Velocity grid (resampled): nx=', nx_smesh, ' nz=', nz_smesh
end block

! ============================================================
! Intersections
! ============================================================
n_intersections = n_interfaces
allocate(intersection(n_intersections))
do n = 1, n_intersections
   call intersection2d_defaults(intersection(n))
   intersection(n)%id = n
end do

do n = 1, n_intersections
   call find_intersection2d(intersection(n), intrface(n), pgrid)
   write(*,'(A,I0,A,I8)') '  Intersection ', n, ' nnode=', intersection(n)%nnode
end do

! ============================================================
! Regions: water(1), sediment(2), crust(3), mantle(4), combined_all(5),
!          combined_sed_crust(6)
! Region 5 = sediment+crust+mantle combined (seafloor to bottom)
!            Used for Pg (diving wave) — global first arrival
! Region 6 = sediment+crust combined (seafloor to Moho)
!            Used for PmP (Moho reflection) up-going wave
! ============================================================
n_regions = 6
allocate(region(n_regions))

do n = 1, 4
   call region2d_defaults(region(n))
   region(n)%id = n
   region(n)%ivgrid = 1
   region(n)%itop => intersection(n)
   region(n)%ibot => intersection(n+1)

   if (region(n)%itop%nnode == 0 .or. region(n)%ibot%nnode == 0) then
      write(*,'(A,I0,A)') '  Region ', n, ' does not exist'
      region(n)%nnode = 0
      cycle
   endif

   call define_region2d(region(n), region(n)%itop, region(n)%ibot, pgrid)

   write(*,'(A,I0,A,I8)') '  Region ', n, ' nnode=', region(n)%nnode
end do

! Region 5: sediment+crust+mantle combined (seafloor to bottom)
if (associated(intersection(2)%rright_node_id)) deallocate(intersection(2)%rright_node_id)
if (associated(intersection(5)%rleft_node_id)) deallocate(intersection(5)%rleft_node_id)

call region2d_defaults(region(5))
region(5)%id = 5
region(5)%ivgrid = 1
region(5)%itop => intersection(2)  ! seafloor
region(5)%ibot => intersection(5)  ! bottom

call define_region2d(region(5), region(5)%itop, region(5)%ibot, pgrid)

write(*,'(A,I0,A,I8)') '  Region 5 (combined all) nnode=', region(5)%nnode

intersection(2)%regright => region(5)
intersection(5)%regleft => region(5)

! Region 6: sediment+crust combined (seafloor to Moho)
! Used for PmP Moho reflection up-going wave
call region2d_defaults(region(6))
region(6)%id = 6
region(6)%ivgrid = 1
region(6)%itop => intersection(2)  ! seafloor
region(6)%ibot => intersection(4)  ! Moho

call define_region2d(region(6), region(6)%itop, region(6)%ibot, pgrid)

write(*,'(A,I0,A,I8)') '  Region 6 (seafloor->Moho) nnode=', region(6)%nnode

write(*,'(A)') '  Calling tag_regular_nodes2d...'
call tag_regular_nodes2d(pgrid)
write(*,'(A)') '  Calling velocities_to_grid2d...'
call velocities_to_grid2d(pgrid)
write(*,'(A)') '  Calling velocities_to_intersection2d...'
do n = 1, n_intersections
   call velocities_to_intersection2d(intersection(n))
end do
write(*,'(A)') '  Calling velocities_to_region2d...'
do n = 1, n_regions
   if (region(n)%nnode > 0) then
      write(*,'(A,I0)') '    region ', n
      call velocities_to_region2d(region(n), pgrid)
   endif
end do

write(*,'(A)') '  Velocities assigned'

do n = 2, n_regions
   if (region(n)%nnode > 0) then
      call tag_active_vgrid_nodes2d(region(n))
   endif
end do

! Set remaining intersection region pointers
intersection(1)%regleft => null()
intersection(1)%regright => region(1)
intersection(3)%regleft => region(2)
intersection(3)%regright => region(3)
intersection(4)%regleft => region(3)
intersection(4)%regright => region(4)
intersection(5)%regright => null()

write(*,'(A)') '  Model setup complete'
write(*,*)

call system_clock(t2, count_rate, count_max)

! ============================================================
! Sources and receivers
! ============================================================
n_sources = nobs_geom
n_receivers = nrec_total

allocate(source(n_sources))
do n = 1, n_sources
   call source2d_defaults(source(n))
   source(n)%id = n
   source(n)%x = obs_x(n)
   source(n)%z = obs_z(n)
end do

allocate(receiver(n_receivers))
allocate(water_corr(n_receivers))
allocate(water_time(n_receivers))
allocate(rec_x_orig(n_receivers))
allocate(rec_z_orig(n_receivers))
allocate(in_water(n_receivers))
allocate(entry_x_best(n_receivers))
allocate(entry_z_best(n_receivers))
allocate(entry_x_arr(n_receivers, 10))
allocate(entry_z_arr(n_receivers, 10))
entry_x_best = 0.0_dp
entry_z_best = 0.0_dp
entry_x_arr = 0.0_dp
entry_z_arr = 0.0_dp

do n = 1, n_receivers
   call receiver2d_defaults(receiver(n))
   receiver(n)%id = n
   rec_x_orig(n) = rec_x(n)
   rec_z_orig(n) = rec_z(n)
   receiver(n)%x = rec_x(n)
   z_seafloor = interpolate_iface_z(intrface(2), rec_x(n))
   if (rec_z(n) < z_seafloor) then
      in_water(n) = .true.
      receiver(n)%z = z_seafloor
      water_corr(n) = (z_seafloor - rec_z(n)) / vwater
   else
      in_water(n) = .false.
      receiver(n)%z = rec_z(n)
      water_corr(n) = 0.0_dp
   endif
   water_time(n) = water_corr(n)
end do

write(*,'(A,I0)') '  Sources:   ', n_sources
write(*,'(A,I0)') '  Receivers: ', n_receivers

! ============================================================
! Compute offset thresholds for PsedP and PmP based on model geometry
! ============================================================
! Average sediment thickness below seafloor
avg_sed_thick = 0.0_dp
do n = 1, intrface(3)%nx
   z_seafloor = interpolate_iface_z(intrface(2), intrface(3)%x(n))
   avg_sed_thick = avg_sed_thick + (intrface(3)%z_nodes(n) - z_seafloor)
end do
avg_sed_thick = avg_sed_thick / dble(intrface(3)%nx)

! Average Moho depth below seafloor
avg_moho_depth = 0.0_dp
do n = 1, intrface(4)%nx
   z_seafloor = interpolate_iface_z(intrface(2), intrface(4)%x(n))
   avg_moho_depth = avg_moho_depth + (intrface(4)%z_nodes(n) - z_seafloor)
end do
avg_moho_depth = avg_moho_depth / dble(intrface(4)%nx)

! PsedP: max offset approx 6*sediment thickness (ray angle < ~80 deg from vertical)
psedp_max_offset = 6.0_dp * avg_sed_thick
! PmP: min offset = 3 km, max offset = 30 km
pmp_min_offset = 3.0_dp
pmp_max_offset = 30.0_dp

write(*,'(A,F8.2,A)') '  Avg sediment thickness: ', avg_sed_thick, ' km'
write(*,'(A,F8.2,A)') '  Avg Moho depth:         ', avg_moho_depth, ' km'
write(*,'(A,F8.2,A)') '  PsedP max offset:       ', psedp_max_offset, ' km'
write(*,'(A,F8.2,A,F8.2,A)') '  PmP offset range:       ', pmp_min_offset, ' - ', pmp_max_offset, ' km'
write(*,'(A)') '  Pg_ss offset range:     0 - 50 km'

! ============================================================
! Define 5 P-wave phases:
!   1 = Pg     (diving wave - global first arrival)
!   2 = PsedP  (sediment bottom reflection)
!   3 = PmP    (Moho reflection)
!   4 = Pg_ss  (Pg sea-surface multiple: solid up -> water up -> reflect -> water down)
!   5 = Pg_ss2 (water-first multiple: water up -> surface reflect -> water down -> solid diving -> water up)
!
! Time fields:
!   tf1:  Point-source FMM in region 5 (combined_all) - Pg
!   tf2:  Point-source FMM in region 2 (sediment) - PsedP down
!   tf3:  Reflect at sed_bot, sweep back through region 2 - PsedP up
!   tf4:  Point-source FMM in region 6 (seafloor->Moho) - PmP down
!   tf5:  Reflect at Moho, sweep P-wave back through region 6 - PmP up
!   tf6:  Transfer tf1 to seafloor, FMM in water up to sea surface  - Pg_ss
!   tf7:  Transfer tf6 to sea surface, reflect, FMM down in water   - Pg_ss
!   tf8:  Transfer tf7 to seafloor, reflect, FMM water up to receiver - Pg_ss
!   tf9:  Point-source FMM in water (region 1), sweep up from seafloor - Pg_ss2
!   tf10: Transfer tf9 to sea surface, reflect, FMM water down to seafloor - Pg_ss2
!   tf11: Transfer tf10 to seafloor, FMM in solid (region 5) - Pg_ss2
!   tf12: Transfer tf11 to seafloor, FMM water up to receiver - Pg_ss2
! ============================================================

n_paths_per_source = 5

do ns = 1, n_sources
   s => source(ns)

   ! Place source on seafloor
   z_seafloor = interpolate_iface_z(intrface(2), s%x)
   s%z = z_seafloor
   s%on_interface = .true.
   s%topint_id = 2
   s%botint_id = 2
   s%region_id = 5  ! source is in region 5 (seafloor->bottom) for Pg

   s%n_paths = n_paths_per_source
   allocate(s%path(n_paths_per_source))

   do n = 1, n_paths_per_source
      call path2d_defaults(s%path(n))
      s%path(n)%id = n
   end do

   ! Phase 1: Pg - diving wave (global first arrival)
   ! Point-source FMM in combined_all region (seafloor to bottom)
   s%path(1)%n_tf = 1
   allocate(s%path(1)%sequence(2))
   allocate(s%path(1)%tf_sequence(1))
   allocate(s%path(1)%vtype_sequence(1))
   s%path(1)%sequence(1) = 0   ! source
   s%path(1)%sequence(2) = 2   ! seafloor
   s%path(1)%vtype_sequence(1) = VTYPE_P
   s%path(1)%valid = .true.
   s%path(1)%used = .true.

   ! Phase 2: PsedP - sediment bottom reflection
   s%path(2)%n_tf = 2
   allocate(s%path(2)%sequence(4))
   allocate(s%path(2)%tf_sequence(2))
   allocate(s%path(2)%vtype_sequence(2))
   s%path(2)%sequence(1) = 0   ! source
   s%path(2)%sequence(2) = 3   ! sed_bot
   s%path(2)%sequence(3) = 3   ! sed_bot (reflected)
   s%path(2)%sequence(4) = 0   ! seafloor
   s%path(2)%vtype_sequence(1) = VTYPE_P
   s%path(2)%vtype_sequence(2) = VTYPE_P
   s%path(2)%valid = .true.
   s%path(2)%used = .true.

   ! Phase 3: PmP - Moho reflection
   ! tf4: source FMM in region 6 (seafloor->Moho) - down to Moho
   ! tf5: reflect at Moho, sweep back through region 6 - up to seafloor
   s%path(3)%n_tf = 2
   allocate(s%path(3)%sequence(4))
   allocate(s%path(3)%tf_sequence(2))
   allocate(s%path(3)%vtype_sequence(2))
   s%path(3)%sequence(1) = 0   ! source
   s%path(3)%sequence(2) = 4   ! Moho
   s%path(3)%sequence(3) = 4   ! Moho (reflected)
   s%path(3)%sequence(4) = 0   ! seafloor
   s%path(3)%vtype_sequence(1) = VTYPE_P
   s%path(3)%vtype_sequence(2) = VTYPE_P
   s%path(3)%valid = .true.
   s%path(3)%used = .true.

   ! Phase 4: Pg_ss — Pg sea-surface + seafloor double multiple
   ! Path: Source → Seafloor → Sea Surface → Seafloor → Sea Surface (receiver)
   ! tf6: transfer tf1 to seafloor, FMM water up to surface
   ! tf7: transfer tf6 to surface, reflect at sea surface, FMM water down to seafloor
   ! tf8: transfer tf7 to seafloor, reflect at seafloor, FMM water up to receiver
   s%path(4)%n_tf = 4
   allocate(s%path(4)%sequence(8))
   allocate(s%path(4)%tf_sequence(4))
   allocate(s%path(4)%vtype_sequence(4))
   s%path(4)%sequence(1) = 0   ! source
   s%path(4)%sequence(2) = 2   ! seafloor
   s%path(4)%sequence(3) = 2   ! seafloor
   s%path(4)%sequence(4) = 1   ! sea surface
   s%path(4)%sequence(5) = 1   ! sea surface
   s%path(4)%sequence(6) = 2   ! seafloor
   s%path(4)%sequence(7) = 2   ! seafloor
   s%path(4)%sequence(8) = 2   ! seafloor (istart for last section)
   s%path(4)%tf_sequence(1) = 1   ! solid path (tf1)
   s%path(4)%tf_sequence(2) = 6   ! water up 1st
   s%path(4)%tf_sequence(3) = 7   ! water down (sea surface reflection)
   s%path(4)%tf_sequence(4) = 8   ! water up 2nd (seafloor reflection)
   s%path(4)%vtype_sequence(1) = VTYPE_P
   s%path(4)%vtype_sequence(2) = VTYPE_P
   s%path(4)%vtype_sequence(3) = VTYPE_P
   s%path(4)%vtype_sequence(4) = VTYPE_P
   s%path(4)%valid = .true.
   s%path(4)%used = .true.

   ! Phase 5: Pg_ss2 — water-first multiple: source → water up → sea surface reflect → water down → solid diving → water up → receiver
   ! tf9:  point-source FMM in water (region 1), sweep up from seafloor
   ! tf10: transfer tf9 to sea surface, reflect, FMM water down to seafloor
   ! tf11: transfer tf10 to seafloor, FMM in solid (region 5)
   ! tf12: transfer tf11 to seafloor, FMM water up to receiver
   s%path(5)%n_tf = 4
   allocate(s%path(5)%sequence(8))
   allocate(s%path(5)%tf_sequence(4))
   allocate(s%path(5)%vtype_sequence(4))
   s%path(5)%sequence(1) = 0   ! source
   s%path(5)%sequence(2) = 1   ! sea surface (end of water up)
   s%path(5)%sequence(3) = 1   ! sea surface (start of water down)
   s%path(5)%sequence(4) = 2   ! seafloor (end of water down)
   s%path(5)%sequence(5) = 2   ! seafloor (start of solid)
   s%path(5)%sequence(6) = 2   ! seafloor (end of solid)
   s%path(5)%sequence(7) = 2   ! seafloor (start of water up)
   s%path(5)%sequence(8) = 2   ! seafloor (istart for last section)
   s%path(5)%tf_sequence(1) = 9   ! water up (source → sea surface)
   s%path(5)%tf_sequence(2) = 10  ! water down (sea surface → seafloor)
   s%path(5)%tf_sequence(3) = 11  ! solid diving (seafloor → seafloor)
   s%path(5)%tf_sequence(4) = 12  ! water up (seafloor → receiver)
   s%path(5)%vtype_sequence(1) = VTYPE_P
   s%path(5)%vtype_sequence(2) = VTYPE_P
   s%path(5)%vtype_sequence(3) = VTYPE_P
   s%path(5)%vtype_sequence(4) = VTYPE_P
   s%path(5)%valid = .true.
   s%path(5)%used = .true.

   write(*,'(A,I0,A,F8.2,A,F8.2,A,F8.2)') '  Source ', ns, ' at x=', s%x, ' z=', s%z, &
      ' seafloor=', z_seafloor

   ntf_max = 12
   allocate(s%time_field(ntf_max))
   do n = 1, ntf_max
      call time_field2d_defaults(s%time_field(n))
      s%time_field(n)%id = n
   end do

end do

! Assign receivers to phases based on source-receiver offset
do n = 1, n_receivers
   iobs = obs_id(n)
   if (rec_phase(n) >= 1 .and. rec_phase(n) <= 5) then
      ! User explicitly specified phase
      receiver(n)%n_rays = 1
      allocate(receiver(n)%ray(1))
      call ray2d_defaults(receiver(n)%ray(1))
      receiver(n)%ray(1)%source_id = iobs
      receiver(n)%ray(1)%raypath_id = rec_phase(n)
      receiver(n)%ray(1)%source => source(iobs)
   else
      ! Auto-assign phases based on source-receiver offset
      src_x = obs_x(iobs)
      rec_offset = abs(rec_x(n) - src_x)

      ! Count active phases using if-else logic
      ! Pg: always active (0-100 km)
      receiver(n)%n_rays = 1  ! Pg always

      ! PsedP: 3 km to psedp_max_offset
      if (rec_offset >= 3.0_dp .and. rec_offset <= psedp_max_offset) then
         receiver(n)%n_rays = receiver(n)%n_rays + 1
      end if

      ! PmP: pmp_min_offset to pmp_max_offset
      if (rec_offset >= pmp_min_offset .and. rec_offset <= pmp_max_offset) then
         receiver(n)%n_rays = receiver(n)%n_rays + 1
      end if

      ! Pg_ss: 0-50 km
      if (rec_offset <= 50.0_dp) then
         receiver(n)%n_rays = receiver(n)%n_rays + 1
      end if

      allocate(receiver(n)%ray(receiver(n)%n_rays))

      m = 1
      call ray2d_defaults(receiver(n)%ray(m))
      receiver(n)%ray(m)%source_id = iobs
      receiver(n)%ray(m)%raypath_id = 1   ! Pg
      receiver(n)%ray(m)%source => source(iobs)

      if (rec_offset >= 3.0_dp .and. rec_offset <= psedp_max_offset) then
         m = m + 1
         call ray2d_defaults(receiver(n)%ray(m))
         receiver(n)%ray(m)%source_id = iobs
         receiver(n)%ray(m)%raypath_id = 2   ! PsedP
         receiver(n)%ray(m)%source => source(iobs)
      end if

      if (rec_offset >= pmp_min_offset .and. rec_offset <= pmp_max_offset) then
         m = m + 1
         call ray2d_defaults(receiver(n)%ray(m))
         receiver(n)%ray(m)%source_id = iobs
         receiver(n)%ray(m)%raypath_id = 3   ! PmP
         receiver(n)%ray(m)%source => source(iobs)
      end if

      if (rec_offset <= 50.0_dp) then
         m = m + 1
         call ray2d_defaults(receiver(n)%ray(m))
         receiver(n)%ray(m)%source_id = iobs
         receiver(n)%ray(m)%raypath_id = 4   ! Pg_ss
         receiver(n)%ray(m)%source => source(iobs)
      end if
   endif
end do

write(*,'(A)') '  Sources and receivers initialized'
write(*,*)

call system_clock(t3, count_rate, count_max)

! ============================================================
! Forward modeling
! ============================================================
write(*,'(A)') '  Starting forward modeling...'
write(*,*)

global_source_counter = 0
n_valid_pg = 0
n_valid_psedp = 0
n_valid_pmp = 0
n_valid_pgss = 0
n_valid_pgss2 = 0

do ns = 1, n_sources

   s => source(ns)

   write(*,'(A,I0,A,F8.2,A,F8.2)') '  Source ', ns, ' at x=', s%x, ' z=', s%z

   ! Reset all regions
   do n = 1, n_regions
      reg => region(n)
      reg%loaded_tf => null()
      if (reg%nnode == 0) cycle
      reg%node_status = NODE_FAR
      reg%arrivaltime = huge_time
      reg%time_gradient = 0.0_dp
   end do

   ! ============================================================
   ! tf1: Pg - Point-source FMM in region 5 (seafloor->bottom)
   ! Global first-arrival: includes crust and mantle diving waves
   ! ============================================================
   reg => region(5)
   if (reg%nnode > 0) then
      call initialize_region_from_source2d(reg, s%x, s%z, 0.0_dp, VTYPE_P)
      call propagate2d(reg, VTYPE_P)

      s%n_time_fields = 1
      s%time_field(1)%reg => reg
      s%time_field(1)%vtype = VTYPE_P
      s%time_field(1)%istart => intersection(2)
      s%time_field(1)%inonstart => intersection(5)
      allocate(s%time_field(1)%arrivaltime(reg%nnode))
      allocate(s%time_field(1)%time_gradient(2, reg%nnode))
      s%time_field(1)%arrivaltime(:) = reg%arrivaltime(:)
      s%time_field(1)%time_gradient(:,:) = reg%time_gradient(:,:)

      s%path(1)%tf_sequence(1) = 1  ! Pg

      n_valid_fmm = count(reg%arrivaltime < huge_time*0.5)
      write(*,'(A,I8,A,I8)') '    tf1 (Pg) valid nodes = ', n_valid_fmm, ' / ', reg%nnode
   endif

   ! ============================================================
   ! tf2: PsedP down - Point-source FMM in region 2 (sediment)
   ! ============================================================
   reg => region(2)
   if (reg%nnode > 0) then
      call initialize_region_from_source2d(reg, s%x, s%z, 0.0_dp, VTYPE_P)
      call propagate2d(reg, VTYPE_P)

      s%n_time_fields = 2
      s%time_field(2)%reg => reg
      s%time_field(2)%vtype = VTYPE_P
      s%time_field(2)%istart => intersection(2)
      s%time_field(2)%inonstart => intersection(3)
      allocate(s%time_field(2)%arrivaltime(reg%nnode))
      allocate(s%time_field(2)%time_gradient(2, reg%nnode))
      s%time_field(2)%arrivaltime(:) = reg%arrivaltime(:)
      s%time_field(2)%time_gradient(:,:) = reg%time_gradient(:,:)

      s%path(2)%tf_sequence(1) = 2  ! PsedP: source -> sed_bot

      n_valid_fmm = count(reg%arrivaltime < huge_time*0.5)
      write(*,'(A,I8,A,I8)') '    tf2 (PsedP down) valid nodes = ', n_valid_fmm, ' / ', reg%nnode
   endif

   ! ============================================================
   ! tf3: PsedP up - reflect at sed_bot, sweep back through region 2
   ! ============================================================
   isec => intersection(3)  ! sed_bot
   if (isec%nnode > 0 .and. s%n_time_fields >= 2) then
      reg => region(2)
      reg%arrivaltime = s%time_field(2)%arrivaltime
      reg%time_gradient = s%time_field(2)%time_gradient

      call transfer_times_to_intersection(isec, s%time_field(2))

      reg => region(2)
      call reflect_gradient_at_interface(isec, reg, VTYPE_P)
      reg%node_status = NODE_FAR
      reg%arrivaltime = huge_time
      reg%time_gradient = 0.0_dp

      call sweep_region_from_interface2d(reg, isec, VTYPE_P, s)
      call propagate2d(reg, VTYPE_P)

      s%n_time_fields = 3
      s%time_field(3)%reg => reg
      s%time_field(3)%vtype = VTYPE_P
      s%time_field(3)%istart => intersection(3)
      s%time_field(3)%inonstart => intersection(2)
      allocate(s%time_field(3)%arrivaltime(reg%nnode))
      allocate(s%time_field(3)%time_gradient(2, reg%nnode))
      s%time_field(3)%arrivaltime(:) = reg%arrivaltime(:)
      s%time_field(3)%time_gradient(:,:) = reg%time_gradient(:,:)

      s%path(2)%tf_sequence(2) = 3  ! PsedP: sed_bot -> seafloor

      n_valid_fmm = count(reg%arrivaltime < huge_time*0.5)
      write(*,'(A,I8,A,I8)') '    tf3 (PsedP up) valid nodes = ', n_valid_fmm, ' / ', reg%nnode
   endif

   ! ============================================================
   ! tf4: PmP down - Point-source FMM in region 6 (seafloor->Moho)
   ! This ensures the down-going wavefield only includes paths above
   ! the Moho (no mantle diving waves), preventing PmP rays from
   ! crossing the Moho during ray tracing.
   ! ============================================================
   reg => region(6)
   if (reg%nnode > 0) then
      call initialize_region_from_source2d(reg, s%x, s%z, 0.0_dp, VTYPE_P)
      call propagate2d(reg, VTYPE_P)

      s%n_time_fields = 4
      s%time_field(4)%reg => reg
      s%time_field(4)%vtype = VTYPE_P
      s%time_field(4)%istart => intersection(2)
      s%time_field(4)%inonstart => intersection(4)
      allocate(s%time_field(4)%arrivaltime(reg%nnode))
      allocate(s%time_field(4)%time_gradient(2, reg%nnode))
      s%time_field(4)%arrivaltime(:) = reg%arrivaltime(:)
      s%time_field(4)%time_gradient(:,:) = reg%time_gradient(:,:)

      s%path(3)%tf_sequence(1) = 4  ! PmP: source -> Moho

      n_valid_fmm = count(reg%arrivaltime < huge_time*0.5)
      write(*,'(A,I8,A,I8)') '    tf4 (PmP down) valid nodes = ', n_valid_fmm, ' / ', reg%nnode
   endif

   ! ============================================================
   ! tf5: PmP up - extract Moho times from tf4 (region 6), reflect, sweep
   ! back through region 6 (seafloor->Moho). Since tf4 was computed in
   ! region 6 (no mantle paths), the reflected rays cannot cross the Moho.
   ! ============================================================
   isec => intersection(4)  ! Moho
   if (isec%nnode > 0 .and. s%n_time_fields >= 4) then
      ! Allocate intersection times if needed
      if (.not. associated(isec%arrivaltime)) then
         allocate(isec%arrivaltime(isec%nnode))
         allocate(isec%time_gradient(2, isec%nnode))
      endif

      ! Transfer Moho arrival times from tf4 (region 6, above Moho only)
      call transfer_times_to_intersection(isec, s%time_field(4))

      ! Sweep reflected wave from Moho through region 6 (seafloor->Moho only)
      reg => region(6)
      call reflect_gradient_at_interface(isec, reg, VTYPE_P)
      reg%node_status = NODE_FAR
      reg%arrivaltime = huge_time
      reg%time_gradient = 0.0_dp

      call sweep_region_from_interface2d(reg, isec, VTYPE_P, s)

      ! Propagate only within region 6 (above Moho)
      call propagate2d(reg, VTYPE_P)

      s%n_time_fields = 5
      s%time_field(5)%reg => reg
      s%time_field(5)%vtype = VTYPE_P
      s%time_field(5)%istart => intersection(4)
      s%time_field(5)%inonstart => intersection(2)
      allocate(s%time_field(5)%arrivaltime(reg%nnode))
      allocate(s%time_field(5)%time_gradient(2, reg%nnode))
      s%time_field(5)%arrivaltime(:) = reg%arrivaltime(:)
      s%time_field(5)%time_gradient(:,:) = reg%time_gradient(:,:)

      s%path(3)%tf_sequence(2) = 5  ! PmP: Moho -> seafloor

      n_valid_fmm = count(reg%arrivaltime < huge_time*0.5)
      write(*,'(A,I8,A,I8)') '    tf5 (PmP up) valid nodes = ', n_valid_fmm, ' / ', reg%nnode
   endif

   ! ============================================================
   ! tf6: Pg_ss water-up - Transfer tf1 to seafloor, FMM in water up to surface
   ! ============================================================
   isec => intersection(2)  ! seafloor
   if (isec%nnode > 0 .and. s%n_time_fields >= 1) then
      call transfer_times_to_intersection(isec, s%time_field(1))

      reg => region(1)
      if (associated(reg)) then
         call reflect_gradient_at_interface(isec, reg, VTYPE_P)
         reg%node_status = NODE_FAR
         reg%arrivaltime = huge_time
         reg%time_gradient = 0.0_dp
         call sweep_region_from_interface2d(reg, isec, VTYPE_P, s)
         call propagate2d(reg, VTYPE_P)

         ! Fix top boundary gradient: dtdz must be < 0 for upward FMM
         do ix = 1, reg%grid%nx
            m_top = reg%region_rnode_id(ix, 1)
            m_below = reg%region_rnode_id(ix, 2)
            if (m_top > 0 .and. m_below > 0) then
               if (reg%time_gradient(2, m_top) > 0.0_dp .or. &
                   abs(reg%time_gradient(2, m_top)) < 0.01_dp) then
                  reg%time_gradient(1:2, m_top) = reg%time_gradient(1:2, m_below)
               endif
            endif
         end do

         s%n_time_fields = 6
         s%time_field(6)%reg => reg
         s%time_field(6)%vtype = VTYPE_P
         s%time_field(6)%istart => intersection(2)
         s%time_field(6)%inonstart => intersection(1)
         allocate(s%time_field(6)%arrivaltime(reg%nnode))
         allocate(s%time_field(6)%time_gradient(2, reg%nnode))
         s%time_field(6)%arrivaltime(:) = reg%arrivaltime(:)
         s%time_field(6)%time_gradient(:,:) = reg%time_gradient(:,:)

         n_valid_fmm = count(reg%arrivaltime < huge_time*0.5_dp)
         write(*,'(A,I8,A,I8)') '    tf6 (Pg_ss water up) valid nodes = ', n_valid_fmm, ' / ', reg%nnode
      endif
   endif

   ! ============================================================
   ! tf7: Pg_ss water-down - Transfer tf6 to sea surface, reflect, FMM down
   ! ============================================================
   isec => intersection(1)  ! sea surface
   if (isec%nnode > 0 .and. s%n_time_fields >= 6) then
      call transfer_times_to_intersection(isec, s%time_field(6))

      reg => region(1)
      if (associated(reg)) then
         call reflect_gradient_at_interface(isec, reg, VTYPE_P)
         reg%node_status = NODE_FAR
         reg%arrivaltime = huge_time
         reg%time_gradient = 0.0_dp
         call sweep_region_from_interface2d(reg, isec, VTYPE_P, s)
         call propagate2d(reg, VTYPE_P)

         ! Fix top boundary gradient: dtdz must be > 0 for downward FMM
         do ix = 1, reg%grid%nx
            m_top = reg%region_rnode_id(ix, 1)
            m_below = reg%region_rnode_id(ix, 2)
            if (m_top > 0 .and. m_below > 0) then
               if (reg%time_gradient(2, m_top) < 0.0_dp .or. &
                   abs(reg%time_gradient(2, m_top)) < 0.01_dp) then
                  reg%time_gradient(1:2, m_top) = reg%time_gradient(1:2, m_below)
               endif
            endif
         end do

         s%n_time_fields = 7
         s%time_field(7)%reg => reg
         s%time_field(7)%vtype = VTYPE_P
         s%time_field(7)%istart => intersection(1)
         s%time_field(7)%inonstart => intersection(2)
         allocate(s%time_field(7)%arrivaltime(reg%nnode))
         allocate(s%time_field(7)%time_gradient(2, reg%nnode))
         s%time_field(7)%arrivaltime(:) = reg%arrivaltime(:)
         s%time_field(7)%time_gradient(:,:) = reg%time_gradient(:,:)

         n_valid_fmm = count(reg%arrivaltime < huge_time*0.5_dp)
         write(*,'(A,I8,A,I8)') '    tf7 (Pg_ss water down) valid nodes = ', n_valid_fmm, ' / ', reg%nnode
      endif
   endif

   ! ============================================================
   ! tf8: Pg_ss water-up 2nd - Transfer tf7 to seafloor, reflect, FMM up to receiver
   ! ============================================================
   isec => intersection(2)  ! seafloor
   if (isec%nnode > 0 .and. s%n_time_fields >= 7) then
      call transfer_times_to_intersection(isec, s%time_field(7))

      reg => region(1)
      if (associated(reg)) then
         call reflect_gradient_at_interface(isec, reg, VTYPE_P)
         reg%node_status = NODE_FAR
         reg%arrivaltime = huge_time
         reg%time_gradient = 0.0_dp
         call sweep_region_from_interface2d(reg, isec, VTYPE_P, s)
         call propagate2d(reg, VTYPE_P)

         ! Fix top boundary gradient: dtdz must be < 0 for upward FMM
         do ix = 1, reg%grid%nx
            m_top = reg%region_rnode_id(ix, 1)
            m_below = reg%region_rnode_id(ix, 2)
            if (m_top > 0 .and. m_below > 0) then
               if (reg%time_gradient(2, m_top) > 0.0_dp .or. &
                   abs(reg%time_gradient(2, m_top)) < 0.01_dp) then
                  reg%time_gradient(1:2, m_top) = reg%time_gradient(1:2, m_below)
               endif
            endif
         end do

         s%n_time_fields = 8
         s%time_field(8)%reg => reg
         s%time_field(8)%vtype = VTYPE_P
         s%time_field(8)%istart => intersection(2)
         s%time_field(8)%inonstart => intersection(1)
         allocate(s%time_field(8)%arrivaltime(reg%nnode))
         allocate(s%time_field(8)%time_gradient(2, reg%nnode))
         s%time_field(8)%arrivaltime(:) = reg%arrivaltime(:)
         s%time_field(8)%time_gradient(:,:) = reg%time_gradient(:,:)

         s%path(4)%tf_sequence(1) = 1  ! solid path (tf1)
         s%path(4)%tf_sequence(2) = 6  ! water up 1st
         s%path(4)%tf_sequence(3) = 7  ! water down (sea surface reflection)
         s%path(4)%tf_sequence(4) = 8  ! water up 2nd (seafloor reflection)

         n_valid_fmm = count(reg%arrivaltime < huge_time*0.5_dp)
         write(*,'(A,I8,A,I8)') '    tf8 (Pg_ss water up 2nd) valid nodes = ', n_valid_fmm, ' / ', reg%nnode
      endif
   endif

   ! ============================================================
   ! tf9: Pg_ss2 water-up — Point-source FMM in water from seafloor
   ! ============================================================
   reg => region(1)
   if (associated(reg)) then
      call initialize_region_from_source2d(reg, s%x, s%z, 0.0_dp, VTYPE_P)
      call propagate2d(reg, VTYPE_P)

      ! Fix top boundary gradient: dtdz must be < 0 for upward FMM
      do ix = 1, reg%grid%nx
         m_top = reg%region_rnode_id(ix, 1)
         m_below = reg%region_rnode_id(ix, 2)
         if (m_top > 0 .and. m_below > 0) then
            if (reg%time_gradient(2, m_top) > 0.0_dp .or. &
                abs(reg%time_gradient(2, m_top)) < 0.01_dp) then
               reg%time_gradient(1:2, m_top) = reg%time_gradient(1:2, m_below)
            endif
         endif
      end do

      s%n_time_fields = 9
      s%time_field(9)%reg => reg
      s%time_field(9)%vtype = VTYPE_P
      s%time_field(9)%istart => intersection(2)
      s%time_field(9)%inonstart => intersection(1)
      allocate(s%time_field(9)%arrivaltime(reg%nnode))
      allocate(s%time_field(9)%time_gradient(2, reg%nnode))
      s%time_field(9)%arrivaltime(:) = reg%arrivaltime(:)
      s%time_field(9)%time_gradient(:,:) = reg%time_gradient(:,:)

      n_valid_fmm = count(reg%arrivaltime < huge_time*0.5_dp)
      write(*,'(A,I8,A,I8)') '    tf9 (Pg_ss2 water up) valid nodes = ', n_valid_fmm, ' / ', reg%nnode
   endif

   ! ============================================================
   ! tf10: Pg_ss2 water-down — Transfer tf9 to sea surface, reflect, FMM down
   ! ============================================================
   isec => intersection(1)  ! sea surface
   if (isec%nnode > 0 .and. s%n_time_fields >= 9) then
      call transfer_times_to_intersection(isec, s%time_field(9))

      reg => region(1)
      if (associated(reg)) then
         call reflect_gradient_at_interface(isec, reg, VTYPE_P)
         reg%node_status = NODE_FAR
         reg%arrivaltime = huge_time
         reg%time_gradient = 0.0_dp
         call sweep_region_from_interface2d(reg, isec, VTYPE_P, s)
         call propagate2d(reg, VTYPE_P)

         ! Fix top boundary gradient: dtdz must be > 0 for downward FMM
         do ix = 1, reg%grid%nx
            m_top = reg%region_rnode_id(ix, 1)
            m_below = reg%region_rnode_id(ix, 2)
            if (m_top > 0 .and. m_below > 0) then
               if (reg%time_gradient(2, m_top) < 0.0_dp .or. &
                   abs(reg%time_gradient(2, m_top)) < 0.01_dp) then
                  reg%time_gradient(1:2, m_top) = reg%time_gradient(1:2, m_below)
               endif
            endif
         end do

         s%n_time_fields = 10
         s%time_field(10)%reg => reg
         s%time_field(10)%vtype = VTYPE_P
         s%time_field(10)%istart => intersection(1)
         s%time_field(10)%inonstart => intersection(2)
         allocate(s%time_field(10)%arrivaltime(reg%nnode))
         allocate(s%time_field(10)%time_gradient(2, reg%nnode))
         s%time_field(10)%arrivaltime(:) = reg%arrivaltime(:)
         s%time_field(10)%time_gradient(:,:) = reg%time_gradient(:,:)

         n_valid_fmm = count(reg%arrivaltime < huge_time*0.5_dp)
         write(*,'(A,I8,A,I8)') '    tf10 (Pg_ss2 water down) valid nodes = ', n_valid_fmm, ' / ', reg%nnode
      endif
   endif

   ! ============================================================
   ! tf11: Pg_ss2 solid diving — Transfer tf10 to seafloor, FMM in solid
   ! ============================================================
   isec => intersection(2)  ! seafloor
   if (isec%nnode > 0 .and. s%n_time_fields >= 10) then
      call transfer_times_to_intersection(isec, s%time_field(10))

      reg => region(5)  ! solid crust
      if (associated(reg)) then
         call reflect_gradient_at_interface(isec, reg, VTYPE_P)
         reg%node_status = NODE_FAR
         reg%arrivaltime = huge_time
         reg%time_gradient = 0.0_dp
         call sweep_region_from_interface2d(reg, isec, VTYPE_P, s)
         call propagate2d(reg, VTYPE_P)

         s%n_time_fields = 11
         s%time_field(11)%reg => reg
         s%time_field(11)%vtype = VTYPE_P
         s%time_field(11)%istart => intersection(2)
         allocate(s%time_field(11)%arrivaltime(reg%nnode))
         allocate(s%time_field(11)%time_gradient(2, reg%nnode))
         s%time_field(11)%arrivaltime(:) = reg%arrivaltime(:)
         s%time_field(11)%time_gradient(:,:) = reg%time_gradient(:,:)

         n_valid_fmm = count(reg%arrivaltime < huge_time*0.5_dp)
         write(*,'(A,I8,A,I8)') '    tf11 (Pg_ss2 solid) valid nodes = ', n_valid_fmm, ' / ', reg%nnode
      endif
   endif

   ! ============================================================
   ! tf12: Pg_ss2 water-up — Transfer tf11 to seafloor, FMM up to receiver
   ! ============================================================
   isec => intersection(2)  ! seafloor
   if (isec%nnode > 0 .and. s%n_time_fields >= 11) then
      call transfer_times_to_intersection(isec, s%time_field(11))

      reg => region(1)  ! water
      if (associated(reg)) then
         call reflect_gradient_at_interface(isec, reg, VTYPE_P)
         reg%node_status = NODE_FAR
         reg%arrivaltime = huge_time
         reg%time_gradient = 0.0_dp
         call sweep_region_from_interface2d(reg, isec, VTYPE_P, s)
         call propagate2d(reg, VTYPE_P)

         ! Fix top boundary gradient: dtdz must be < 0 for upward FMM
         do ix = 1, reg%grid%nx
            m_top = reg%region_rnode_id(ix, 1)
            m_below = reg%region_rnode_id(ix, 2)
            if (m_top > 0 .and. m_below > 0) then
               if (reg%time_gradient(2, m_top) > 0.0_dp .or. &
                   abs(reg%time_gradient(2, m_top)) < 0.01_dp) then
                  reg%time_gradient(1:2, m_top) = reg%time_gradient(1:2, m_below)
               endif
            endif
         end do

         s%n_time_fields = 12
         s%time_field(12)%reg => reg
         s%time_field(12)%vtype = VTYPE_P
         s%time_field(12)%istart => intersection(2)
         s%time_field(12)%inonstart => intersection(1)
         allocate(s%time_field(12)%arrivaltime(reg%nnode))
         allocate(s%time_field(12)%time_gradient(2, reg%nnode))
         s%time_field(12)%arrivaltime(:) = reg%arrivaltime(:)
         s%time_field(12)%time_gradient(:,:) = reg%time_gradient(:,:)

         s%path(5)%tf_sequence(1) = 9   ! water up
         s%path(5)%tf_sequence(2) = 10  ! water down
         s%path(5)%tf_sequence(3) = 11  ! solid diving
         s%path(5)%tf_sequence(4) = 12  ! water up to receiver

         n_valid_fmm = count(reg%arrivaltime < huge_time*0.5_dp)
         write(*,'(A,I8,A,I8)') '    tf12 (Pg_ss2 water up) valid nodes = ', n_valid_fmm, ' / ', reg%nnode
      endif
   endif

   ! ============================================================
   ! Trace rays
   ! ============================================================
   do n = 1, n_receivers
      do m = 1, receiver(n)%n_rays
         if (receiver(n)%ray(m)%source_id /= ns) cycle
         path_id = receiver(n)%ray(m)%raypath_id

         ! Check if all required time fields have been computed.
         ! Note: Cannot rely on s%n_time_fields alone because it may skip values
         ! when a region is pinched out (e.g., sediment layer nnode=0).
         block
            logical :: all_tf_ready
            integer :: itf_check, k_tf
            all_tf_ready = .true.
            do k_tf = 1, s%path(path_id)%n_tf
               itf_check = s%path(path_id)%tf_sequence(k_tf)
               if (itf_check < 1 .or. itf_check > s%n_time_fields) then
                  all_tf_ready = .false.
                  exit
               end if
               if (.not. associated(s%time_field(itf_check)%arrivaltime)) then
                  all_tf_ready = .false.
                  exit
               end if
            end do
            if (s%path(path_id)%valid .and. all_tf_ready) then
               ! Use the time field of the last section in the path
               itf = s%path(path_id)%tf_sequence(s%path(path_id)%n_tf)

            ! For Pg (path_id==1) in water: check if direct water path is faster
            if (path_id == 1 .and. in_water(n)) then
               dist_direct = sqrt((rec_x_orig(n) - s%x)**2 + (rec_z_orig(n) - s%z)**2)
               t_direct_water = dist_direct / vwater

               ! Find the best Pg entry point and compute Pg travel time
               call find_optimal_entry(s, ns, n, path_id, itf)
               entry_x_best(n) = receiver(n)%x
               entry_z_best(n) = receiver(n)%z
               entry_x_arr(n, path_id) = receiver(n)%x
               entry_z_arr(n, path_id) = receiver(n)%z

               ! Get Pg travel time at the entry point
               tf_ptr => s%time_field(itf)
               if (associated(tf_ptr%reg%arrivaltime)) then
                  call interpolate_arrivaltime2d(tf_ptr%reg, receiver(n)%x, receiver(n)%z, t_pg_fmm)
               else
                  t_pg_fmm = huge_time
               endif
               t_pg_total = t_pg_fmm + water_time(n)

               if (t_direct_water <= t_pg_total) then
                  ! Direct water path is faster: use it as Pg
                  receiver(n)%ray(m)%valid = .true.
                  receiver(n)%ray(m)%direct_water = .true.
                  receiver(n)%ray(m)%receiver_time = t_direct_water
                  water_time(n) = 0.0_dp  ! already included in direct path
                  n_valid_pg = n_valid_pg + 1
                  ! Restore receiver position
                  receiver(n)%x = rec_x_orig(n)
                  receiver(n)%z = interpolate_iface_z(intrface(2), rec_x_orig(n))
                  cycle
               endif

               ! Pg through solid is faster: trace ray normally
               call trace_ray_from_receiver2d(receiver(n), s, receiver(n)%ray(m))

               if (in_water(n)) then
                  receiver(n)%x = rec_x_orig(n)
                  receiver(n)%z = interpolate_iface_z(intrface(2), rec_x_orig(n))
               endif

               if (receiver(n)%ray(m)%valid) then
                  receiver(n)%ray(m)%receiver_time = receiver(n)%ray(m)%receiver_time &
                       + water_time(n)
                  n_valid_pg = n_valid_pg + 1
               endif
            else if (path_id == 4) then
               ! Pg_ss: 4-section path ends at sea surface (receiver position).
               ! Receiver was moved to seafloor at init — restore to original
               ! water-column position so tf8 can be interpolated directly.
               ! No water_time needed — tf8 already includes the full path.
               receiver(n)%x = rec_x_orig(n)
               receiver(n)%z = rec_z_orig(n)
               call trace_ray_from_receiver2d(receiver(n), s, receiver(n)%ray(m))

               if (receiver(n)%ray(m)%valid) then
                  n_valid_pgss = n_valid_pgss + 1
               endif
            else if (path_id == 5) then
               ! Pg_ss2: 4-section path ends at sea surface (receiver position).
               ! Same handling as Pg_ss — restore receiver to original position.
               receiver(n)%x = rec_x_orig(n)
               receiver(n)%z = rec_z_orig(n)
               call trace_ray_from_receiver2d(receiver(n), s, receiver(n)%ray(m))

               if (receiver(n)%ray(m)%valid) then
                  n_valid_pgss2 = n_valid_pgss2 + 1
               endif
            else
               ! Non-Pg, non-Pg_ss or non-water receiver: original logic
               if (in_water(n)) then
                  call find_optimal_entry(s, ns, n, path_id, itf)
                  entry_x_best(n) = receiver(n)%x
                  entry_z_best(n) = receiver(n)%z
                  entry_x_arr(n, path_id) = receiver(n)%x
                  entry_z_arr(n, path_id) = receiver(n)%z
               endif

               call trace_ray_from_receiver2d(receiver(n), s, receiver(n)%ray(m))

               if (in_water(n)) then
                  receiver(n)%x = rec_x_orig(n)
                  receiver(n)%z = interpolate_iface_z(intrface(2), rec_x_orig(n))
               endif

               if (receiver(n)%ray(m)%valid) then
                  receiver(n)%ray(m)%receiver_time = receiver(n)%ray(m)%receiver_time &
                       + water_time(n)
                  receiver(n)%ray(m)%water_time = water_time(n)
                  if (path_id == 1) n_valid_pg = n_valid_pg + 1
                  if (path_id == 2) n_valid_psedp = n_valid_psedp + 1
                  if (path_id == 3) n_valid_pmp = n_valid_pmp + 1
               endif
            endif
         else
            receiver(n)%ray(m)%valid = .false.
         endif
         end block
      end do
   end do

   write(*,'(A)') '  Sources done'

   if (write_rays) call write_rays_for_source(s, ns)

   do n = 1, s%n_time_fields
      if (associated(s%time_field(n)%arrivaltime)) deallocate(s%time_field(n)%arrivaltime)
      if (associated(s%time_field(n)%time_gradient)) deallocate(s%time_field(n)%time_gradient)
   end do

   write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)') '  Source ', ns, ' done. Pg=', n_valid_pg, &
      ' PsedP=', n_valid_psedp, ' PmP=', n_valid_pmp, ' Pg_ss=', n_valid_pgss, &
      ' Pg_ss2=', n_valid_pgss2

end do

call system_clock(t4, count_rate, count_max)

write(*,*)
write(*,'(A)') '  Forward modeling complete'
write(*,'(A,F12.2,A)') '  Init time:       ', dble(t2-t1)/dble(count_rate), ' sec'
write(*,'(A,F12.2,A)') '  Setup time:      ', dble(t3-t2)/dble(count_rate), ' sec'
write(*,'(A,F12.2,A)') '  Propagation time:', dble(t4-t3)/dble(count_rate), ' sec'
write(*,'(A,I0)') '  Valid Pg rays:    ', n_valid_pg
write(*,'(A,I0)') '  Valid PsedP rays: ', n_valid_psedp
write(*,'(A,I0)') '  Valid PmP rays:   ', n_valid_pmp
write(*,'(A,I0)') '  Valid Pg_ss rays: ', n_valid_pgss
write(*,'(A,I0)') '  Valid Pg_ss2 rays: ', n_valid_pgss2

call write_arrivals2d()

write(*,*)
write(*,'(A)') '  FMT2D finished'

contains

! ============================================================
! Helper subroutines
! ============================================================

function interpolate_iface_z(iface, px) result(z_val)
  type(Tinterface2d) :: iface
  real(kind=dp) :: px, z_val
  integer :: ixx
  real(kind=dp) :: w

  if (iface%nx < 2) then
     z_val = iface%z_nodes(1)
     return
  endif

  ixx = floor((px - iface%x(1)) / iface%dx) + 1
  if (ixx < 1) then
     z_val = iface%z_nodes(1)
     return
  endif
  if (ixx >= iface%nx) then
     z_val = iface%z_nodes(iface%nx)
     return
  endif

  w = (px - iface%x(ixx)) / iface%dx
  if (w < 0.0_dp) w = 0.0_dp
  if (w > 1.0_dp) w = 1.0_dp

  z_val = (1.0_dp - w)*iface%z_nodes(ixx) + w*iface%z_nodes(ixx+1)
end function interpolate_iface_z


subroutine transfer_times_to_intersection(isec, tf)
  type(Tintersection2d) :: isec
  type(Ttime_field2d) :: tf
  integer :: n, m
  logical :: copied

  if (.not. associated(isec%arrivaltime)) then
     allocate(isec%arrivaltime(isec%nnode))
     allocate(isec%time_gradient(2, isec%nnode))
  endif

  isec%arrivaltime = huge_time
  isec%time_gradient = 0.0_dp

  if (.not. associated(tf%reg)) return

  copied = .false.

  ! Use per-region boundary node IDs (safe from overwrite by later define_region2d calls)
  if (isec%id == tf%reg%ibot%id) then
     copied = .true.
     do n = 1, isec%nnode
        m = 0
        if (associated(tf%reg%ibot_rnode_id)) m = tf%reg%ibot_rnode_id(n)
         if (m > 0 .and. m <= tf%reg%nnode) then
            if (associated(tf%arrivaltime) .and. m <= size(tf%arrivaltime)) then
               isec%arrivaltime(n) = tf%arrivaltime(m)
            else
               isec%arrivaltime(n) = tf%reg%arrivaltime(m)
            endif
            if (associated(tf%time_gradient) .and. m <= size(tf%time_gradient, 2)) then
               isec%time_gradient(1:2, n) = tf%time_gradient(1:2, m)
            else
               isec%time_gradient(1:2, n) = tf%reg%time_gradient(1:2, m)
            endif
         endif
     end do
  else if (isec%id == tf%reg%itop%id) then
     copied = .true.
     do n = 1, isec%nnode
        m = 0
        if (associated(tf%reg%itop_rnode_id)) m = tf%reg%itop_rnode_id(n)
         if (m > 0 .and. m <= tf%reg%nnode) then
            if (associated(tf%arrivaltime) .and. m <= size(tf%arrivaltime)) then
               isec%arrivaltime(n) = tf%arrivaltime(m)
            else
               isec%arrivaltime(n) = tf%reg%arrivaltime(m)
            endif
            if (associated(tf%time_gradient) .and. m <= size(tf%time_gradient, 2)) then
               isec%time_gradient(1:2, n) = tf%time_gradient(1:2, m)
            else
               isec%time_gradient(1:2, n) = tf%reg%time_gradient(1:2, m)
            endif
         endif
     end do
  endif

  if (.not. copied) then
     do n = 1, isec%nnode
        call sample_region_time(tf%reg, isec%x(n), isec%z(n), &
             isec%arrivaltime(n), isec%time_gradient(1:2,n))
     end do
  endif

end subroutine transfer_times_to_intersection


subroutine sample_region_time(reg, px, pz, atime, grad)
  type(Tregion2d) :: reg
  real(kind=dp) :: px, pz, atime, grad(2)
  type(Tpropagation_grid2d), pointer :: grid
  integer :: ix, iz, n1, n2, n3, n4, nn
  real(kind=dp) :: w1, w2

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

  n1 = reg%region_rnode_id(ix, iz)
  n2 = reg%region_rnode_id(ix+1, iz)
  n3 = reg%region_rnode_id(ix, iz+1)
  n4 = reg%region_rnode_id(ix+1, iz+1)

  if (n1 > 0 .and. n2 > 0 .and. n3 > 0 .and. n4 > 0) then
     atime = (1.0_dp-w1)*(1.0_dp-w2)*reg%arrivaltime(n1) + &
          w1*(1.0_dp-w2)*reg%arrivaltime(n2) + &
          (1.0_dp-w1)*w2*reg%arrivaltime(n3) + &
          w1*w2*reg%arrivaltime(n4)
     grad = (1.0_dp-w1)*(1.0_dp-w2)*reg%time_gradient(1:2,n1) + &
          w1*(1.0_dp-w2)*reg%time_gradient(1:2,n2) + &
          (1.0_dp-w1)*w2*reg%time_gradient(1:2,n3) + &
          w1*w2*reg%time_gradient(1:2,n4)
  else
     nn = nearest_region_node(reg, px, pz)
     if (nn > 0) then
        atime = reg%arrivaltime(nn)
        grad = reg%time_gradient(1:2, nn)
     else
        atime = huge_time
        grad = 0.0_dp
     endif
  endif
end subroutine sample_region_time


function nearest_region_node(reg, px, pz) result(nn)
  type(Tregion2d) :: reg
  real(kind=dp) :: px, pz
  integer :: nn, i
  real(kind=dp) :: dist, distmin

  nn = 0
  distmin = 1.0e100_dp
  do i = 1, reg%nnode
     dist = (px - reg%x(i))**2 + (pz - reg%z(i))**2
     if (dist < distmin) then
        distmin = dist
        nn = i
     endif
  end do
end function nearest_region_node


subroutine reflect_gradient_at_interface(isec, target_reg, target_vtype)
  type(Tintersection2d) :: isec
  type(Tregion2d) :: target_reg
  integer :: target_vtype
  integer :: n, target_side
  real(kind=dp) :: norm(2), h, grad_perp, grad_par(2)
  real(kind=dp) :: v_target, p_total, p_par2, p_norm2

  do n = 1, isec%nnode
     call interface_normal2d(isec%x(n), isec, norm(1), norm(2), h)

     if (isec%arrivaltime(n) >= huge_time*0.5_dp) then
        isec%time_gradient(1:2,n) = 0.0_dp
        cycle
     endif

     grad_perp = dot_product(isec%time_gradient(1:2, n), norm)
     grad_par = isec%time_gradient(1:2, n) - grad_perp * norm

     target_side = interface_target_side(isec, target_reg, grad_perp)
     v_target = interface_target_velocity(isec, target_reg, target_vtype, n, target_side)
     if (v_target <= 0.01_dp) then
        isec%arrivaltime(n) = huge_time
        isec%time_gradient(1:2,n) = 0.0_dp
        cycle
     endif

     p_total = 1.0_dp / v_target
     p_par2 = sum(grad_par**2)
     p_norm2 = p_total*p_total - p_par2

     if (p_norm2 < -1.0e-10_dp) then
        isec%arrivaltime(n) = huge_time
        isec%time_gradient(1:2,n) = 0.0_dp
     else
        if (p_norm2 < 0.0_dp) p_norm2 = 0.0_dp
        isec%time_gradient(1:2,n) = grad_par + &
             target_side * sqrt(p_norm2) * norm
     endif
  end do

end subroutine reflect_gradient_at_interface


function interface_target_side(isec, target_reg, grad_perp) result(side)
  type(Tintersection2d) :: isec
  type(Tregion2d) :: target_reg
  real(kind=dp) :: grad_perp
  integer :: side

  if (associated(target_reg%itop)) then
     if (target_reg%itop%id == isec%id) then
        side = 1
        return
     endif
  endif
  if (associated(target_reg%ibot)) then
     if (target_reg%ibot%id == isec%id) then
        side = -1
        return
     endif
  endif

  if (grad_perp < 0.0_dp) then
     side = -1
  else
     side = 1
  endif
end function interface_target_side


function interface_target_velocity(isec, target_reg, vtype, inode, side) result(vel)
  type(Tintersection2d) :: isec
  type(Tregion2d) :: target_reg
  integer :: vtype, inode, side
  real(kind=dp) :: vel
  integer :: m

  vel = 0.0_dp
  if (side > 0) then
     if (associated(isec%vel_bot)) vel = isec%vel_bot(vtype, inode)
  else
     if (associated(isec%vel_top)) vel = isec%vel_top(vtype, inode)
  endif

  if (vel > 0.01_dp) return

  m = 0
  if (associated(target_reg%itop)) then
     if (target_reg%itop%id == isec%id .and. associated(target_reg%itop_rnode_id)) &
          m = target_reg%itop_rnode_id(inode)
  endif
  if (m == 0 .and. associated(target_reg%ibot)) then
     if (target_reg%ibot%id == isec%id .and. associated(target_reg%ibot_rnode_id)) &
          m = target_reg%ibot_rnode_id(inode)
  endif

  if (m > 0 .and. m <= target_reg%nnode) vel = target_reg%velocity(vtype, m)
end function interface_target_velocity


subroutine write_arrivals2d()
  integer :: iunit, iunit2, ircv, iray, nvalid
  real(kind=dp) :: t_arr, t_solid

  iunit = 40
  open(unit=iunit, file=trim(arrivals_file), status='replace', action='write')

  nvalid = 0
  do ircv = 1, n_receivers
     do iray = 1, receiver(ircv)%n_rays
        if (receiver(ircv)%ray(iray)%valid) nvalid = nvalid + 1
     end do
  end do

  iunit2 = 41
  open(unit=iunit2, file=trim(obs_file), status='replace', action='write')
  write(iunit2, '(I0)') nvalid

  do ircv = 1, n_receivers
     do iray = 1, receiver(ircv)%n_rays
        if (receiver(ircv)%ray(iray)%valid) then
           t_arr = receiver(ircv)%ray(iray)%receiver_time
           t_solid = t_arr - receiver(ircv)%ray(iray)%water_time
        else
           t_arr = -1.0_dp
           t_solid = -1.0_dp
        endif
        write(iunit, '(3I8,2F15.6,F12.6)') ircv, &
             receiver(ircv)%ray(iray)%source_id, &
             receiver(ircv)%ray(iray)%raypath_id, &
             t_arr, t_solid, receiver(ircv)%ray(iray)%water_time
        if (receiver(ircv)%ray(iray)%valid) then
           write(iunit2, '(3I8,F15.6)') ircv, &
                receiver(ircv)%ray(iray)%source_id, &
                receiver(ircv)%ray(iray)%raypath_id, &
                t_arr
        endif
     end do
  end do

  close(iunit)
  close(iunit2)

  ! Write entry points file for inversion (entry_x, entry_z, water_time per ray)
  open(unit=42, file='entry_points.dat', status='replace', action='write')
  write(42, '(I0)') nvalid
  do ircv = 1, n_receivers
     do iray = 1, receiver(ircv)%n_rays
        if (receiver(ircv)%ray(iray)%valid) then
           write(42, '(3I8,2F15.6,F12.6)') ircv, &
                receiver(ircv)%ray(iray)%source_id, &
                receiver(ircv)%ray(iray)%raypath_id, &
                entry_x_arr(ircv, receiver(ircv)%ray(iray)%raypath_id), &
                entry_z_arr(ircv, receiver(ircv)%ray(iray)%raypath_id), &
                receiver(ircv)%ray(iray)%water_time
        endif
     end do
  end do
  close(42)

  write(*,'(A,A,A,I0,A)') '  Observed times written to ', trim(obs_file), ' (', nvalid, ' valid)'
  write(*,'(A,A)') '  Arrivals written to ', trim(arrivals_file)
  write(*,'(A)') '  Entry points written to entry_points.dat'

end subroutine write_arrivals2d


subroutine write_rays_for_source(s, src_id)
  type(Tsource2d), target :: s
  integer :: src_id
  integer :: n, m, k, iunit, isec
  real(kind=dp) :: x_entry, z_sf, x_rcv, z_rcv
  integer :: npts_water
  character(len=64) :: fname

  write(fname, '(A,A,I0,A)') trim(rays_prefix), '_src', src_id, '.dat'
  iunit = 50 + src_id
  open(unit=iunit, file=trim(fname), status='replace', action='write')

  do n = 1, n_receivers
     do m = 1, receiver(n)%n_rays
        if (receiver(n)%ray(m)%source_id /= src_id) cycle
        if (.not. receiver(n)%ray(m)%valid) cycle

        write(iunit, '(A,I0,1X,I0,1X,I0,1X,F12.6)') '> ', src_id, n, &
             receiver(n)%ray(m)%raypath_id, &
             receiver(n)%ray(m)%receiver_time

        ! Direct water path: receiver -> source straight line
        if (receiver(n)%ray(m)%direct_water) then
           x_rcv = rec_x_orig(n)
           z_rcv = rec_z_orig(n)
           ! Recompute dist_direct for this receiver (avoid stale global value)
           dist_direct = sqrt((x_rcv - s%x)**2 + (z_rcv - s%z)**2)
           npts_water = max(3, nint(dist_direct / 0.1_dp) + 1)
           do k = 1, npts_water
              write(iunit, '(F12.6,1X,F12.6)') &
                   x_rcv + (s%x - x_rcv) * dble(k-1) / dble(npts_water-1), &
                   z_rcv + (s%z - z_rcv) * dble(k-1) / dble(npts_water-1)
           end do
           write(iunit, *)
           cycle
        endif

        x_entry = rec_x_orig(n)
        z_sf = interpolate_iface_z(intrface(2), rec_x_orig(n))
        if (in_water(n)) then
           ! Get the entry point from the ray's last section (closest to receiver).
           ! This avoids the bug where entry_x_best(n) is overwritten by later rays
           ! for the same receiver.
           isec = receiver(n)%ray(m)%nsections
           if (associated(receiver(n)%ray(m)%section(isec)%point)) then
              x_entry = receiver(n)%ray(m)%section(isec)%point(1,1)
              z_sf = receiver(n)%ray(m)%section(isec)%point(2,1)
           else
              x_entry = entry_x_best(n)
              z_sf = entry_z_best(n)
           endif
        endif
        x_rcv = rec_x_orig(n)
        z_rcv = rec_z_orig(n)

        if (z_sf > z_rcv + 0.001_dp) then
           write(iunit, '(F12.6,1X,F12.6)') x_rcv, z_rcv
           npts_water = max(2, nint(abs(z_sf - z_rcv) / 0.05_dp) + 1)
           do k = 1, npts_water
              write(iunit, '(F12.6,1X,F12.6)') &
                   x_rcv + (x_entry - x_rcv) * dble(k) / dble(npts_water), &
                   z_rcv + (z_sf - z_rcv) * dble(k) / dble(npts_water)
           end do
        endif

        do isec = receiver(n)%ray(m)%nsections, 1, -1
           if (.not. associated(receiver(n)%ray(m)%section(isec)%point)) cycle
           do k = 1, receiver(n)%ray(m)%section(isec)%npoints
              write(iunit, '(F12.6,1X,F12.6)') &
                   receiver(n)%ray(m)%section(isec)%point(1,k), &
                   receiver(n)%ray(m)%section(isec)%point(2,k)
           end do
        end do
        write(iunit, *)
     end do
  end do

  close(iunit)
end subroutine write_rays_for_source


subroutine find_optimal_entry(s, src_id, ircv, path_id, itf)
  type(Tsource2d), target :: s
  integer :: src_id, ircv, path_id, itf
  integer :: i, nseafloor, nrefine, iref
  real(kind=dp) :: x_test, z_sf, t_fmm, t_water, t_total, t_best
  real(kind=dp) :: dist_water, dist_src, x_best, z_best, wt_best
  real(kind=dp) :: x_left, x_right, dx_refine
  real(kind=dp) :: rx, rz, min_src_dist
  type(Ttime_field2d), pointer :: tf
  type(Tregion2d), pointer :: tf_reg

  rx = rec_x_orig(ircv)
  rz = rec_z_orig(ircv)

  tf => s%time_field(itf)
  if (.not. associated(tf%reg)) return

  tf_reg => tf%reg

  if (associated(tf%arrivaltime)) then
     tf_reg%arrivaltime = tf%arrivaltime
     tf_reg%time_gradient = tf%time_gradient
  endif

  nseafloor = intrface(2)%nx

  ! Minimum distance from source: skip seafloor points too close to the source
  ! where the FMM gradient is zero or near-zero (causes ray tracing failure)
  min_src_dist = 3.0_dp * pgrid%dx

  t_best = huge_time
  x_best = rx
  z_best = interpolate_iface_z(intrface(2), rx)
  wt_best = (z_best - rz) / vwater

  do i = 1, nseafloor
     x_test = intrface(2)%x(i)
     z_sf = intrface(2)%z_nodes(i)

     ! Skip points too close to the source (gradient is zero there)
     dist_src = sqrt((x_test - s%x)**2 + (z_sf - s%z)**2)
     if (dist_src < min_src_dist) cycle

     call interpolate_arrivaltime2d(tf_reg, x_test, z_sf, t_fmm)
     if (t_fmm > huge_time * 0.5) cycle

     dist_water = sqrt((x_test - rx)**2 + (z_sf - rz)**2)
     t_water = dist_water / vwater
     t_total = t_fmm + t_water

     if (t_total < t_best) then
        t_best = t_total
        x_best = x_test
        z_best = z_sf
        wt_best = t_water
     endif
  end do

  if (x_best <= intrface(2)%x(1) .or. x_best >= intrface(2)%x(nseafloor)) then
     receiver(ircv)%x = x_best
     receiver(ircv)%z = z_best
     water_time(ircv) = wt_best
     return
  endif

  do i = 1, nseafloor - 1
     if (x_best >= intrface(2)%x(i) .and. x_best <= intrface(2)%x(i+1)) then
        x_left = intrface(2)%x(i)
        x_right = intrface(2)%x(i+1)
        exit
     endif
  end do

  nrefine = 10
  dx_refine = (x_right - x_left) / dble(nrefine)
  do iref = 1, nrefine - 1
     x_test = x_left + dble(iref) * dx_refine
     z_sf = interpolate_iface_z(intrface(2), x_test)

     ! Skip points too close to the source
     dist_src = sqrt((x_test - s%x)**2 + (z_sf - s%z)**2)
     if (dist_src < min_src_dist) cycle

     call interpolate_arrivaltime2d(tf_reg, x_test, z_sf, t_fmm)
     if (t_fmm > huge_time * 0.5) cycle

     dist_water = sqrt((x_test - rx)**2 + (z_sf - rz)**2)
     t_water = dist_water / vwater
     t_total = t_fmm + t_water

     if (t_total < t_best) then
        t_best = t_total
        x_best = x_test
        z_best = z_sf
        wt_best = t_water
     endif
  end do

  receiver(ircv)%x = x_best
  receiver(ircv)%z = z_best
  water_time(ircv) = wt_best

end subroutine find_optimal_entry


end program fmt2d