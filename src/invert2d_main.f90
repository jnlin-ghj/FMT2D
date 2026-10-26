program invert2d

use type_definitions2d
use global_variables2d
use mod_rays2d
use read_tomo2d
use mod_fmt2d
use mod_frechet2d
use mod_inversion2d

implicit none

integer :: n, m, i, j, k, ns, ipath, itf, step, iz, ix
integer :: ntf, vtype, prev_tf
integer :: t1, t2, t3, t4, count_rate, count_max
integer :: path_id, nrec_total, iobs, idx
integer :: nobs_geom, nx_smesh, nz_smesh, nx_refl, nx_sed
integer :: ntf_max, n_paths_per_source
integer :: n_valid_pg, n_valid_psedp, n_valid_pmp
integer :: nn, n_valid_top, n_alive_sweep, n_valid_fmm
integer :: n_updated, mm, icell, ixx, izz, i1, i2, i3, ii, jj

real(kind=dp) :: vwater, vair, z_seafloor, z_sedbot, z_moho, t_arrival
real(kind=dp) :: dtdx, dtdz, dtdz_sq, dx_if, v_loc, dtdz_orig
real(kind=dp) :: pz_rel_node_wm, z_sf_node_wm, wx_z_node_wm

real(kind=dp), allocatable :: xpos(:), topo(:), zpos_rel(:), vel(:,:)
real(kind=dp), allocatable :: refl_x(:), refl_z(:)
real(kind=dp), allocatable :: sed_x(:), sed_z(:)
real(kind=dp), allocatable :: obs_x(:), obs_z(:)
real(kind=dp), allocatable :: rec_x(:), rec_z(:)
integer, allocatable :: obs_id(:), rec_phase(:)
real(kind=dp), allocatable :: geom_obs_time(:), geom_uncertainty(:)

real(kind=dp), allocatable :: water_corr(:)

type(Tsource2d), pointer :: s
type(Treceiver2d), pointer :: rec
type(Tray2d), pointer :: ray
type(Tregion2d), pointer :: reg
type(Tintersection2d), pointer :: isec

logical :: do_frechet_derivatives
integer :: niter_inv, lsqr_max_iter_in
real(kind=dp) :: lam_damp_in, lam_smooth_in, inv_crit_chi_in, target_pert_in
real(kind=dp) :: lam_damp_iface_in, lam_smooth_iface_in
real(kind=dp) :: refl_weight_in, phase_wt_pmp_in, phase_wt_pgss_in, phase_wt_psed_in
real(kind=dp) :: corr_x_in, corr_z_in, target_dd_in
real(kind=dp) :: iface_step_ratio_in, max_rel_pert_in
real(kind=dp) :: corr_x_iface_in, iface_dx_target_in
  integer :: alt_inv_mode_in
character(len=256) :: smesh_vp_file, refl_file, sed_file
character(len=256) :: geom_file, obs_file

call system_clock(t1, count_rate, count_max)

write(*,'(A)') '============================================'
write(*,'(A)') '  INVERT2D - Standalone Inversion Program'
write(*,'(A)') '============================================'
write(*,*)

! Read inversion parameters
open(1, file='invert2d.in', status='old', action='read')
read(1, '(A)') smesh_vp_file
read(1, '(A)') refl_file
read(1, '(A)') sed_file
read(1, '(A)') geom_file
obs_file = ''
read(1, *) niter_inv
read(1, *) lam_damp_in
read(1, *) lam_smooth_in
read(1, *) inv_crit_chi_in
read(1, *) lsqr_max_iter_in
read(1, *) target_pert_in
read(1, *) lam_damp_iface_in
read(1, *) lam_smooth_iface_in
read(1, *) refl_weight_in
read(1, *) phase_wt_pmp_in
read(1, *) phase_wt_pgss_in
read(1, *) phase_wt_psed_in
read(1, *) corr_x_in
read(1, *) corr_z_in
! Initialize optional parameters before conditional reads
target_dd_in = 0.0_dp
iface_step_ratio_in = 0.0_dp
max_rel_pert_in = 0.0_dp
alt_inv_mode_in = -1
corr_x_iface_in = corr_x_in
iface_dx_target_in = 0.0_dp
read(1, *, end=82) smooth_type, target_dd_in, iface_step_ratio_in, max_rel_pert_in
read(1, *, end=82) alt_inv_mode_in
read(1, *, end=82) corr_x_iface_in
read(1, *, end=82) iface_dx_target_in
82 continue
if (iface_step_ratio_in <= 0.0_dp) iface_step_ratio_in = 0.5_dp
if (max_rel_pert_in <= 0.0_dp) max_rel_pert_in = 0.30_dp
if (alt_inv_mode_in < 0 .or. alt_inv_mode_in > 2) alt_inv_mode_in = -1
if (corr_x_iface_in <= 0.0_dp) corr_x_iface_in = corr_x_in
close(1)

write(*,'(A,A)') '  smesh_vp file:    ', trim(smesh_vp_file)
write(*,'(A,A)') '  moho file:        ', trim(refl_file)
write(*,'(A,A)') '  sed_bottom file:  ', trim(sed_file)
write(*,'(A,A)') '  geom file:        ', trim(geom_file)
write(*,'(A,A)') '  obs_times file:   ', trim(obs_file)
write(*,'(A,I4)')     '  niter:            ', niter_inv
write(*,'(A,F8.4)')   '  lam_damp:         ', lam_damp_in
write(*,'(A,F8.4)')   '  lam_smooth:       ', lam_smooth_in
write(*,'(A,F8.4)')   '  inv_crit_chi:     ', inv_crit_chi_in
write(*,'(A,I6)')     '  lsqr_max_iter:    ', lsqr_max_iter_in
write(*,'(A,F8.4)')   '  target_pert:      ', target_pert_in
write(*,'(A,F8.4)')   '  lam_damp_iface:   ', lam_damp_iface_in
write(*,'(A,F8.4)')   '  lam_smooth_iface: ', lam_smooth_iface_in
write(*,'(A,F8.4)')   '  refl_weight:      ', refl_weight_in
write(*,'(A,F8.4)')   '  phase_wt_pmp:     ', phase_wt_pmp_in
write(*,'(A,F8.4)')   '  phase_wt_pgss:    ', phase_wt_pgss_in
write(*,'(A,F8.4)')   '  phase_wt_psed:    ', phase_wt_psed_in
write(*,'(A,F8.4)')   '  corr_x:           ', corr_x_in
write(*,'(A,F8.4)')   '  corr_z:           ', corr_z_in
write(*,'(A,I4)')     '  smooth_type:      ', smooth_type
write(*,'(A,F8.4)')   '  target_dd:        ', target_dd_in
write(*,'(A,F8.4)')   '  iface_step_ratio: ', iface_step_ratio_in
write(*,'(A,F8.4)')   '  max_rel_pert:     ', max_rel_pert_in
write(*,'(A,I4)')     '  alt_inv_mode:     ', alt_inv_mode_in
write(*,'(A,F8.4)')   '  corr_x_iface:     ', corr_x_iface_in
write(*,*)

! Read velocity mesh, reflection interface, sediment interface, geometry
call read_smesh(smesh_vp_file, nx_smesh, nz_smesh, vwater, xpos, topo, zpos_rel, vel, vair)
vwater_global = vwater
vair_global = vair
call read_refl(refl_file, nx_refl, refl_x, refl_z)
call read_refl(sed_file, nx_sed, sed_x, sed_z)
call read_geom(geom_file, nobs_geom, obs_x, obs_z, nrec_total, rec_x, rec_z, obs_id, rec_phase, geom_obs_time, geom_uncertainty)

write(*,*)
write(*,'(A)') '  Building FMT2D model from input...'

! ============================================================
! 5 interfaces: 1=water_surface, 2=seafloor, 3=sed_bottom, 4=Moho, 5=bottom
! ============================================================
n_interfaces = 5
allocate(intrface(n_interfaces))

do n = 1, n_interfaces
   call interface2d_defaults(intrface(n))
   intrface(n)%id = n
end do

! Interface 1: water surface (z=0)
intrface(1)%nx = nx_smesh
intrface(1)%nnode = nx_smesh
intrface(1)%dx = xpos(2) - xpos(1)
intrface(1)%x0 = xpos(1)
allocate(intrface(1)%x(nx_smesh), intrface(1)%z_nodes(nx_smesh))
allocate(intrface(1)%z_nodes_init(nx_smesh))
intrface(1)%x = xpos
intrface(1)%z_nodes = 0.0_dp
intrface(1)%z_nodes_init = 0.0_dp

! Interface 2: seafloor
intrface(2)%nx = nx_smesh
intrface(2)%nnode = nx_smesh
intrface(2)%dx = xpos(2) - xpos(1)
intrface(2)%x0 = xpos(1)
allocate(intrface(2)%x(nx_smesh), intrface(2)%z_nodes(nx_smesh))
allocate(intrface(2)%z_nodes_init(nx_smesh))
intrface(2)%x = xpos
intrface(2)%z_nodes = topo
intrface(2)%z_nodes_init = topo

! Interface 3: sediment bottom
intrface(3)%nx = nx_sed
intrface(3)%nnode = nx_sed
intrface(3)%dx = sed_x(2) - sed_x(1)
intrface(3)%x0 = sed_x(1)
allocate(intrface(3)%x(nx_sed), intrface(3)%z_nodes(nx_sed))
allocate(intrface(3)%z_nodes_init(nx_sed))
intrface(3)%x = sed_x
intrface(3)%z_nodes = sed_z
intrface(3)%z_nodes_init = sed_z

! Interface 4: Moho
intrface(4)%nx = nx_refl
intrface(4)%nnode = nx_refl
intrface(4)%dx = refl_x(2) - refl_x(1)
intrface(4)%x0 = refl_x(1)
allocate(intrface(4)%x(nx_refl), intrface(4)%z_nodes(nx_refl))
allocate(intrface(4)%z_nodes_init(nx_refl))
intrface(4)%x = refl_x
intrface(4)%z_nodes = refl_z
intrface(4)%z_nodes_init = refl_z

! Interface 5: bottom
intrface(5)%nx = nx_smesh
intrface(5)%nnode = nx_smesh
intrface(5)%dx = xpos(2) - xpos(1)
intrface(5)%x0 = xpos(1)
allocate(intrface(5)%x(nx_smesh), intrface(5)%z_nodes(nx_smesh))
allocate(intrface(5)%z_nodes_init(nx_smesh))
intrface(5)%x = xpos
intrface(5)%z_nodes = zpos_rel(nz_smesh) + maxval(topo)
intrface(5)%z_nodes_init = zpos_rel(nz_smesh) + maxval(topo)

write(*,'(A)') '  Interfaces set up'

! ============================================================
! Velocity grid (P-wave)
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

! Save original grid info before resampling
nz_orig = nz_smesh
allocate(z_orig(nz_smesh))
z_orig = zpos_rel

! Resample velocity to propagation grid
block
  integer :: nz_new, iz_new, iz_old, ix_new
  real(kind=dp) :: dz_new, z_val, wz_old
  real(kind=dp), allocatable :: zpos_new(:), vel_new(:,:)
  real(kind=dp) :: pz_rel_node, z_sf_node, wx_z_node
  integer :: ivx_z_node

  dz_new = pgrid%dz
  nz_new = pgrid%nz

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
end do

! ============================================================
! Regions: water(1), sediment(2), crust(3), mantle(4), combined_all(5),
!          combined_sed_crust(6)
! Region 5 = sediment+crust+mantle combined (seafloor to bottom)
!            Used for Pg (diving wave)
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
      region(n)%nnode = 0
      cycle
   endif

   call define_region2d(region(n), region(n)%itop, region(n)%ibot, pgrid)
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

! Update intersection region pointers for region 5
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

call tag_regular_nodes2d(pgrid)
call velocities_to_grid2d(pgrid)

do n = 1, n_intersections
   call velocities_to_intersection2d(intersection(n))
end do
do n = 1, n_regions
   if (region(n)%nnode > 0) then
      call velocities_to_region2d(region(n), pgrid)
   endif
end do

do n = 2, n_regions
   if (region(n)%nnode > 0) then
      call tag_active_vgrid_nodes2d(region(n))
   endif
end do

! Mask out water layer from active velocity grid nodes
do ixx = 1, nx_smesh
  i = ixx
  if (i < 1) i = 1
  if (i >= nx_smesh) i = nx_smesh - 1
  wx_z_node_wm = dble(ixx - i)
  if (wx_z_node_wm < 0.0_dp) wx_z_node_wm = 0.0_dp
  if (wx_z_node_wm > 1.0_dp) wx_z_node_wm = 1.0_dp
  z_sf_node_wm = (1.0_dp - wx_z_node_wm) * topo(i) + wx_z_node_wm * topo(i + 1)
  do izz = 1, nz_smesh
    pz_rel_node_wm = vgrid(1, VTYPE_P)%z(izz) - z_sf_node_wm
    if (pz_rel_node_wm < 0.0_dp) then
      vgrid(1, VTYPE_P)%active(VTYPE_P, ixx, izz) = .false.
    end if
  end do
end do

write(*,'(A,I8)') '  vg active nodes (after water mask) = ', count(vgrid(1, VTYPE_P)%active)

! Set remaining intersection region pointers
intersection(1)%regleft => null()
intersection(1)%regright => region(1)
! intersection(2): regleft stays region(1), regright already set to region(5)
intersection(3)%regleft => region(2)
intersection(3)%regright => region(3)
intersection(4)%regleft => region(3)
intersection(4)%regright => region(4)
! intersection(5): regleft already set to region(5)
intersection(5)%regright => null()

write(*,'(A)') '  Model setup complete'

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
entry_x_best = 0.0_dp
entry_z_best = 0.0_dp
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

! ============================================================
! Define 5 phases for P-wave inversion:
!   1 = Pg    (diving wave — global first arrival in region 5)
!   2 = PsedP (sediment bottom reflection)
!   3 = PmP   (Moho reflection)
!   4 = Pg_ss (Pg sea-surface multiple)
!   5 = Pg_ss2 (water-first multiple)
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

   ! Phase 1: Pg — diving wave (global first arrival)
   s%path(1)%n_tf = 1
   allocate(s%path(1)%sequence(2))
   allocate(s%path(1)%tf_sequence(1))
   allocate(s%path(1)%vtype_sequence(1))
   s%path(1)%sequence(1) = 0   ! source
   s%path(1)%sequence(2) = 2   ! seafloor
   s%path(1)%tf_sequence(1) = 1  ! tf1: Pg
   s%path(1)%vtype_sequence(1) = VTYPE_P
   s%path(1)%valid = .true.
   s%path(1)%used = .true.

   ! Phase 2: PsedP — sediment bottom reflection
   s%path(2)%n_tf = 2
   allocate(s%path(2)%sequence(4))
   allocate(s%path(2)%tf_sequence(2))
   allocate(s%path(2)%vtype_sequence(2))
   s%path(2)%sequence(1) = 0   ! source
   s%path(2)%sequence(2) = 3   ! sed_bot
   s%path(2)%sequence(3) = 3   ! sed_bot (reflected)
   s%path(2)%sequence(4) = 0   ! seafloor
   s%path(2)%tf_sequence(1) = 2  ! tf2: PsedP down
   s%path(2)%tf_sequence(2) = 3  ! tf3: PsedP up
   s%path(2)%vtype_sequence(1) = VTYPE_P
   s%path(2)%vtype_sequence(2) = VTYPE_P
   s%path(2)%valid = .true.
   s%path(2)%used = .true.

   ! Phase 3: PmP — Moho reflection
   s%path(3)%n_tf = 2
   allocate(s%path(3)%sequence(4))
   allocate(s%path(3)%tf_sequence(2))
   allocate(s%path(3)%vtype_sequence(2))
   s%path(3)%sequence(1) = 0   ! source
   s%path(3)%sequence(2) = 4   ! Moho
   s%path(3)%sequence(3) = 4   ! Moho (reflected)
   s%path(3)%sequence(4) = 0   ! seafloor
   s%path(3)%tf_sequence(1) = 4  ! tf4: PmP down
   s%path(3)%tf_sequence(2) = 5  ! tf5: PmP up
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

   ! Phase 5: Pg_ss2 — water-first multiple
   ! Path: Source → water up → sea surface reflect → water down → solid diving → water up → receiver
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
   s%path(5)%tf_sequence(1) = 9   ! water up
   s%path(5)%tf_sequence(2) = 10  ! water down
   s%path(5)%tf_sequence(3) = 11  ! solid diving
   s%path(5)%tf_sequence(4) = 12  ! water up to receiver
   s%path(5)%vtype_sequence(1) = VTYPE_P
   s%path(5)%vtype_sequence(2) = VTYPE_P
   s%path(5)%vtype_sequence(3) = VTYPE_P
   s%path(5)%vtype_sequence(4) = VTYPE_P
   s%path(5)%valid = .true.
   s%path(5)%used = .true.

   ntf_max = 12
   allocate(s%time_field(ntf_max))
   do n = 1, ntf_max
      call time_field2d_defaults(s%time_field(n))
      s%time_field(n)%id = n
   end do

   s%n_time_fields = 12

   ! tf1: Pg — Point-source FMM in region 5 (seafloor->bottom)
   if (region(5)%nnode > 0) then
      s%time_field(1)%reg => region(5)
      s%time_field(1)%vtype = VTYPE_P
      s%time_field(1)%istart => intersection(2)
      s%time_field(1)%inonstart => intersection(5)
      allocate(s%time_field(1)%arrivaltime(region(5)%nnode))
      allocate(s%time_field(1)%time_gradient(2, region(5)%nnode))
      s%time_field(1)%arrivaltime = huge_time
      s%time_field(1)%time_gradient = 0.0_dp
   end if

   ! tf2: PsedP down — Point-source FMM in region 2 (sediment)
   if (region(2)%nnode > 0) then
      s%time_field(2)%reg => region(2)
      s%time_field(2)%vtype = VTYPE_P
      s%time_field(2)%istart => intersection(2)
      s%time_field(2)%inonstart => intersection(3)
      allocate(s%time_field(2)%arrivaltime(region(2)%nnode))
      allocate(s%time_field(2)%time_gradient(2, region(2)%nnode))
      s%time_field(2)%arrivaltime = huge_time
      s%time_field(2)%time_gradient = 0.0_dp
   end if

   ! tf3: PsedP up — reflect at sed_bot, sweep back through region 2
   if (region(2)%nnode > 0) then
      s%time_field(3)%reg => region(2)
      s%time_field(3)%vtype = VTYPE_P
      s%time_field(3)%istart => intersection(3)
      s%time_field(3)%inonstart => intersection(2)
      allocate(s%time_field(3)%arrivaltime(region(2)%nnode))
      allocate(s%time_field(3)%time_gradient(2, region(2)%nnode))
      s%time_field(3)%arrivaltime = huge_time
      s%time_field(3)%time_gradient = 0.0_dp
   end if

   ! tf4: PmP down — Point-source FMM in region 6 (seafloor->Moho)
   if (region(6)%nnode > 0) then
      s%time_field(4)%reg => region(6)
      s%time_field(4)%vtype = VTYPE_P
      s%time_field(4)%istart => intersection(2)
      s%time_field(4)%inonstart => intersection(4)
      allocate(s%time_field(4)%arrivaltime(region(6)%nnode))
      allocate(s%time_field(4)%time_gradient(2, region(6)%nnode))
      s%time_field(4)%arrivaltime = huge_time
      s%time_field(4)%time_gradient = 0.0_dp
   end if

   ! tf5: PmP up — reflect at Moho, sweep back through region 6 (seafloor->Moho)
   if (region(6)%nnode > 0) then
      s%time_field(5)%reg => region(6)
      s%time_field(5)%vtype = VTYPE_P
      s%time_field(5)%istart => intersection(4)
      s%time_field(5)%inonstart => intersection(2)
      allocate(s%time_field(5)%arrivaltime(region(6)%nnode))
      allocate(s%time_field(5)%time_gradient(2, region(6)%nnode))
      s%time_field(5)%arrivaltime = huge_time
      s%time_field(5)%time_gradient = 0.0_dp
   end if

   ! tf6: Pg_ss water-up — sweep upward through region 1 (water)
   if (region(1)%nnode > 0) then
      s%time_field(6)%reg => region(1)
      s%time_field(6)%vtype = VTYPE_P
      s%time_field(6)%istart => intersection(2)
      s%time_field(6)%inonstart => intersection(1)
      allocate(s%time_field(6)%arrivaltime(region(1)%nnode))
      allocate(s%time_field(6)%time_gradient(2, region(1)%nnode))
      s%time_field(6)%arrivaltime = huge_time
      s%time_field(6)%time_gradient = 0.0_dp
   end if

   ! tf7: Pg_ss water-down — sweep downward through region 1 (water)
   if (region(1)%nnode > 0) then
      s%time_field(7)%reg => region(1)
      s%time_field(7)%vtype = VTYPE_P
      s%time_field(7)%istart => intersection(1)
      s%time_field(7)%inonstart => intersection(2)
      allocate(s%time_field(7)%arrivaltime(region(1)%nnode))
      allocate(s%time_field(7)%time_gradient(2, region(1)%nnode))
      s%time_field(7)%arrivaltime = huge_time
      s%time_field(7)%time_gradient = 0.0_dp
   end if

   ! tf8: Pg_ss water-up 2nd — sweep upward through region 1 (water)
   if (region(1)%nnode > 0) then
      s%time_field(8)%reg => region(1)
      s%time_field(8)%vtype = VTYPE_P
      s%time_field(8)%istart => intersection(2)
      s%time_field(8)%inonstart => intersection(1)
      allocate(s%time_field(8)%arrivaltime(region(1)%nnode))
      allocate(s%time_field(8)%time_gradient(2, region(1)%nnode))
      s%time_field(8)%arrivaltime = huge_time
      s%time_field(8)%time_gradient = 0.0_dp
   end if

   ! tf9: Pg_ss2 water-up — point-source FMM in water
   if (region(1)%nnode > 0) then
      s%time_field(9)%reg => region(1)
      s%time_field(9)%vtype = VTYPE_P
      s%time_field(9)%istart => intersection(2)
      s%time_field(9)%inonstart => intersection(1)
      allocate(s%time_field(9)%arrivaltime(region(1)%nnode))
      allocate(s%time_field(9)%time_gradient(2, region(1)%nnode))
      s%time_field(9)%arrivaltime = huge_time
      s%time_field(9)%time_gradient = 0.0_dp
   end if

   ! tf10: Pg_ss2 water-down — sweep downward through water
   if (region(1)%nnode > 0) then
      s%time_field(10)%reg => region(1)
      s%time_field(10)%vtype = VTYPE_P
      s%time_field(10)%istart => intersection(1)
      s%time_field(10)%inonstart => intersection(2)
      allocate(s%time_field(10)%arrivaltime(region(1)%nnode))
      allocate(s%time_field(10)%time_gradient(2, region(1)%nnode))
      s%time_field(10)%arrivaltime = huge_time
      s%time_field(10)%time_gradient = 0.0_dp
   end if

   ! tf11: Pg_ss2 solid diving — FMM in solid (region 5)
   if (region(5)%nnode > 0) then
      s%time_field(11)%reg => region(5)
      s%time_field(11)%vtype = VTYPE_P
      s%time_field(11)%istart => intersection(2)
      allocate(s%time_field(11)%arrivaltime(region(5)%nnode))
      allocate(s%time_field(11)%time_gradient(2, region(5)%nnode))
      s%time_field(11)%arrivaltime = huge_time
      s%time_field(11)%time_gradient = 0.0_dp
   end if

   ! tf12: Pg_ss2 water-up — sweep upward through water to receiver
   if (region(1)%nnode > 0) then
      s%time_field(12)%reg => region(1)
      s%time_field(12)%vtype = VTYPE_P
      s%time_field(12)%istart => intersection(2)
      s%time_field(12)%inonstart => intersection(1)
      allocate(s%time_field(12)%arrivaltime(region(1)%nnode))
      allocate(s%time_field(12)%time_gradient(2, region(1)%nnode))
      s%time_field(12)%arrivaltime = huge_time
      s%time_field(12)%time_gradient = 0.0_dp
   end if

end do

! Assign receivers to phases based on rec_phase
! Phase mapping: 1=Pg, 2=PsedP, 3=PmP, 4=Pg_ss, 5=Pg_ss2
do n = 1, n_receivers
   iobs = obs_id(n)
   if (rec_phase(n) == 1) then
      ! Pg
      receiver(n)%n_rays = 1
      allocate(receiver(n)%ray(1))
      call ray2d_defaults(receiver(n)%ray(1))
      receiver(n)%ray(1)%source_id = iobs
      receiver(n)%ray(1)%raypath_id = 1
      receiver(n)%ray(1)%source => source(iobs)
   else if (rec_phase(n) == 2) then
      ! PsedP
      receiver(n)%n_rays = 1
      allocate(receiver(n)%ray(1))
      call ray2d_defaults(receiver(n)%ray(1))
      receiver(n)%ray(1)%source_id = iobs
      receiver(n)%ray(1)%raypath_id = 2
      receiver(n)%ray(1)%source => source(iobs)
   else if (rec_phase(n) == 3) then
      ! PmP
      receiver(n)%n_rays = 1
      allocate(receiver(n)%ray(1))
      call ray2d_defaults(receiver(n)%ray(1))
      receiver(n)%ray(1)%source_id = iobs
      receiver(n)%ray(1)%raypath_id = 3
      receiver(n)%ray(1)%source => source(iobs)
   else if (rec_phase(n) == 4) then
      ! Pg_ss
      receiver(n)%n_rays = 1
      allocate(receiver(n)%ray(1))
      call ray2d_defaults(receiver(n)%ray(1))
      receiver(n)%ray(1)%source_id = iobs
      receiver(n)%ray(1)%raypath_id = 4
      receiver(n)%ray(1)%source => source(iobs)
   else if (rec_phase(n) == 5) then
      ! Pg_ss2
      receiver(n)%n_rays = 1
      allocate(receiver(n)%ray(1))
      call ray2d_defaults(receiver(n)%ray(1))
      receiver(n)%ray(1)%source_id = iobs
      receiver(n)%ray(1)%raypath_id = 5
      receiver(n)%ray(1)%source => source(iobs)
   endif
   ! Set observed travel time from geom file
   if (receiver(n)%n_rays > 0) then
      receiver(n)%ray(1)%obs_time = geom_obs_time(n)
   end if
end do

write(*,'(A)') '  Sources and receivers initialized'
write(*,*)

call system_clock(t2, count_rate, count_max)

write(*,'(A)') '  Starting inversion...'
write(*,*)

call initialize_inversion2d(do_frechet_derivatives)

if (do_frechet_derivatives) then

   ! Only call read_observed_times2d if obs_file is specified.
   ! Otherwise, obs_time is already set from the geom file in the receiver loop above.
   if (len_trim(obs_file) > 0) then
      call read_observed_times2d(obs_file)
   end if

   ! Read pre-computed entry points from forward model (for consistent water_time)
   ! call read_entry_points_file disabled - use dynamic entry points

   call lsqr_inversion2d(niter_inv, lam_damp_in, lam_smooth_in, &
     inv_crit_chi_in, lsqr_max_iter_in, target_pert_in, &
     lam_damp_iface_in, lam_smooth_iface_in, &
     refl_weight_in, phase_wt_pmp_in, phase_wt_pgss_in, phase_wt_psed_in, corr_x_in, corr_z_in, &
     target_dd_in, iface_step_ratio_in, max_rel_pert_in, alt_inv_mode_in, corr_x_iface_in)

endif

call system_clock(t3, count_rate, count_max)

write(*,*)
write(*,'(A)') '  Inversion complete'
write(*,'(A,F12.2,A)') '  Total time: ', dble(t3-t1)/dble(count_rate), ' sec'
write(*,'(A)') '  INVERT2D finished'

contains

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

end program invert2d