module mod_inversion2d
use mod_fmt2d
use mod_rays2d
use mod_frechet2d

implicit none

integer, save :: inv_nobs = 0
real(kind=dp), allocatable, save :: inv_obs_times(:)
real(kind=dp), allocatable, save :: saved_vel(:)
logical, save :: obs_initialized = .false.
real(kind=dp), save :: inv_crit_chi = 0.1_dp

  ! Pre-computed entry points from forward model (for consistent water_time)
  logical, save :: use_entry_points_file = .false.
  integer, save :: n_entry_points = 0
  real(kind=dp), allocatable, save :: entry_point_data(:,:)  ! (5, n_entry_points)
  ! col 1=ircv, 2=src_id, 3=path_id, 4=entry_x, 5=entry_z, 6=water_time

! Pre-computed Gaussian smoothing sparse matrices (CSR format)
! R_h: horizontal smoothing, R_v: vertical smoothing
! R = R_v * R_h (applied as separable)
logical, save :: gs_built = .false.
integer, save :: gs_nparm = 0
integer, save :: gs_h_nnz = 0
integer, save :: gs_v_nnz = 0
integer, allocatable, save :: gs_h_rowptr(:)
integer, allocatable, save :: gs_h_colind(:)
real(kind=dp), allocatable, save :: gs_h_values(:)
integer, allocatable, save :: gs_v_rowptr(:)
integer, allocatable, save :: gs_v_colind(:)
real(kind=dp), allocatable, save :: gs_v_values(:)

! When .true., run_forward2d skips the Frechet derivative computation for each
! ray. Used by line-search RMS checks, which only need receiver travel times.
logical, save :: skip_frechet_in_forward = .false.

! Pre-computed interface Gaussian smoothing sparse matrix (1D horizontal, R_h only)
! Used by smooth_type=3 (Hybrid): velocity uses Laplacian, interface uses Gaussian R_h
logical, save :: gsi_built = .false.
integer, save :: gsi_nparm = 0
integer, save :: gsi_h_nnz = 0
integer, allocatable, save :: gsi_h_rowptr(:)
integer, allocatable, save :: gsi_h_colind(:)
real(kind=dp), allocatable, save :: gsi_h_values(:)

interface
   subroutine velocities_to_grid2d(pgrid)
     use type_definitions2d
     type(Tpropagation_grid2d),target :: pgrid
   end subroutine
end interface

interface
   subroutine velocities_to_intersection2d(isec)
     use type_definitions2d
     type(Tintersection2d) :: isec
   end subroutine
end interface

interface
   subroutine velocities_to_region2d(reg,pgrid)
     use type_definitions2d
     type(Tregion2d) :: reg
     type(Tpropagation_grid2d),target :: pgrid
   end subroutine
end interface

interface
   subroutine transfer_times_to_intersection(isec, tf)
     use type_definitions2d
     type(Tintersection2d) :: isec
     type(Ttime_field2d) :: tf
   end subroutine
end interface

interface
   subroutine reflect_gradient_at_interface(isec, target_reg, target_vtype)
     use type_definitions2d
     type(Tintersection2d) :: isec
     type(Tregion2d) :: target_reg
     integer :: target_vtype
   end subroutine
end interface

interface
   subroutine fit_gradient_at_node2d(reg, centernode, dtdx, dtdz)
     use type_definitions2d
     type(Tregion2d) :: reg
     integer :: centernode
     real(kind=dp) :: dtdx, dtdz
   end subroutine
end interface

interface
   subroutine tag_regular_nodes2d(pgrid)
     use type_definitions2d
     type(Tpropagation_grid2d) :: pgrid
   end subroutine
end interface

private
public :: lsqr_inversion2d, update_model2d, compute_rms2d, write_model2d, &
          read_frechet_data2d, read_observed_times2d, &
          sparse_Ax, sparse_ATx, cgnr_solver, &
          compute_path_integral_time2d, &
          write_residuals2d, write_ray_density2d, &
          read_entry_points_file, get_entry_point

contains

subroutine lsqr_inversion2d(niter, lam_damp_in, lam_smooth_in, &
  inv_crit_chi_in, lsqr_max_iter_in, target_pert_in, &
  lam_damp_iface_in, lam_smooth_iface_in, &
  refl_weight_in, phase_wt_pmp_in, phase_wt_pgss_in, phase_wt_psed_in, corr_x_in, corr_z_in, &
  target_dd_in, iface_step_ratio_in, max_rel_pert_in, alt_inv_mode_in, corr_x_iface_in)
  implicit none
  integer, intent(in) :: niter, lsqr_max_iter_in
  real(kind=dp), intent(in) :: lam_damp_in, lam_smooth_in, inv_crit_chi_in
  real(kind=dp), intent(in) :: target_pert_in
  real(kind=dp), intent(in) :: lam_damp_iface_in, lam_smooth_iface_in
  real(kind=dp), intent(in) :: refl_weight_in, phase_wt_pmp_in, phase_wt_pgss_in, phase_wt_psed_in
  real(kind=dp), intent(in) :: corr_x_in, corr_z_in
  real(kind=dp), intent(in) :: target_dd_in
  real(kind=dp), intent(in) :: iface_step_ratio_in, max_rel_pert_in
  integer, intent(in) :: alt_inv_mode_in
  real(kind=dp), intent(in) :: corr_x_iface_in

  integer :: iter, ndata, nparm, nnz, ndata_total, nnz_data
  integer :: alt_inv_mode  ! 0=joint, 1=vel-only, 2=iface-only
  integer :: ircv, iray, k, i, j, idata, n_outliers, n_active
  integer, allocatable :: gindx(:), row_ptr(:)
  real(kind=dp), allocatable :: gval(:), residual(:)
  real(kind=dp), allocatable :: dmodel_total(:)
  real(kind=dp), allocatable :: dmodel_total_sum(:)
  real(kind=dp) :: rms_orig, prev_rms, rms_new
  real(kind=dp) :: res_check
  real(kind=dp), allocatable :: mvscale(:)
  real(kind=dp), allocatable :: mvscale_saved(:)
  real(kind=dp), allocatable :: data_wt(:)
  real(kind=dp) :: avg_dmv, max_dv
  real(kind=dp) :: avg_dmd, target_dd, wd_best_iface
  real(kind=dp) :: wd_best, wdv_iface_val

  ! Best model tracking
  real(kind=dp) :: best_rms
  integer :: best_iter
  integer :: ls_prev_idx   ! ladder index accepted in the previous iteration
  real(kind=dp), allocatable :: best_vel(:,:)
  real(kind=dp), allocatable :: best_iface_z(:)
  logical :: use_best_model

  integer :: ig, ix, iz, idx, vtype, nx_vg, nz_vg
  integer :: ii, iiface, ip
  integer :: n_vel_parms, n_iface_parms
  type(Tvelocity_grid2d), pointer :: vg
  type(Tinterface2d), pointer :: ifc

  logical, allocatable :: node_hit(:)
  integer :: n_hit

  real(kind=dp), allocatable :: cm(:), ecmi(:)

  nparm = n_inv_parms
  if (nparm <= 0) return

  prev_rms = huge_time
  best_rms = huge_time
  best_iter = 0
  ls_prev_idx = 1
  use_best_model = .false.

  inv_crit_chi = inv_crit_chi_in

  allocate(mvscale(nparm))
  allocate(dmodel_total(nparm))
  allocate(dmodel_total_sum(nparm))
  dmodel_total = 0.0_dp
  dmodel_total_sum = 0.0_dp
  mvscale = 1.0_dp

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx; nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0) then
          mvscale(idx) = 1.0_dp / max(vg%velocity_init(vtype, ix, iz), 1.0_dp)
        end if
      end do
    end do
  end do
  do ii = 1, n_inv_iface
    iiface = ifaces_to_be_inv(ii)
    ifc => intrface(iiface)
    do ix = 1, ifc%nx
      idx = ifc%start_index + ix
      if (idx >= 1 .and. idx <= nparm) then
        mvscale(idx) = refl_weight_in
      end if
    end do
  end do

  n_vel_parms = 0
  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    nx_vg = vg%nx; nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        if (vg%inv_idx(ix, iz) > 0) n_vel_parms = n_vel_parms + 1
      end do
    end do
  end do
  n_iface_parms = nparm - n_vel_parms

  allocate(cm(nparm), ecmi(nparm))
  cm = 1.0_dp
  ecmi = 0.0_dp

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx; nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0) then
          cm(idx) = (0.05_dp * vg%velocity_init(vtype, ix, iz))**2
          ecmi(idx) = lam_damp_in / cm(idx)
        end if
      end do
    end do
  end do
  do ii = 1, n_inv_iface
    iiface = ifaces_to_be_inv(ii)
    ifc => intrface(iiface)
    do ix = 1, ifc%nx
      idx = ifc%start_index + ix
      if (idx >= 1 .and. idx <= nparm) then
        cm(idx) = 1.0_dp
        ecmi(idx) = lam_damp_iface_in / cm(idx)
      end if
    end do
  end do

  ! Open inversion log file
  open(unit=99, file='invert_log_latest.txt', status='replace', action='write')

  ! Write initial residuals and ray density (before any inversion)
  call run_forward2d()
  call write_residuals2d('residuals_initial.dat')
  call write_ray_density2d('ray_density_initial.dat')

  do iter = 1, niter

    ! Alternating inversion mode:
    ! 0 = joint (velocity + interface)
    ! 1 = velocity only (interface fixed)
    ! 2 = interface only (velocity fixed)
    ! Alternating strategy: 3 velocity-only iterations per 1 interface-only iteration
    if (n_inv_iface > 0) then
      if (alt_inv_mode_in == -1) then
        ! Auto: alternating mode (3 vel + 1 iface)
        if (mod(iter, 4) == 0) then
          alt_inv_mode = 2  ! every 4th iteration: interface-only update
        else
          alt_inv_mode = 1  ! other iterations: velocity-only update
        end if
      else
        ! User-specified mode: 0=joint, 1=vel-only, 2=iface-only
        alt_inv_mode = alt_inv_mode_in
      end if
    else
      alt_inv_mode = 0  ! no interface, always velocity-only
    end if

    write(*,'(A,I4,A,I1)') '  Iteration ', iter, ' alt_inv_mode=', alt_inv_mode
    if (alt_inv_mode == 1) write(*,'(A)') '    -> Velocity-only update'
    if (alt_inv_mode == 2) write(*,'(A)') '    -> Interface-only update'

    do ircv = 1, n_receivers
      do iray = 1, receiver(ircv)%n_rays
        receiver(ircv)%ray(iray)%valid = .true.
      end do
    end do

    call run_forward2d()

    ndata = 0
    nnz = 0
    n_outliers = 0
    do ircv = 1, n_receivers
      do iray = 1, receiver(ircv)%n_rays
        if (receiver(ircv)%ray(iray)%valid .and. &
            receiver(ircv)%ray(iray)%obs_time >= 0.0_dp) then
          res_check = receiver(ircv)%ray(iray)%obs_time - &
                      receiver(ircv)%ray(iray)%receiver_time
          if (abs(res_check) > inv_crit_chi) then
            n_outliers = n_outliers + 1
            cycle
          end if
          ndata = ndata + 1
          nnz = nnz + receiver(ircv)%ray(iray)%n_pdev
        end if
      end do
    end do

    write(*,'(A,I6,A,I4,A,F6.2)') '  ndata=', ndata, ' n_outliers=', n_outliers, &
      ' crit_chi=', inv_crit_chi

    if (ndata == 0) then
      write(*,'(A,I4,A)') 'Iteration ', iter, ' no valid data'
      exit
    end if

    allocate(gindx(nnz), gval(nnz), residual(ndata), row_ptr(ndata+1))
    allocate(data_wt(ndata))

    call read_frechet_data2d(ndata_total, nnz_data, gindx, gval, residual, row_ptr)

    ! Debug: count interface Frechet derivatives
    if (n_inv_iface > 0) then
      n_active = 0
      do j = n_vel_parms + 1, nparm
        if (any(gindx == j)) n_active = n_active + 1
      end do
      write(*,'(A,I6,A,I6)') '  Interface params with Frechet: ', n_active, ' / ', n_iface_parms
    end if

    rms_orig = compute_rms2d(ndata, residual)

    idata = 0
    do ircv = 1, n_receivers
      do iray = 1, receiver(ircv)%n_rays
        if (.not. receiver(ircv)%ray(iray)%valid) cycle
        if (receiver(ircv)%ray(iray)%obs_time < 0.0_dp) cycle
        res_check = receiver(ircv)%ray(iray)%obs_time - receiver(ircv)%ray(iray)%receiver_time
        if (abs(res_check) > inv_crit_chi) cycle

        idata = idata + 1
        data_wt(idata) = 1.0_dp
        if (receiver(ircv)%ray(iray)%raypath_id == 2) then  ! PsedP weight
          data_wt(idata) = data_wt(idata) * phase_wt_psed_in
        else if (receiver(ircv)%ray(iray)%raypath_id == 3) then  ! PmP weight
          data_wt(idata) = data_wt(idata) * phase_wt_pmp_in
        else if (receiver(ircv)%ray(iray)%raypath_id == 4 .or. &
                 receiver(ircv)%ray(iray)%raypath_id == 5) then
          data_wt(idata) = data_wt(idata) * phase_wt_pgss_in  ! Pg_ss/Pg_ss2 weight
        end if
      end do
    end do

    ! Apply alternating inversion: zero out mvscale for fixed parameters
    ! This effectively removes them from the LSQR solve
    allocate(mvscale_saved(nparm))
    mvscale_saved(1:nparm) = mvscale(1:nparm)
    if (alt_inv_mode == 1) then
      ! Velocity-only: fix interface parameters
      do j = n_vel_parms + 1, nparm
        mvscale(j) = 0.0_dp
      end do
    else if (alt_inv_mode == 2) then
      ! Interface-only: fix velocity parameters
      do j = 1, n_vel_parms
        mvscale(j) = 0.0_dp
      end do
    end if

    do i = 1, ndata
      do k = row_ptr(i), row_ptr(i+1) - 1
        j = gindx(k)
        if (j >= 1 .and. j <= nparm) then
          gval(k) = gval(k) * mvscale(j) * data_wt(i)
        end if
      end do
      residual(i) = residual(i) * data_wt(i)
    end do

    write(*,'(A,F14.8)') '    RMS_orig=', rms_orig

    allocate(node_hit(nparm))
    node_hit = .false.
    n_hit = 0
    do i = 1, ndata
      do k = row_ptr(i), row_ptr(i+1) - 1
        j = gindx(k)
        if (j >= 1 .and. j <= nparm) then
          if (.not. node_hit(j)) then
            node_hit(j) = .true.
            n_hit = n_hit + 1
          end if
        end if
      end do
    end do
    write(*,'(A,I6,A,I6)') '    Ray-covered nodes: ', n_hit, ' / ', nparm

    if (iter <= 3) then
      call fd_verify_frechet()
    end if

    block
      ! ============================================================
      ! LSQR solver with auto_damping (tomo2d approach)
      ! Key change: instead of scaling dm to match target_pert,
      ! we search for the right damping weight so that LSQR
      ! naturally produces avg_dmv ≈ target_dv
      ! ============================================================
      real(kind=dp) :: w_smooth, w_smooth_iface
      real(kind=dp), allocatable :: gtg_diag(:)
      real(kind=dp), allocatable :: smooth_work(:), smooth_out(:)
      real(kind=dp) :: target_dv
      real(kind=dp) :: gtg_min, gtg_max, gtg_med
      real(kind=dp), allocatable :: gtg_sorted(:)
      real(kind=dp) :: eps_min

      ! LSQR variables
      integer :: ndata_aug, lsqr_iter, max_lsqr_iter
      real(kind=dp), allocatable :: u(:), v_vec(:), w_vec(:), z_vec(:)
      real(kind=dp), allocatable :: Av_aug(:), Atu_aug(:), Av_data(:)
      real(kind=dp) :: rnorm_ls, bnorm_ls, test2_ls

      ! Model update variables
      real(kind=dp) :: trial_s, trial_v, cur_s_v, init_v
      real(kind=dp) :: step_len, ls_step, rms_try, rms_best
      integer :: ls_idx, ls_idx_ok
      real(kind=dp), allocatable :: saved_vel(:,:)
      real(kind=dp), allocatable :: saved_iface_z(:)
      real(kind=dp) :: max_rel_pert
      integer :: ip, ii, iiface
      integer :: max_nx_vg, max_nz_vg
      real(kind=dp) :: dz
      real(kind=dp) :: sum_dz2, step_iface
      real(kind=dp) :: rms_after_iface
      integer :: n_dz_nonzero, ix_mid

      ! Auto-damping variables
      real(kind=dp) :: wd1, wd2, wd_rts, wd_xl
      real(kind=dp) :: fl, f_val, dx_wd, swap_wd
      real(kind=dp) :: ave_dmv0, ddm
      integer :: j_sec, secant_itermax
      real(kind=dp), allocatable :: dm_trial(:)

      max_rel_pert = max_rel_pert_in
      w_smooth = lam_smooth_in
      w_smooth_iface = lam_smooth_iface_in
      max_lsqr_iter = min(lsqr_max_iter_in, 2000)

      ! Recompute mvscale from velocity_init (not current velocity)
      do ig = 1, n_inv_vgrid
        vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
        vtype = vgrids_to_be_inv(ig)%vtype
        nx_vg = vg%nx; nz_vg = vg%nz
        do iz = 1, nz_vg
          do ix = 1, nx_vg
            idx = vg%inv_idx(ix, iz)
            if (idx > 0) then
              mvscale(idx) = 1.0_dp / max(vg%velocity_init(vtype, ix, iz), 1.0_dp)
            end if
          end do
        end do
      end do

      ! Compute gtg diagonal for active node detection
      allocate(gtg_diag(nparm))
      gtg_diag = 0.0_dp
      do i = 1, ndata
        do k = row_ptr(i), row_ptr(i+1) - 1
          j = gindx(k)
          if (j >= 1 .and. j <= nparm) then
            gtg_diag(j) = gtg_diag(j) + gval(k) * gval(k)
          end if
        end do
      end do

      allocate(gtg_sorted(nparm))
      gtg_sorted = gtg_diag
      call sort_array(gtg_sorted, nparm)
      gtg_min = minval(gtg_diag, mask=gtg_diag > 0.0_dp)
      gtg_max = maxval(gtg_diag)
      gtg_med = gtg_sorted(nparm/2)
      write(*,'(A,ES10.3,A,ES10.3,A,ES10.3)') '    gtg_diag: min=', gtg_min, ' max=', gtg_max, ' med=', gtg_med
      deallocate(gtg_sorted)

      eps_min = max(1.0e-6_dp, 1.0e-4_dp * gtg_max)

      n_active = 0
      do j = 1, nparm
        if (j <= n_vel_parms .and. gtg_diag(j) > eps_min * 0.01_dp) then
          n_active = n_active + 1
        end if
      end do
      write(*,'(A,I6,A,I6)') '    Active vel nodes: ', n_active, ' / ', n_vel_parms

      target_dv = target_pert_in
      if (target_dv <= 0.0_dp) target_dv = 0.03_dp

      target_dd = target_dd_in
      if (target_dd <= 0.0_dp) target_dd = 0.5_dp  ! default 0.5 km

      allocate(smooth_work(nparm), smooth_out(nparm))
      ! Augmented system layout depends on smooth_type:
      ! smooth_type=1 (Laplacian): [G; w*L; wdv*I] => ndata_aug = ndata + nparm + nparm
      ! smooth_type=2 (Gaussian):  [G; w*R_h; w*R_v; wdv*I] => ndata_aug = ndata + 2*nparm + nparm
      ! smooth_type=3 (Hybrid):    [G; w*L_vel + w_iface*R_h_iface; wdv*I] => ndata_aug = ndata + nparm + nparm
      if (smooth_type == 1 .or. smooth_type == 3) then
        ndata_aug = ndata + nparm + nparm
      else
        ndata_aug = ndata + 2*nparm + nparm
      end if
      allocate(u(ndata_aug), v_vec(nparm), w_vec(nparm), z_vec(nparm))
      allocate(Av_aug(ndata_aug), Atu_aug(nparm), Av_data(ndata))
      allocate(dm_trial(nparm))

      ! Pre-compute Gaussian smoothing sparse matrices (once for all iterations)
      if (w_smooth > 0.0_dp .and. smooth_type == 2) then
        call build_gauss_smooth_matrix(nparm, corr_x_in, corr_z_in, mvscale)
      end if
      if (w_smooth_iface > 0.0_dp .and. (smooth_type == 2 .or. smooth_type == 3)) then
        call build_iface_gauss_smooth_matrix(nparm, n_vel_parms, corr_x_iface_in, mvscale)
      end if

      ! ============================================================
      ! Damping approach:
      ! If user specified lam_damp > 0, use fixed damping.
      ! Otherwise, use auto-damping (following tomo2d):
      ! 1. Solve without damping first
      ! 2. If avg_dmv <= target_dv, no damping needed
      ! 3. Otherwise, search for damping weight that gives avg_dmv ≈ target_dv
      ! ============================================================

      if (lam_damp_in > 0.0_dp) then
        ! User-specified fixed damping
        if (n_inv_iface > 0 .and. lam_damp_iface_in > 0.0_dp) then
          wdv_iface_val = lam_damp_iface_in
          wd_best_iface = lam_damp_iface_in
        else
          wdv_iface_val = lam_damp_in
          wd_best_iface = lam_damp_in
        end if
        call run_lsqr(lam_damp_in, wdv_iface_val, w_smooth, w_smooth_iface, ndata, nparm, nnz, gindx, gval, row_ptr, &
          residual, gtg_diag, eps_min, n_vel_parms, ndata_aug, max_lsqr_iter, &
          mvscale, u, v_vec, w_vec, z_vec, Av_aug, Atu_aug, Av_data, &
          smooth_work, smooth_out, dmodel_total, dmodel_total_sum, corr_x_in, corr_z_in, &
          smooth_type, lsqr_iter, rnorm_ls, bnorm_ls, test2_ls)
        dm_trial = dmodel_total + dmodel_total_sum
        avg_dmv = calc_avg_dmv(dm_trial, gtg_diag, eps_min, nparm, n_vel_parms)
        if (n_iface_parms > 0) then
          avg_dmd = calc_avg_dmd(dm_trial, mvscale, nparm, n_vel_parms)
        else
          avg_dmd = 0.0_dp
        end if
        wd_best = lam_damp_in
        write(*,'(A,F8.3,A,F8.3,A,F8.3,A)') '    Fixed damping: wd=', lam_damp_in, &
          ' avg_dmv=', avg_dmv*100.0_dp, '% avg_dmd=', avg_dmd, 'km'
      else
      ! Auto-damping (following tomo2d)
      ! Step 1: Solve without any damping
      call run_lsqr(0.0_dp, 0.0_dp, w_smooth, w_smooth_iface, ndata, nparm, nnz, gindx, gval, row_ptr, &
        residual, gtg_diag, eps_min, n_vel_parms, ndata_aug, max_lsqr_iter, &
        mvscale, u, v_vec, w_vec, z_vec, Av_aug, Atu_aug, Av_data, &
        smooth_work, smooth_out, dmodel_total, dmodel_total_sum, corr_x_in, corr_z_in, &
        smooth_type, lsqr_iter, rnorm_ls, bnorm_ls, test2_ls)

      dm_trial = dmodel_total + dmodel_total_sum
      avg_dmv = calc_avg_dmv(dm_trial, gtg_diag, eps_min, nparm, n_vel_parms)
      if (n_iface_parms > 0) then
        avg_dmd = calc_avg_dmd(dm_trial, mvscale, nparm, n_vel_parms)
      else
        avg_dmd = 0.0_dp
      end if
      write(*,'(A,F8.3,A,F8.3,A)') '    No-damping: avg_dmv=', avg_dmv*100.0_dp, &
        '% avg_dmd=', avg_dmd, 'km'

      ! Step 2: Auto-damping for interface (depth) parameters
      wd_best_iface = 0.0_dp
      if (n_iface_parms > 0 .and. avg_dmd > target_dd) then
        write(*,'(A)') '    Searching for interface damping...'
        secant_itermax = 10
        ddm = target_dd * 0.3_dp

        wd1 = 1.0_dp
        wd2 = 100.0_dp

        ! Evaluate at wd1
        call run_lsqr(0.0_dp, wd1, w_smooth, w_smooth_iface, ndata, nparm, nnz, gindx, gval, row_ptr, &
          residual, gtg_diag, eps_min, n_vel_parms, ndata_aug, max_lsqr_iter, &
          mvscale, u, v_vec, w_vec, z_vec, Av_aug, Atu_aug, Av_data, &
          smooth_work, smooth_out, dm_trial, dmodel_total_sum, corr_x_in, corr_z_in, &
          smooth_type, lsqr_iter, rnorm_ls, bnorm_ls, test2_ls)
        dm_trial = dm_trial + dmodel_total_sum
        fl = calc_avg_dmd(dm_trial, mvscale, nparm, n_vel_parms) - target_dd
        write(*,'(A,F8.3,A,F8.3,A)') '    wdd1=', wd1, ' avg_dmd=', (fl+target_dd), 'km'

        ! Evaluate at wd2
        call run_lsqr(0.0_dp, wd2, w_smooth, w_smooth_iface, ndata, nparm, nnz, gindx, gval, row_ptr, &
          residual, gtg_diag, eps_min, n_vel_parms, ndata_aug, max_lsqr_iter, &
          mvscale, u, v_vec, w_vec, z_vec, Av_aug, Atu_aug, Av_data, &
          smooth_work, smooth_out, dm_trial, dmodel_total_sum, corr_x_in, corr_z_in, &
          smooth_type, lsqr_iter, rnorm_ls, bnorm_ls, test2_ls)
        dm_trial = dm_trial + dmodel_total_sum
        f_val = calc_avg_dmd(dm_trial, mvscale, nparm, n_vel_parms) - target_dd
        write(*,'(A,F8.3,A,F8.3,A)') '    wdd2=', wd2, ' avg_dmd=', (f_val+target_dd), 'km'

        ! Ensure wd1 gives smaller |f| than wd2 (swap if needed)
        if (abs(fl) < abs(f_val)) then
          wd_rts = wd1
          wd_xl = wd2
          swap_wd = fl; fl = f_val; f_val = swap_wd
        else
          wd_xl = wd1
          wd_rts = wd2
        end if

        ! Secant/bisection iteration
        do j_sec = 1, secant_itermax
          dx_wd = (wd_xl - wd_rts) * f_val / (f_val - fl)
          if (wd_rts + dx_wd <= 0.0_dp) then
            wd_xl = wd_rts; fl = f_val
            wd_rts = wd_rts * 0.5_dp
          else if (wd_rts + dx_wd > 1.0e7_dp) then
            wd_xl = wd_rts; fl = f_val
            wd_rts = wd_rts + (1.0e7_dp - wd_rts) * 0.5_dp
          else
            wd_xl = wd_rts; fl = f_val
            wd_rts = wd_rts + dx_wd
          end if

          call run_lsqr(0.0_dp, wd_rts, w_smooth, w_smooth_iface, ndata, nparm, nnz, gindx, gval, row_ptr, &
            residual, gtg_diag, eps_min, n_vel_parms, ndata_aug, max_lsqr_iter, &
            mvscale, u, v_vec, w_vec, z_vec, Av_aug, Atu_aug, Av_data, &
            smooth_work, smooth_out, dm_trial, dmodel_total_sum, corr_x_in, corr_z_in, &
            smooth_type, lsqr_iter, rnorm_ls, bnorm_ls, test2_ls)
          dm_trial = dm_trial + dmodel_total_sum
          f_val = calc_avg_dmd(dm_trial, mvscale, nparm, n_vel_parms) - target_dd
          write(*,'(A,I2,A,F10.2,A,F8.3,A)') '    secant_d ', j_sec, &
            ' wdd=', wd_rts, ' avg_dmd=', (f_val+target_dd), 'km'

          if (wd_rts > 1.0e7_dp) exit
          if (abs(f_val) < ddm) exit
        end do

        wd_best_iface = wd_rts
        write(*,'(A,F10.2,A)') '    Interface damping: wdd=', wd_best_iface
      else if (n_iface_parms > 0) then
        write(*,'(A)') '    No interface damping needed (avg_dmd <= target_dd)'
      end if

      ! Step 3: Auto-damping for velocity parameters (with wdv_iface=wd_best_iface)
      if (avg_dmv <= target_dv) then
        ! No velocity damping needed
        wd_best = 0.0_dp
        write(*,'(A)') '    No velocity damping needed (avg_dmv <= target)'
      else
        ! Search for appropriate velocity damping weight
        write(*,'(A)') '    Searching for velocity damping...'
        secant_itermax = 10
        ddm = target_dv * 0.3_dp

        wd1 = 1.0_dp
        wd2 = 100.0_dp

        ! Evaluate at wd1
        call run_lsqr(wd1, wd_best_iface, w_smooth, w_smooth_iface, ndata, nparm, nnz, gindx, gval, row_ptr, &
          residual, gtg_diag, eps_min, n_vel_parms, ndata_aug, max_lsqr_iter, &
          mvscale, u, v_vec, w_vec, z_vec, Av_aug, Atu_aug, Av_data, &
          smooth_work, smooth_out, dm_trial, dmodel_total_sum, corr_x_in, corr_z_in, &
          smooth_type, lsqr_iter, rnorm_ls, bnorm_ls, test2_ls)
        dm_trial = dm_trial + dmodel_total_sum
        fl = calc_avg_dmv(dm_trial, gtg_diag, eps_min, nparm, n_vel_parms) - target_dv
        write(*,'(A,F8.3,A,F8.3,A)') '    wdv1=', wd1, ' avg_dmv=', (fl+target_dv)*100.0_dp, '%'

        ! Evaluate at wd2
        call run_lsqr(wd2, wd_best_iface, w_smooth, w_smooth_iface, ndata, nparm, nnz, gindx, gval, row_ptr, &
          residual, gtg_diag, eps_min, n_vel_parms, ndata_aug, max_lsqr_iter, &
          mvscale, u, v_vec, w_vec, z_vec, Av_aug, Atu_aug, Av_data, &
          smooth_work, smooth_out, dm_trial, dmodel_total_sum, corr_x_in, corr_z_in, &
          smooth_type, lsqr_iter, rnorm_ls, bnorm_ls, test2_ls)
        dm_trial = dm_trial + dmodel_total_sum
        f_val = calc_avg_dmv(dm_trial, gtg_diag, eps_min, nparm, n_vel_parms) - target_dv
        write(*,'(A,F8.3,A,F8.3,A)') '    wdv2=', wd2, ' avg_dmv=', (f_val+target_dv)*100.0_dp, '%'

        ! Ensure wd1 gives smaller |f| than wd2 (swap if needed)
        if (abs(fl) < abs(f_val)) then
          wd_rts = wd1
          wd_xl = wd2
          swap_wd = fl; fl = f_val; f_val = swap_wd
        else
          wd_xl = wd1
          wd_rts = wd2
        end if

        ! Secant/bisection iteration
        do j_sec = 1, secant_itermax
          dx_wd = (wd_xl - wd_rts) * f_val / (f_val - fl)
          if (wd_rts + dx_wd <= 0.0_dp) then
            wd_xl = wd_rts; fl = f_val
            wd_rts = wd_rts * 0.5_dp
          else if (wd_rts + dx_wd > 1.0e7_dp) then
            wd_xl = wd_rts; fl = f_val
            wd_rts = wd_rts + (1.0e7_dp - wd_rts) * 0.5_dp
          else
            wd_xl = wd_rts; fl = f_val
            wd_rts = wd_rts + dx_wd
          end if

          call run_lsqr(wd_rts, wd_best_iface, w_smooth, w_smooth_iface, ndata, nparm, nnz, gindx, gval, row_ptr, &
            residual, gtg_diag, eps_min, n_vel_parms, ndata_aug, max_lsqr_iter, &
            mvscale, u, v_vec, w_vec, z_vec, Av_aug, Atu_aug, Av_data, &
            smooth_work, smooth_out, dm_trial, dmodel_total_sum, corr_x_in, corr_z_in, &
            smooth_type, lsqr_iter, rnorm_ls, bnorm_ls, test2_ls)
          dm_trial = dm_trial + dmodel_total_sum
          f_val = calc_avg_dmv(dm_trial, gtg_diag, eps_min, nparm, n_vel_parms) - target_dv
          write(*,'(A,I2,A,F10.2,A,F8.3,A)') '    secant_v ', j_sec, &
            ' wdv=', wd_rts, ' avg_dmv=', (f_val+target_dv)*100.0_dp, '%'

          if (wd_rts > 1.0e7_dp) exit
          if (abs(f_val) < ddm) exit
        end do

        wd_best = wd_rts
      end if  ! end velocity auto-damping

      ! Final solve with best damping weights
      wdv_iface_val = wd_best_iface
      call run_lsqr(wd_best, wdv_iface_val, w_smooth, w_smooth_iface, ndata, nparm, nnz, gindx, gval, row_ptr, &
        residual, gtg_diag, eps_min, n_vel_parms, ndata_aug, max_lsqr_iter, &
        mvscale, u, v_vec, w_vec, z_vec, Av_aug, Atu_aug, Av_data, &
        smooth_work, smooth_out, dmodel_total, dmodel_total_sum, corr_x_in, corr_z_in, &
        smooth_type, lsqr_iter, rnorm_ls, bnorm_ls, test2_ls)

      avg_dmv = calc_avg_dmv(dmodel_total + dmodel_total_sum, gtg_diag, eps_min, nparm, n_vel_parms)
      if (n_iface_parms > 0) then
        avg_dmd = calc_avg_dmd(dmodel_total + dmodel_total_sum, mvscale, nparm, n_vel_parms)
      else
        avg_dmd = 0.0_dp
      end if
      end if  ! end if (lam_damp > 0) ... else ...

      write(*,'(A,ES10.3,A,ES10.3,A,F8.3,A,A,F8.3,A)') '    Final: wdv=', wd_best, &
        ' wdd=', wd_best_iface, ' avg_dmv=', avg_dmv*100.0_dp, '%', ' avg_dmd=', avg_dmd, 'km'

      ! Zero out inactive nodes
      do j = 1, nparm
        if (j <= n_vel_parms .and. gtg_diag(j) < eps_min * 0.01_dp) then
          dmodel_total(j) = 0.0_dp
        end if
      end do

      ! Restore mvscale and zero out dmodel_total for fixed parameters (alternating inversion)
      mvscale(1:nparm) = mvscale_saved(1:nparm)
      deallocate(mvscale_saved)
      if (alt_inv_mode == 1) then
        ! Velocity-only: zero out interface update
        do j = n_vel_parms + 1, nparm
          dmodel_total(j) = 0.0_dp
        end do
      else if (alt_inv_mode == 2) then
        ! Interface-only: zero out velocity update
        do j = 1, n_vel_parms
          dmodel_total(j) = 0.0_dp
        end do
      end if

      ! Apply dm_cap (now 30% instead of 10%)
      do i = 1, nparm
        if (i <= n_vel_parms) then
          if (abs(dmodel_total(i)) > max_rel_pert) then
            dmodel_total(i) = sign(max_rel_pert, dmodel_total(i))
          end if
        end if
      end do

      ! ============================================================
      ! Model update with simple line search (step 1.0 or 0.5)
      ! Like tomo2d but with minimal line search for stability
      ! ============================================================
      ! Compute max grid dimensions for saved_vel allocation
      max_nx_vg = 0
      max_nz_vg = 0
      do ig = 1, n_inv_vgrid
        vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
        if (vg%nx > max_nx_vg) max_nx_vg = vg%nx
        if (vg%nz > max_nz_vg) max_nz_vg = vg%nz
      end do
      if (.not. allocated(saved_vel)) then
        allocate(saved_vel(max_nx_vg, max_nz_vg))
      else if (size(saved_vel,1) < max_nx_vg .or. size(saved_vel,2) < max_nz_vg) then
        deallocate(saved_vel)
        allocate(saved_vel(max_nx_vg, max_nz_vg))
      end if

      ! Save current velocity
      do ig = 1, n_inv_vgrid
        vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
        vtype = vgrids_to_be_inv(ig)%vtype
        nx_vg = vg%nx; nz_vg = vg%nz
        do iz = 1, nz_vg
          do ix = 1, nx_vg
            saved_vel(ix, iz) = vg%velocity(vtype, ix, iz)
          end do
        end do
      end do

      ! Save current interface depths
      if (n_inv_iface > 0) then
        allocate(saved_iface_z(n_iface_parms))
        ip = 0
        do ii = 1, n_inv_iface
          iiface = ifaces_to_be_inv(ii)
          ifc => intrface(iiface)
          do ix = 1, ifc%nx
            ip = ip + 1
            saved_iface_z(ip) = ifc%z_nodes(ix)
          end do
        end do
      end if

      ! Line search with warm start: begin at the step accepted in the
      ! previous iteration and walk down the ladder (1.0, 0.75, 0.5, 0.25,
      ! 0.125, 0.0625) until a step lowers RMS below rms_orig. Smaller steps
      ! are never retried with a larger one, so once the ladder has moved
      ! down it stays there. rms_new = rms_orig means every remaining step
      ! failed.
      step_len = 0.0_dp
      rms_new = rms_orig
      rms_best = rms_orig
      ls_idx_ok = 1

      do ls_idx = ls_prev_idx, 6
        if (ls_idx == 1) then
          ls_step = 1.0_dp
        else if (ls_idx == 2) then
          ls_step = 0.75_dp
        else if (ls_idx == 3) then
          ls_step = 0.5_dp
        else if (ls_idx == 4) then
          ls_step = 0.25_dp
        else if (ls_idx == 5) then
          ls_step = 0.125_dp
        else
          ls_step = 0.0625_dp
        end if
        call restore_model_state(saved_vel, saved_iface_z, n_iface_parms)
        call apply_model_step(ls_step, ls_step*iface_step_ratio_in, dmodel_total, mvscale, nparm, n_vel_parms, &
          saved_iface_z, n_iface_parms)
        call compute_rms_after_forward(rms_try)
        write(*,'(A,F6.3,A,F10.6)') '    step=', ls_step, ' RMS=', rms_try
        if (rms_try < rms_best) then
          rms_best = rms_try
          step_len = ls_step
          ls_idx_ok = ls_idx
        end if
        ! Early exit as soon as a step reduces RMS.
        if (rms_try < rms_orig) exit
      end do

      ! Remember the accepted step (0.0 keeps the previous value)
      if (step_len > 0.0_dp) ls_prev_idx = ls_idx_ok

      ! Apply the best step found
      if (step_len > 0.0_dp) then
        call restore_model_state(saved_vel, saved_iface_z, n_iface_parms)
        call apply_model_step(step_len, step_len*iface_step_ratio_in, dmodel_total, mvscale, nparm, n_vel_parms, &
          saved_iface_z, n_iface_parms)
        rms_new = rms_best
      else
        ! No step reduced RMS, revert
        call restore_model_state(saved_vel, saved_iface_z, n_iface_parms)
        rms_new = rms_orig
        write(*,'(A)') '    All steps failed, reverting to original model'
      end if

      ! Scale dmodel_total: velocity uses step_len, interface uses step_iface
      step_iface = step_len * iface_step_ratio_in
      if (n_inv_iface == 0) step_iface = 0.0_dp
      do j = 1, nparm
        if (j <= n_vel_parms) then
          dmodel_total(j) = step_len * dmodel_total(j)
        else
          if (n_inv_iface > 0) then
            dmodel_total(j) = step_iface * dmodel_total(j)
          else
            dmodel_total(j) = step_len * dmodel_total(j)
          end if
        end if
      end do
      dmodel_total_sum(1:nparm) = dmodel_total_sum(1:nparm) + dmodel_total(1:nparm)

      deallocate(gtg_diag, smooth_work, smooth_out)
      deallocate(u, v_vec, w_vec, z_vec, Av_aug, Atu_aug, Av_data)
      deallocate(dm_trial)
    end block

    deallocate(gindx, gval, residual, row_ptr)
    deallocate(data_wt)
    deallocate(node_hit)

    write(*,'(A,I4,A,I6,A,F14.8,A,F14.8,A,ES10.3,A,ES10.3,A,ES10.3)') 'Iteration ', iter, &
      ' ndata=', ndata, ' RMS=', rms_orig, ' RMS_new=', rms_new, &
      ' avg_dmv=', avg_dmv, ' wdv=', wd_best, ' wdd=', wd_best_iface
    write(99,'(A,I4,A,I6,A,F14.8,A,F14.8,A,ES10.3,A,ES10.3,A,ES10.3)') 'Iteration ', iter, &
      ' ndata=', ndata, ' RMS=', rms_orig, ' RMS_new=', rms_new, &
      ' avg_dmv=', avg_dmv, ' wdv=', wd_best, ' wdd=', wd_best_iface
    flush(99)

    ! Track best model (lowest RMS_new)
    if (rms_new < best_rms) then
      best_rms = rms_new
      best_iter = iter
      ! Save current velocity and interface
      if (.not. allocated(best_vel)) then
        ! Use first vgrid dimensions for allocation
        vg => vgrid(vgrids_to_be_inv(1)%igrid, vgrids_to_be_inv(1)%vtype)
        allocate(best_vel(vg%nx, vg%nz))
        if (n_inv_iface > 0) allocate(best_iface_z(n_iface_parms))
      end if
      do ig = 1, n_inv_vgrid
        vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
        vtype = vgrids_to_be_inv(ig)%vtype
        nx_vg = vg%nx; nz_vg = vg%nz
        do iz = 1, nz_vg
          do ix = 1, nx_vg
            best_vel(ix, iz) = vg%velocity(vtype, ix, iz)
          end do
        end do
      end do
      if (n_inv_iface > 0) then
        ip = 0
        do ii = 1, n_inv_iface
          iiface = ifaces_to_be_inv(ii)
          ifc => intrface(iiface)
          do ix = 1, ifc%nx
            ip = ip + 1
            best_iface_z(ip) = ifc%z_nodes(ix)
          end do
        end do
      end if
      write(*,'(A,I4,A,F14.8)') '    ** New best model at iter ', iter, ' RMS=', rms_new
    end if

    if (rms_orig < 1.0e-4_dp) exit

    if (iter > 1 .and. prev_rms > 0.0_dp) then
      if (prev_rms > rms_orig .and. (prev_rms - rms_orig) / prev_rms < 1.0e-5_dp .and. iter > 5) exit
    end if
    prev_rms = rms_orig

  end do

  ! Restore best model if it's better than the final one
  if (allocated(best_vel) .and. best_iter > 0 .and. best_iter < iter) then
    block
      integer :: ip2, ii2, iiface2, ix2, iz2, ig2
      write(*,'(A,I4,A,F14.8,A,I4)') '  Restoring best model from iter ', best_iter, &
        ' RMS=', best_rms, ' (final iter=', iter, ')'
      do ig2 = 1, n_inv_vgrid
        vg => vgrid(vgrids_to_be_inv(ig2)%igrid, vgrids_to_be_inv(ig2)%vtype)
        vtype = vgrids_to_be_inv(ig2)%vtype
        nx_vg = vg%nx; nz_vg = vg%nz
        do iz2 = 1, nz_vg
          do ix2 = 1, nx_vg
            vg%velocity(vtype, ix2, iz2) = best_vel(ix2, iz2)
          end do
        end do
      end do
      if (n_inv_iface > 0 .and. allocated(best_iface_z)) then
        ip2 = 0
        do ii2 = 1, n_inv_iface
          iiface2 = ifaces_to_be_inv(ii2)
          ifc => intrface(iiface2)
          do ix2 = 1, ifc%nx
            ip2 = ip2 + 1
            ifc%z_nodes(ix2) = best_iface_z(ip2)
          end do
        end do
        call rebuild_geometry2d()
      end if
    end block
    deallocate(best_vel)
    if (allocated(best_iface_z)) deallocate(best_iface_z)
  end if

  deallocate(mvscale)
  deallocate(dmodel_total)
  deallocate(dmodel_total_sum)
  deallocate(cm, ecmi)

  close(99)

  call write_model2d(iter)

  ! Write final residuals and ray density (after inversion)
  call run_forward2d()
  call write_residuals2d('residuals_final.dat')
  call write_ray_density2d('ray_density_final.dat')

end subroutine lsqr_inversion2d


subroutine subspace_solver_2d(ndata, nparm, nnz, gindx, gval, &
  row_ptr, b, epsv, etav, epsilon, mvscale, dmodel_total, delta_m, subdim)

  implicit none
  integer, intent(in) :: ndata, nparm, nnz, subdim
  integer, intent(in) :: gindx(nnz), row_ptr(ndata+1)
  real(kind=dp), intent(in) :: gval(nnz), b(ndata)
  real(kind=dp), intent(in) :: epsv, etav, epsilon
  real(kind=dp), intent(in) :: mvscale(nparm), dmodel_total(nparm)
  real(kind=dp), intent(out) :: delta_m(nparm)

  integer :: i, j, k, actual_subdim
  real(kind=dp), allocatable :: gamma(:), a(:,:), ga(:), gta(:)
  real(kind=dp), allocatable :: mati(:,:), r1mat(:), r2mat(:)
  real(kind=dp), allocatable :: smv(:), dem(:)
  real(kind=dp) :: rsum, rnorm, ecmi_j, rnorm_data, rnorm_damp, rnorm_smooth

  allocate(gamma(nparm), a(nparm, subdim))
  a = 0.0_dp
  gamma = 0.0_dp

  call sparse_ATx(ndata, nparm, nnz, gindx, gval, row_ptr, b, gamma)

  rnorm_data = sqrt(dot_product(gamma, gamma))

  do j = 1, nparm
    gamma(j) = gamma(j) - epsv * dmodel_total(j)
  end do

  rnorm_damp = 0.0_dp
  do j = 1, nparm
    rnorm_damp = rnorm_damp + (epsv * dmodel_total(j))**2
  end do
  rnorm_damp = sqrt(rnorm_damp)

  if (etav > 0.0_dp) then
    allocate(smv(nparm), dem(nparm))
    smv = dmodel_total
    call apply_laplacian_smooth_2d(nparm, smv, dem)
    smv = 0.0_dp
    call apply_laplacian_smooth_2d(nparm, dem, smv)
    rnorm_smooth = sqrt(dot_product(smv, smv)) * etav
    do j = 1, nparm
      gamma(j) = gamma(j) - etav * smv(j)
    end do
    deallocate(smv, dem)
  else
    rnorm_smooth = 0.0_dp
  end if

  rnorm = sqrt(dot_product(gamma, gamma))
  write(*,'(A,ES12.4)') '    |gradient|=', rnorm
  write(*,'(A,ES12.4,A,ES12.4,A,ES12.4)') '    |data|=', rnorm_data, &
    ' |damp|=', rnorm_damp, ' |smooth|=', rnorm_smooth

  if (rnorm < 1.0e-30_dp) then
    delta_m = 0.0_dp
    deallocate(gamma, a)
    return
  end if

  do j = 1, nparm
    a(j, 1) = gamma(j)
  end do

  rsum = sqrt(dot_product(a(:,1), a(:,1)))
  if (rsum > 0.0_dp) a(:,1) = a(:,1) / rsum

  allocate(ga(ndata), gta(nparm))
  if (etav > 0.0_dp) then
    allocate(smv(nparm), dem(nparm))
  end if

  actual_subdim = subdim
  if (subdim > 1) then
    do i = 2, subdim
      call sparse_Ax(ndata, nparm, nnz, gindx, gval, row_ptr, a(:,i-1), ga)
      call sparse_ATx(ndata, nparm, nnz, gindx, gval, row_ptr, ga, gta)

      do j = 1, nparm
        a(j, i) = gta(j) + epsilon * a(j, i-1)
      end do

      if (etav > 0.0_dp) then
        smv = a(:, i-1)
        call apply_laplacian_smooth_2d(nparm, smv, dem)
        smv = 0.0_dp
        call apply_laplacian_smooth_2d(nparm, dem, smv)
        do j = 1, nparm
          a(j, i) = a(j, i) + etav * smv(j)
        end do
      end if

      rsum = sqrt(dot_product(a(:,i), a(:,i)))
      if (rsum > 0.0_dp) then
        a(:,i) = a(:,i) / rsum
      else
        actual_subdim = i - 1
        exit
      end if
    end do
  end if

  call svd_orthogonalize(a, nparm, actual_subdim)

  allocate(mati(actual_subdim, actual_subdim))
  allocate(r1mat(actual_subdim), r2mat(actual_subdim))
  mati = 0.0_dp

  do i = 1, actual_subdim
    call sparse_Ax(ndata, nparm, nnz, gindx, gval, row_ptr, a(:,i), ga)
    call sparse_ATx(ndata, nparm, nnz, gindx, gval, row_ptr, ga, gta)

    if (etav > 0.0_dp) then
      smv = a(:, i)
      call apply_laplacian_smooth_2d(nparm, smv, dem)
      smv = 0.0_dp
      call apply_laplacian_smooth_2d(nparm, dem, smv)
    end if

    do j = 1, nparm
      rsum = gta(j) + epsv * a(j, i)
      if (etav > 0.0_dp) then
        rsum = rsum + etav * smv(j)
      end if
      do k = 1, actual_subdim
        mati(k, i) = mati(k, i) + rsum * a(j, k)
      end do
    end do
  end do

  if (etav > 0.0_dp) then
    deallocate(smv, dem)
  end if

  r1mat = 0.0_dp
  do i = 1, actual_subdim
    do j = 1, nparm
      r1mat(i) = r1mat(i) + a(j, i) * gamma(j)
    end do
  end do

  call solve_small_matrix(mati, r1mat, r2mat, actual_subdim)

  delta_m = 0.0_dp
  do i = 1, actual_subdim
    do j = 1, nparm
      delta_m(j) = delta_m(j) + a(j, i) * r2mat(i)
    end do
  end do

  write(*,'(A,I3,A,ES12.4)') '    subspace dim=', actual_subdim, &
    ' |delta_m|=', sqrt(dot_product(delta_m, delta_m))

  deallocate(gamma, a, ga, gta, mati, r1mat, r2mat)

end subroutine subspace_solver_2d


subroutine apply_laplacian_smooth_2d(nparm, vin, vout)
  implicit none
  integer, intent(in) :: nparm
  real(kind=dp), intent(in) :: vin(nparm)
  real(kind=dp), intent(out) :: vout(nparm)

  integer :: ig, ix, iz, vtype, nx_vg, nz_vg, idx, jdx
  integer :: n_neighbors
  type(Tvelocity_grid2d), pointer :: vg
  integer :: ii, iiface, nx_if
  type(Tinterface2d), pointer :: ifc

  vout = 0.0_dp

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx <= 0) cycle
        n_neighbors = 0
        if (ix > 1) then
          jdx = vg%inv_idx(ix-1, iz)
          if (jdx > 0) then
            vout(idx) = vout(idx) + vin(jdx)
            n_neighbors = n_neighbors + 1
          end if
        end if
        if (ix < nx_vg) then
          jdx = vg%inv_idx(ix+1, iz)
          if (jdx > 0) then
            vout(idx) = vout(idx) + vin(jdx)
            n_neighbors = n_neighbors + 1
          end if
        end if
        if (iz > 1) then
          jdx = vg%inv_idx(ix, iz-1)
          if (jdx > 0) then
            vout(idx) = vout(idx) + vin(jdx)
            n_neighbors = n_neighbors + 1
          end if
        end if
        if (iz < nz_vg) then
          jdx = vg%inv_idx(ix, iz+1)
          if (jdx > 0) then
            vout(idx) = vout(idx) + vin(jdx)
            n_neighbors = n_neighbors + 1
          end if
        end if
        vout(idx) = vout(idx) - dble(n_neighbors) * vin(idx)
      end do
    end do
  end do

  ! Interface smoothing: 1D Laplacian along x
  do ii = 1, n_inv_iface
    iiface = ifaces_to_be_inv(ii)
    ifc => intrface(iiface)
    nx_if = ifc%nx
    do ix = 1, nx_if
      idx = ifc%start_index + ix
      if (idx < 1 .or. idx > nparm) cycle
      if (ix == 1) then
        ! One-sided: forward difference
        jdx = ifc%start_index + ix + 1
        if (jdx <= nparm) then
          vout(idx) = vout(idx) + vin(jdx) - vin(idx)
        end if
      else if (ix == nx_if) then
        ! One-sided: backward difference
        jdx = ifc%start_index + ix - 1
        if (jdx >= 1) then
          vout(idx) = vout(idx) + vin(jdx) - vin(idx)
        end if
      else
        ! Central: second derivative
        jdx = ifc%start_index + ix - 1
        vout(idx) = vout(idx) + vin(jdx)
        jdx = ifc%start_index + ix + 1
        vout(idx) = vout(idx) + vin(jdx)
        vout(idx) = vout(idx) - 2.0_dp * vin(idx)
      end if
    end do
  end do

end subroutine apply_laplacian_smooth_2d


subroutine build_gauss_smooth_matrix(nparm, corr_x, corr_z, mvscale_arr)
  ! Pre-compute Gaussian smoothing sparse matrices in CSR format.
  ! R_h: horizontal smoothing, R_v: vertical smoothing.
  ! R = R_v * R_h (applied as separable).
  ! mvscale is baked into the matrix values.
  implicit none
  integer, intent(in) :: nparm
  real(kind=dp), intent(in) :: corr_x, corr_z
  real(kind=dp), intent(in) :: mvscale_arr(nparm)

  integer :: ig, ix, iz, jx, jz, idx, jdx, k, nnz_est
  integer :: nx_vg, nz_vg, kx_min, kx_max
  real(kind=dp) :: dx_vg, dist_x, dist_z, beta, beta_sum
  real(kind=dp), allocatable :: diag_h(:), diag_v(:)
  type(Tvelocity_grid2d), pointer :: vg

  if (gs_built .and. gs_nparm == nparm) return  ! already built

  ! Clean up previous allocation if any
  if (allocated(gs_h_rowptr)) deallocate(gs_h_rowptr)
  if (allocated(gs_h_colind)) deallocate(gs_h_colind)
  if (allocated(gs_h_values)) deallocate(gs_h_values)
  if (allocated(gs_v_rowptr)) deallocate(gs_v_rowptr)
  if (allocated(gs_v_colind)) deallocate(gs_v_colind)
  if (allocated(gs_v_values)) deallocate(gs_v_values)

  gs_nparm = nparm
  allocate(gs_h_rowptr(nparm+1))
  allocate(gs_v_rowptr(nparm+1))
  allocate(diag_h(nparm), diag_v(nparm))

  ! Initialize entire rowptr to 0 before counting (interface rows stay 0)
  gs_h_rowptr = 0
  gs_v_rowptr = 0

  ! ---- Pass 1: count non-zeros for R_h and R_v ----
  gs_h_rowptr(1) = 1   ! 1-based CSR
  gs_v_rowptr(1) = 1

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    nx_vg = vg%nx
    nz_vg = vg%nz
    dx_vg = vg%dx
    kx_min = -int(3.0_dp * corr_x / dx_vg) - 1
    kx_max = int(3.0_dp * corr_x / dx_vg) + 1

    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx <= 0) cycle

        ! Count R_h non-zeros for this row
        nnz_est = 1  ! diagonal
        do jx = max(1, ix + kx_min), min(nx_vg, ix + kx_max)
          if (jx == ix) cycle
          jdx = vg%inv_idx(jx, iz)
          if (jdx > 0) nnz_est = nnz_est + 1
        end do
        gs_h_rowptr(idx+1) = nnz_est

        ! Count R_v non-zeros for this row
        nnz_est = 1  ! diagonal
        do jz = 1, nz_vg
          if (jz == iz) cycle
          dist_z = vg%z(jz) - vg%z(iz)
          if (abs(dist_z) > 3.0_dp * corr_z) cycle
          jdx = vg%inv_idx(ix, jz)
          if (jdx > 0) nnz_est = nnz_est + 1
        end do
        gs_v_rowptr(idx+1) = nnz_est
      end do
    end do
  end do

  ! Convert counts to CSR row pointers (cumulative sum)
  do idx = 1, nparm
    gs_h_rowptr(idx+1) = gs_h_rowptr(idx) + gs_h_rowptr(idx+1)
    gs_v_rowptr(idx+1) = gs_v_rowptr(idx) + gs_v_rowptr(idx+1)
  end do
  gs_h_nnz = gs_h_rowptr(nparm+1) - 1
  gs_v_nnz = gs_v_rowptr(nparm+1) - 1

  allocate(gs_h_colind(gs_h_nnz), gs_h_values(gs_h_nnz))
  allocate(gs_v_colind(gs_v_nnz), gs_v_values(gs_v_nnz))

  ! ---- Pass 2: fill values ----
  diag_h = 0.0_dp
  diag_v = 0.0_dp

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    nx_vg = vg%nx
    nz_vg = vg%nz
    dx_vg = vg%dx
    kx_min = -int(3.0_dp * corr_x / dx_vg) - 1
    kx_max = int(3.0_dp * corr_x / dx_vg) + 1

    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx <= 0) cycle

        ! Fill R_h row
        k = gs_h_rowptr(idx)
        beta_sum = 0.0_dp
        do jx = max(1, ix + kx_min), min(nx_vg, ix + kx_max)
          if (jx == ix) cycle
          jdx = vg%inv_idx(jx, iz)
          if (jdx <= 0) cycle
          dist_x = vg%x(jx) - vg%x(ix)
          beta = exp(-dist_x**2 / (corr_x**2))
          beta_sum = beta_sum + beta
          gs_h_colind(k) = jdx
          gs_h_values(k) = beta * mvscale_arr(jdx)
          k = k + 1
        end do
        ! Diagonal
        gs_h_colind(k) = idx
        gs_h_values(k) = -beta_sum * mvscale_arr(idx)
        diag_h(idx) = gs_h_values(k)

        ! Fill R_v row
        k = gs_v_rowptr(idx)
        beta_sum = 0.0_dp
        do jz = 1, nz_vg
          if (jz == iz) cycle
          dist_z = vg%z(jz) - vg%z(iz)
          if (abs(dist_z) > 3.0_dp * corr_z) cycle
          jdx = vg%inv_idx(ix, jz)
          if (jdx <= 0) cycle
          beta = exp(-dist_z**2 / (corr_z**2))
          beta_sum = beta_sum + beta
          gs_v_colind(k) = jdx
          gs_v_values(k) = beta * mvscale_arr(jdx)
          k = k + 1
        end do
        ! Diagonal
        gs_v_colind(k) = idx
        gs_v_values(k) = -beta_sum * mvscale_arr(idx)
        diag_v(idx) = gs_v_values(k)
      end do
    end do
  end do

  deallocate(diag_h, diag_v)
  gs_built = .true.

  write(*,'(A,I8,A,I8,A)') '  Gaussian smooth matrix built: R_h nnz=', gs_h_nnz, ', R_v nnz=', gs_v_nnz

end subroutine build_gauss_smooth_matrix


subroutine apply_gauss_smooth_2d(nparm, vin, vout, corr_x, corr_z, mvscale_arr)
  ! Gaussian smoothing using pre-computed sparse matrices.
  ! R = R_v * R_h (separable), with mvscale baked in.
  implicit none
  integer, intent(in) :: nparm
  real(kind=dp), intent(in) :: vin(nparm)
  real(kind=dp), intent(out) :: vout(nparm)
  real(kind=dp), intent(in) :: corr_x, corr_z
  real(kind=dp), intent(in) :: mvscale_arr(nparm)

  integer :: i, k
  real(kind=dp), allocatable :: temp(:)

  if (.not. gs_built .or. gs_nparm /= nparm) then
    ! Fallback: build matrix on first call
    call build_gauss_smooth_matrix(nparm, corr_x, corr_z, mvscale_arr)
  end if

  vout = 0.0_dp
  allocate(temp(nparm))
  temp = 0.0_dp

  ! Step 1: temp = R_h * vin (sparse matrix-vector multiply)
  do i = 1, nparm
    do k = gs_h_rowptr(i), gs_h_rowptr(i+1) - 1
      temp(i) = temp(i) + gs_h_values(k) * vin(gs_h_colind(k))
    end do
  end do

  ! Step 2: vout = R_v * temp (sparse matrix-vector multiply)
  do i = 1, nparm
    do k = gs_v_rowptr(i), gs_v_rowptr(i+1) - 1
      vout(i) = vout(i) + gs_v_values(k) * temp(gs_v_colind(k))
    end do
  end do

  deallocate(temp)

end subroutine apply_gauss_smooth_2d


subroutine apply_gauss_smooth_h(nparm, vin, vout)
  ! Apply horizontal smoothing R_h only (sparse SpMV)
  implicit none
  integer, intent(in) :: nparm
  real(kind=dp), intent(in) :: vin(nparm)
  real(kind=dp), intent(out) :: vout(nparm)

  integer :: i, k

  vout = 0.0_dp
  do i = 1, nparm
    do k = gs_h_rowptr(i), gs_h_rowptr(i+1) - 1
      vout(i) = vout(i) + gs_h_values(k) * vin(gs_h_colind(k))
    end do
  end do

end subroutine apply_gauss_smooth_h


subroutine apply_gauss_smooth_v(nparm, vin, vout)
  ! Apply vertical smoothing R_v only (sparse SpMV)
  implicit none
  integer, intent(in) :: nparm
  real(kind=dp), intent(in) :: vin(nparm)
  real(kind=dp), intent(out) :: vout(nparm)

  integer :: i, k

  vout = 0.0_dp
  do i = 1, nparm
    do k = gs_v_rowptr(i), gs_v_rowptr(i+1) - 1
      vout(i) = vout(i) + gs_v_values(k) * vin(gs_v_colind(k))
    end do
  end do

end subroutine apply_gauss_smooth_v


subroutine apply_gauss_smooth_h_T(nparm, vin, vout)
  ! Apply transposed horizontal smoothing: vout = M_h^T * vin.
  ! Needed for LSQR adjoint: mvscale is baked as column scaling (M = R*D),
  ! so M^T = D*R^T, which differs from M when mvscale is non-uniform.
  implicit none
  integer, intent(in) :: nparm
  real(kind=dp), intent(in) :: vin(nparm)
  real(kind=dp), intent(out) :: vout(nparm)

  integer :: i, k

  vout = 0.0_dp
  do i = 1, nparm
    do k = gs_h_rowptr(i), gs_h_rowptr(i+1) - 1
      vout(gs_h_colind(k)) = vout(gs_h_colind(k)) + gs_h_values(k) * vin(i)
    end do
  end do

end subroutine apply_gauss_smooth_h_T


subroutine apply_gauss_smooth_v_T(nparm, vin, vout)
  ! Apply transposed vertical smoothing: vout = M_v^T * vin (see above).
  implicit none
  integer, intent(in) :: nparm
  real(kind=dp), intent(in) :: vin(nparm)
  real(kind=dp), intent(out) :: vout(nparm)

  integer :: i, k

  vout = 0.0_dp
  do i = 1, nparm
    do k = gs_v_rowptr(i), gs_v_rowptr(i+1) - 1
      vout(gs_v_colind(k)) = vout(gs_v_colind(k)) + gs_v_values(k) * vin(i)
    end do
  end do

end subroutine apply_gauss_smooth_v_T


subroutine build_iface_gauss_smooth_matrix(nparm, n_vel_parms, corr_x_iface, mvscale_arr)
  ! Pre-compute interface Gaussian smoothing sparse matrix (1D horizontal, R_h only).
  ! Only fills rows for interface parameters (indices n_vel_parms+1 : nparm).
  ! mvscale is baked into the matrix values.
  implicit none
  integer, intent(in) :: nparm, n_vel_parms
  real(kind=dp), intent(in) :: corr_x_iface
  real(kind=dp), intent(in) :: mvscale_arr(nparm)

  integer :: ii, iiface, ix, jx, idx, jdx, k, nnz_est, nx_if, kx_min, kx_max
  real(kind=dp) :: dx_if, dist_x, beta, beta_sum
  type(Tinterface2d), pointer :: ifc

  if (gsi_built .and. gsi_nparm == nparm) return
  if (corr_x_iface <= 0.0_dp) return

  if (allocated(gsi_h_rowptr)) deallocate(gsi_h_rowptr)
  if (allocated(gsi_h_colind)) deallocate(gsi_h_colind)
  if (allocated(gsi_h_values)) deallocate(gsi_h_values)

  gsi_nparm = nparm
  allocate(gsi_h_rowptr(nparm+1))
  gsi_h_rowptr = 0
  gsi_h_rowptr(1) = 1   ! 1-based CSR (critical: must be 1, not 0)

  ! Pass 1: count non-zeros for interface rows
  do ii = 1, n_inv_iface
    iiface = ifaces_to_be_inv(ii)
    ifc => intrface(iiface)
    nx_if = ifc%nx
    dx_if = ifc%dx
    kx_min = -int(3.0_dp * corr_x_iface / dx_if) - 1
    kx_max = int(3.0_dp * corr_x_iface / dx_if) + 1
    do ix = 1, nx_if
      idx = ifc%start_index + ix
      if (idx < 1 .or. idx > nparm) cycle
      nnz_est = 1  ! diagonal
      do jx = max(1, ix + kx_min), min(nx_if, ix + kx_max)
        if (jx == ix) cycle
        jdx = ifc%start_index + jx
        if (jdx >= 1 .and. jdx <= nparm) nnz_est = nnz_est + 1
      end do
      gsi_h_rowptr(idx+1) = nnz_est
    end do
  end do

  ! Convert counts to CSR row pointers (cumulative sum)
  do idx = 1, nparm
    gsi_h_rowptr(idx+1) = gsi_h_rowptr(idx) + gsi_h_rowptr(idx+1)
  end do
  gsi_h_nnz = gsi_h_rowptr(nparm+1) - 1

  if (gsi_h_nnz <= 0) then
    gsi_built = .true.
    return
  end if

  allocate(gsi_h_colind(gsi_h_nnz), gsi_h_values(gsi_h_nnz))

  ! Pass 2: fill values
  do ii = 1, n_inv_iface
    iiface = ifaces_to_be_inv(ii)
    ifc => intrface(iiface)
    nx_if = ifc%nx
    dx_if = ifc%dx
    kx_min = -int(3.0_dp * corr_x_iface / dx_if) - 1
    kx_max = int(3.0_dp * corr_x_iface / dx_if) + 1
    do ix = 1, nx_if
      idx = ifc%start_index + ix
      if (idx < 1 .or. idx > nparm) cycle
      k = gsi_h_rowptr(idx)
      beta_sum = 0.0_dp
      do jx = max(1, ix + kx_min), min(nx_if, ix + kx_max)
        if (jx == ix) cycle
        jdx = ifc%start_index + jx
        if (jdx < 1 .or. jdx > nparm) cycle
        dist_x = ifc%x(jx) - ifc%x(ix)
        beta = exp(-dist_x**2 / (corr_x_iface**2))
        beta_sum = beta_sum + beta
        gsi_h_colind(k) = jdx
        gsi_h_values(k) = beta * mvscale_arr(jdx)
        k = k + 1
      end do
      ! Diagonal
      gsi_h_colind(k) = idx
      gsi_h_values(k) = -beta_sum * mvscale_arr(idx)
    end do
  end do

  gsi_built = .true.
  write(*,'(A,I8,A,F8.3)') '  Interface Gaussian smooth matrix built: R_h nnz=', gsi_h_nnz, ' corr_x_iface=', corr_x_iface

end subroutine build_iface_gauss_smooth_matrix


subroutine apply_iface_gauss_smooth_h(nparm, vin, vout)
  ! Apply interface horizontal Gaussian smoothing R_h (sparse SpMV)
  implicit none
  integer, intent(in) :: nparm
  real(kind=dp), intent(in) :: vin(nparm)
  real(kind=dp), intent(out) :: vout(nparm)

  integer :: i, k

  vout = 0.0_dp
  if (.not. gsi_built .or. gsi_nparm /= nparm) return
  do i = 1, nparm
    do k = gsi_h_rowptr(i), gsi_h_rowptr(i+1) - 1
      vout(i) = vout(i) + gsi_h_values(k) * vin(gsi_h_colind(k))
    end do
  end do

end subroutine apply_iface_gauss_smooth_h


subroutine svd_orthogonalize(a, n, subdim)
  implicit none
  integer, intent(in) :: n, subdim
  real(kind=dp), intent(inout) :: a(n, subdim)

  integer :: i, j, k, isum, jj
  integer, allocatable :: wnz(:)
  real(kind=dp), allocatable :: v(:,:), w(:), rv1(:)
  real(kind=dp) :: anorm, c, f, g, h, s, scale, x, y, z
  real(kind=dp), parameter :: wtol = 1.0e-10_dp
  integer :: its, l, m, nm, isw
  integer :: m_min

  m_min = min(n, subdim)
  allocate(v(subdim, subdim), w(subdim), rv1(subdim), wnz(subdim))

  g = 0.0_dp
  scale = 0.0_dp
  anorm = 0.0_dp

  do i = 1, subdim
    l = i + 1
    rv1(i) = scale * g
    g = 0.0_dp
    s = 0.0_dp
    scale = 0.0_dp
    if (i <= n) then
      do k = i, n
        scale = scale + abs(a(k, i))
      end do
      if (scale /= 0.0_dp) then
        do k = i, n
          a(k, i) = a(k, i) / scale
          s = s + a(k, i) * a(k, i)
        end do
        f = a(i, i)
        g = -sign(sqrt(s), f)
        h = f * g - s
        a(i, i) = f - g
        do j = l, subdim
          s = 0.0_dp
          do k = i, n
            s = s + a(k, i) * a(k, j)
          end do
          f = s / h
          do k = i, n
            a(k, j) = a(k, j) + f * a(k, i)
          end do
        end do
        do k = i, n
          a(k, i) = scale * a(k, i)
        end do
      end if
    end if
    w(i) = scale * g
    g = 0.0_dp
    s = 0.0_dp
    scale = 0.0_dp
    if (i <= n .and. i /= subdim) then
      do k = l, subdim
        scale = scale + abs(a(i, k))
      end do
      if (scale /= 0.0_dp) then
        do k = l, subdim
          a(i, k) = a(i, k) / scale
          s = s + a(i, k) * a(i, k)
        end do
        f = a(i, l)
        g = -sign(sqrt(s), f)
        h = f * g - s
        a(i, l) = f - g
        do k = l, subdim
          rv1(k) = a(i, k) / h
        end do
        do j = l, n
          s = 0.0_dp
          do k = l, subdim
            s = s + a(j, k) * a(i, k)
          end do
          do k = l, subdim
            a(j, k) = a(j, k) + s * rv1(k)
          end do
        end do
        do k = l, subdim
          a(i, k) = scale * a(i, k)
        end do
      end if
    end if
    anorm = max(anorm, abs(w(i)) + abs(rv1(i)))
  end do

  do i = subdim, 1, -1
    if (i < subdim) then
      if (g /= 0.0_dp) then
        do j = l, subdim
          v(j, i) = (a(i, j) / a(i, l)) / g
        end do
        do j = l, subdim
          s = 0.0_dp
          do k = l, subdim
            s = s + a(i, k) * v(k, j)
          end do
          do k = l, subdim
            v(k, j) = v(k, j) + s * v(k, i)
          end do
        end do
      end if
      do j = l, subdim
        v(i, j) = 0.0_dp
        v(j, i) = 0.0_dp
      end do
    end if
    v(i, i) = 1.0_dp
    g = rv1(i)
    l = i
  end do

  do i = m_min, 1, -1
    l = i + 1
    g = w(i)
    do j = l, subdim
      a(i, j) = 0.0_dp
    end do
    if (g /= 0.0_dp) then
      g = 1.0_dp / g
      do j = l, subdim
        s = 0.0_dp
        do k = l, n
          s = s + a(k, i) * a(k, j)
        end do
        f = (s / a(i, i)) * g
        do k = i, n
          a(k, j) = a(k, j) + f * a(k, i)
        end do
      end do
      do j = i, n
        a(j, i) = a(j, i) * g
      end do
    else
      do j = i, n
        a(j, i) = 0.0_dp
      end do
    end if
    a(i, i) = a(i, i) + 1.0_dp
  end do

  do k = subdim, 1, -1
    do its = 1, 30
      isw = 0
      do l = k, 1, -1
        nm = l - 1
        if (abs(rv1(l)) + anorm == anorm) then
          isw = 1
          exit
        end if
        if (abs(w(nm)) + anorm == anorm) exit
      end do
      if (isw == 0) then
        c = 0.0_dp
        s = 1.0_dp
        do i = l, k
          f = s * rv1(i)
          rv1(i) = c * rv1(i)
          if (abs(f) + anorm == anorm) exit
          g = w(i)
          h = sqrt(f*f + g*g)
          w(i) = h
          h = 1.0_dp / h
          c = g * h
          s = -f * h
          do j = 1, n
            y = a(j, nm)
            z = a(j, i)
            a(j, nm) = y * c + z * s
            a(j, i) = -y * s + z * c
          end do
        end do
      end if
      z = w(k)
      if (l == k) then
        if (z < 0.0_dp) then
          w(k) = -z
          do j = 1, subdim
            v(j, k) = -v(j, k)
          end do
        end if
        exit
      end if
      x = w(l)
      nm = k - 1
      y = w(nm)
      g = rv1(nm)
      h = rv1(k)
      f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0_dp * h * y)
      g = sqrt(f*f + 1.0_dp)
      f = ((x - z) * (x + z) + h * ((y / (f + sign(g, f))) - h)) / x
      c = 1.0_dp
      s = 1.0_dp
      do j = l, nm
        i = j + 1
        g = rv1(i)
        y = w(i)
        h = s * g
        g = c * g
        z = sqrt(f*f + h*h)
        rv1(j) = z
        c = f / z
        s = h / z
        f = x * c + g * s
        g = -x * s + g * c
        h = y * s
        y = y * c
        do jj = 1, subdim
          x = v(jj, j)
          z = v(jj, i)
          v(jj, j) = x * c + z * s
          v(jj, i) = -x * s + z * c
        end do
        z = sqrt(f*f + h*h)
        w(j) = z
        if (z /= 0.0_dp) then
          z = 1.0_dp / z
          c = f * z
          s = h * z
        end if
        f = c * g + s * y
        x = -s * g + c * y
        do jj = 1, n
          y = a(jj, j)
          z = a(jj, i)
          a(jj, j) = y * c + z * s
          a(jj, i) = -y * s + z * c
        end do
      end do
      rv1(l) = 0.0_dp
      rv1(k) = f
      w(k) = x
    end do
  end do

  isum = 0
  do i = 1, subdim
    if (abs(w(i)) >= wtol) then
      isum = isum + 1
      wnz(isum) = i
    end if
  end do
  if (isum < subdim) then
    write(*,'(A,I3,A,I3)') '    SVD: subspace reduced from ', subdim, ' to ', isum
  end if

  deallocate(v, w, rv1, wnz)

end subroutine svd_orthogonalize


subroutine solve_small_matrix(a, b, x, n)
  implicit none
  integer, intent(in) :: n
  real(kind=dp), intent(inout) :: a(n, n)
  real(kind=dp), intent(in) :: b(n)
  real(kind=dp), intent(out) :: x(n)

  integer :: i, j, k, imax
  integer, allocatable :: indx(:)
  real(kind=dp) :: aamax, dum, sum_val
  real(kind=dp), allocatable :: vv(:), y(:,:)
  real(kind=dp), parameter :: tiny = 1.0e-20_dp

  allocate(indx(n), vv(n), y(n, n))

  y = 0.0_dp
  do i = 1, n
    y(i, i) = 1.0_dp
  end do

  do i = 1, n
    aamax = 0.0_dp
    do j = 1, n
      if (abs(a(i, j)) > aamax) aamax = abs(a(i, j))
    end do
    if (aamax == 0.0_dp) then
      write(*,*) '    Singular matrix in solve_small_matrix'
      x = 0.0_dp
      deallocate(indx, vv, y)
      return
    end if
    vv(i) = 1.0_dp / aamax
  end do

  do j = 1, n
    do i = 1, j - 1
      sum_val = a(i, j)
      do k = 1, i - 1
        sum_val = sum_val - a(i, k) * a(k, j)
      end do
      a(i, j) = sum_val
    end do
    aamax = 0.0_dp
    imax = j
    do i = j, n
      sum_val = a(i, j)
      do k = 1, j - 1
        sum_val = sum_val - a(i, k) * a(k, j)
      end do
      a(i, j) = sum_val
      dum = vv(i) * abs(sum_val)
      if (dum >= aamax) then
        imax = i
        aamax = dum
      end if
    end do
    if (j /= imax) then
      do k = 1, n
        dum = a(imax, k)
        a(imax, k) = a(j, k)
        a(j, k) = dum
      end do
      vv(imax) = vv(j)
    end if
    indx(j) = imax
    if (a(j, j) == 0.0_dp) a(j, j) = tiny
    if (j /= n) then
      dum = 1.0_dp / a(j, j)
      do i = j + 1, n
        a(i, j) = a(i, j) * dum
      end do
    end if
  end do

  do j = 1, n
    call lubksb_local(a, n, indx, y(1, j))
  end do

  do i = 1, n
    do j = 1, n
      a(i, j) = y(i, j)
    end do
  end do

  x = 0.0_dp
  do i = 1, n
    do j = 1, n
      x(i) = x(i) + a(i, j) * b(j)
    end do
  end do

  deallocate(indx, vv, y)

end subroutine solve_small_matrix


subroutine lubksb_local(a, n, indx, b)
  implicit none
  integer, intent(in) :: n, indx(n)
  real(kind=dp), intent(in) :: a(n, n)
  real(kind=dp), intent(inout) :: b(n)

  integer :: i, j, ii, ll
  real(kind=dp) :: sum_val

  ii = 0
  do i = 1, n
    ll = indx(i)
    sum_val = b(ll)
    b(ll) = b(i)
    if (ii /= 0) then
      do j = ii, i - 1
        sum_val = sum_val - a(i, j) * b(j)
      end do
    else if (sum_val /= 0.0_dp) then
      ii = 1
    end if
    b(i) = sum_val
  end do
  do i = n, 1, -1
    sum_val = b(i)
    do j = i + 1, n
      sum_val = sum_val - a(i, j) * b(j)
    end do
    b(i) = sum_val / a(i, i)
  end do

end subroutine lubksb_local


subroutine build_laplacian_smoothing(nparm, n_L_nnz, L_row, L_col, L_val)
  implicit none
  integer, intent(in) :: nparm
  integer, intent(out) :: n_L_nnz
  integer, allocatable, intent(out) :: L_row(:), L_col(:)
  real(kind=dp), allocatable, intent(out) :: L_val(:)

  integer :: ig, ix, iz, idx, vtype, nx_vg, nz_vg
  integer :: n_smooth, n_smooth_nnz, n_neighbors, j, irow, jdx
  type(Tvelocity_grid2d), pointer :: vg

  n_smooth = 0
  n_smooth_nnz = 0
  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx <= 0) cycle
        n_neighbors = 0
        if (ix > 1) then
          if (vg%inv_idx(ix-1, iz) > 0) n_neighbors = n_neighbors + 1
        end if
        if (ix < nx_vg) then
          if (vg%inv_idx(ix+1, iz) > 0) n_neighbors = n_neighbors + 1
        end if
        if (iz > 1) then
          if (vg%inv_idx(ix, iz-1) > 0) n_neighbors = n_neighbors + 1
        end if
        if (iz < nz_vg) then
          if (vg%inv_idx(ix, iz+1) > 0) n_neighbors = n_neighbors + 1
        end if
        if (n_neighbors > 0) then
          n_smooth = n_smooth + 1
          n_smooth_nnz = n_smooth_nnz + 1 + n_neighbors
        end if
      end do
    end do
  end do

  n_L_nnz = n_smooth_nnz
  allocate(L_row(n_L_nnz), L_col(n_L_nnz), L_val(n_L_nnz))

  j = 0
  irow = 0
  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx <= 0) cycle
        n_neighbors = 0
        if (ix > 1) then
          if (vg%inv_idx(ix-1, iz) > 0) n_neighbors = n_neighbors + 1
        end if
        if (ix < nx_vg) then
          if (vg%inv_idx(ix+1, iz) > 0) n_neighbors = n_neighbors + 1
        end if
        if (iz > 1) then
          if (vg%inv_idx(ix, iz-1) > 0) n_neighbors = n_neighbors + 1
        end if
        if (iz < nz_vg) then
          if (vg%inv_idx(ix, iz+1) > 0) n_neighbors = n_neighbors + 1
        end if
        if (n_neighbors > 0) then
          irow = irow + 1
          j = j + 1
          L_row(j) = irow
          L_col(j) = idx
          L_val(j) = -dble(n_neighbors)
          if (ix > 1) then
            jdx = vg%inv_idx(ix-1, iz)
            if (jdx > 0) then
              j = j + 1
              L_row(j) = irow
              L_col(j) = jdx
              L_val(j) = 1.0_dp
            end if
          end if
          if (ix < nx_vg) then
            jdx = vg%inv_idx(ix+1, iz)
            if (jdx > 0) then
              j = j + 1
              L_row(j) = irow
              L_col(j) = jdx
              L_val(j) = 1.0_dp
            end if
          end if
          if (iz > 1) then
            jdx = vg%inv_idx(ix, iz-1)
            if (jdx > 0) then
              j = j + 1
              L_row(j) = irow
              L_col(j) = jdx
              L_val(j) = 1.0_dp
            end if
          end if
          if (iz < nz_vg) then
            jdx = vg%inv_idx(ix, iz+1)
            if (jdx > 0) then
              j = j + 1
              L_row(j) = irow
              L_col(j) = jdx
              L_val(j) = 1.0_dp
            end if
          end if
        end if
      end do
    end do
  end do

  n_L_nnz = j

  write(*,'(A,I8,A,I8)') '  Laplacian smoothing: n_smooth=', n_smooth, ' nnz=', n_L_nnz

end subroutine build_laplacian_smoothing


subroutine cgnr_laplacian_solver(ndata, nparm, nnz, gindx, gval, row_ptr, b, &
  lam_vel, smooth_scale, n_L_nnz, L_row, L_col, L_val, x, max_iter)
  implicit none
  integer, intent(in) :: ndata, nparm, nnz, max_iter
  integer, intent(in) :: gindx(nnz), row_ptr(ndata+1)
  real(kind=dp), intent(in) :: gval(nnz), b(ndata)
  real(kind=dp), intent(in) :: lam_vel, smooth_scale
  integer, intent(in) :: n_L_nnz, L_row(n_L_nnz), L_col(n_L_nnz)
  real(kind=dp), intent(in) :: L_val(n_L_nnz)
  real(kind=dp), intent(inout) :: x(nparm)

  real(kind=dp), allocatable :: r(:), p(:), Ap(:)
  real(kind=dp) :: rnorm, pAp, alpha, beta, rnorm_new
  real(kind=dp) :: rhs_norm
  integer :: iter, i, n_L_rows
  real(kind=dp), allocatable :: Lp(:), LtLp(:)

  n_L_rows = 0
  do i = 1, n_L_nnz
    if (L_row(i) > n_L_rows) n_L_rows = L_row(i)
  end do

  allocate(r(nparm), p(nparm), Ap(nparm))
  allocate(Lp(n_L_rows), LtLp(nparm))

  call sparse_ATx(ndata, nparm, nnz, gindx, gval, row_ptr, b, r)

  r = r + lam_vel**2 * x + smooth_scale**2 * matmul_LtL(nparm, n_L_nnz, L_row, L_col, L_val, x)

  rhs_norm = sqrt(dot_product(r, r))
  if (rhs_norm < 1.0e-30_dp) then
    x = 0.0_dp
    deallocate(r, p, Ap, Lp, LtLp)
    return
  end if

  p = r
  rnorm = dot_product(r, r)

  do iter = 1, max_iter
    call compute_Mp_laplacian(ndata, nparm, nnz, gindx, gval, row_ptr, &
      lam_vel, smooth_scale, n_L_nnz, L_row, L_col, L_val, p, Ap)

    pAp = dot_product(p, Ap)
    if (abs(pAp) < 1.0e-30_dp) exit

    alpha = rnorm / pAp
    x = x + alpha * p
    r = r - alpha * Ap

    rnorm_new = dot_product(r, r)

    if (mod(iter, 50) == 0 .or. iter == max_iter) then
      write(*,'(A,I5,A,F12.6,A,F12.6)') '    CGNR iter=', iter, &
        ' |r|/|rhs|=', sqrt(rnorm_new)/rhs_norm, ' |x|=', sqrt(dot_product(x,x))
    end if

    if (sqrt(rnorm_new) / rhs_norm < 1.0e-3_dp) exit

    beta = rnorm_new / rnorm
    p = r + beta * p
    rnorm = rnorm_new
  end do

  write(*,'(A,I5,A,F12.6,A,F12.6)') '    CGNR: ', iter, &
    ' iterations, |r|/|rhs|=', sqrt(rnorm)/rhs_norm, ' |x|=', sqrt(dot_product(x,x))

  deallocate(r, p, Ap, Lp, LtLp)

contains

  function matmul_LtL(n, nnz_L, L_r, L_c, L_v, vec) result(y)
    integer, intent(in) :: n, nnz_L
    integer, intent(in) :: L_r(nnz_L), L_c(nnz_L)
    real(kind=dp), intent(in) :: L_v(nnz_L), vec(n)
    real(kind=dp) :: y(n), Lvec(size(L_r))

    Lvec = 0.0_dp
    do i = 1, nnz_L
      Lvec(L_r(i)) = Lvec(L_r(i)) + L_v(i) * vec(L_c(i))
    end do

    y = 0.0_dp
    do i = 1, nnz_L
      y(L_c(i)) = y(L_c(i)) + L_v(i) * Lvec(L_r(i))
    end do
  end function matmul_LtL

  subroutine compute_Mp_laplacian(nd, np, nz, gi, gv, rp, lam, ss, nlz, lr, lc, lv, pvec, Mpvec)
    integer, intent(in) :: nd, np, nz, nlz
    integer, intent(in) :: gi(nz), rp(nd+1), lr(nlz), lc(nlz)
    real(kind=dp), intent(in) :: gv(nz), lam, ss, lv(nlz), pvec(np)
    real(kind=dp), intent(out) :: Mpvec(np)

    real(kind=dp), allocatable :: Ap_vec(:), AtAp_vec(:), Lp_vec(:), LtLp_vec(:)
    integer :: n_lr

    n_lr = 0
    do i = 1, nlz
      if (lr(i) > n_lr) n_lr = lr(i)
    end do

    allocate(Ap_vec(nd), AtAp_vec(np), Lp_vec(n_lr), LtLp_vec(np))

    call sparse_Ax(nd, np, nz, gi, gv, rp, pvec, Ap_vec)
    call sparse_ATx(nd, np, nz, gi, gv, rp, Ap_vec, AtAp_vec)

    Lp_vec = 0.0_dp
    do i = 1, nlz
      Lp_vec(lr(i)) = Lp_vec(lr(i)) + lv(i) * pvec(lc(i))
    end do

    LtLp_vec = 0.0_dp
    do i = 1, nlz
      LtLp_vec(lc(i)) = LtLp_vec(lc(i)) + lv(i) * Lp_vec(lr(i))
    end do

    Mpvec = AtAp_vec + lam**2 * pvec + ss**2 * LtLp_vec

    deallocate(Ap_vec, AtAp_vec, Lp_vec, LtLp_vec)
  end subroutine compute_Mp_laplacian

end subroutine cgnr_laplacian_solver


subroutine build_gaussian_smoothing(mvscale, nparm, corr_h, corr_v, &
  n_smooth, nnz_Rh, Rh_indx, Rh_val, Rh_rowptr, &
  nnz_Rv, Rv_indx, Rv_val, Rv_rowptr, smooth_to_parm)

  implicit none
  integer, intent(in) :: nparm
  real(kind=dp), intent(in) :: mvscale(nparm), corr_h, corr_v
  integer, intent(out) :: n_smooth, nnz_Rh, nnz_Rv
  integer, allocatable, intent(out) :: Rh_indx(:), Rh_rowptr(:), Rv_indx(:), Rv_rowptr(:)
  real(kind=dp), allocatable, intent(out) :: Rh_val(:), Rv_val(:)
  integer, allocatable, intent(out) :: smooth_to_parm(:)

  integer :: ig, ix, iz, jx, jz, idx, jdx, vtype, nx_vg, nz_vg
  real(kind=dp) :: dx, dz, beta, beta_sum
  type(Tvelocity_grid2d), pointer :: vg

  n_smooth = 0
  nnz_Rh = 0
  nnz_Rv = 0

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz

    do iz = 2, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx <= 0) cycle

        n_smooth = n_smooth + 1

        do jx = 1, nx_vg
          if (jx == ix) cycle
          jdx = vg%inv_idx(jx, iz)
          if (jdx <= 0) cycle
          dx = abs(vg%x(jx) - vg%x(ix))
          if (dx <= corr_h) nnz_Rh = nnz_Rh + 1
        end do
        nnz_Rh = nnz_Rh + 1

        do jz = 2, nz_vg
          if (jz == iz) cycle
          jdx = vg%inv_idx(ix, jz)
          if (jdx <= 0) cycle
          dz = abs(vg%z(jz) - vg%z(iz))
          if (dz <= corr_v) nnz_Rv = nnz_Rv + 1
        end do
        nnz_Rv = nnz_Rv + 1
      end do
    end do
  end do

  allocate(Rh_indx(nnz_Rh), Rh_val(nnz_Rh), Rh_rowptr(n_smooth+1))
  allocate(Rv_indx(nnz_Rv), Rv_val(nnz_Rv), Rv_rowptr(n_smooth+1))
  allocate(smooth_to_parm(n_smooth))

  n_smooth = 0
  nnz_Rh = 0
  nnz_Rv = 0
  Rh_rowptr(1) = 1
  Rv_rowptr(1) = 1

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz

    do iz = 2, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx <= 0) cycle

        n_smooth = n_smooth + 1
        smooth_to_parm(n_smooth) = idx

        beta_sum = 0.0_dp
        do jx = 1, nx_vg
          if (jx == ix) cycle
          jdx = vg%inv_idx(jx, iz)
          if (jdx <= 0) cycle
          dx = abs(vg%x(jx) - vg%x(ix))
          if (dx <= corr_h) then
            beta = exp(-dx*dx / (corr_h*corr_h))
            nnz_Rh = nnz_Rh + 1
            Rh_indx(nnz_Rh) = jdx
            Rh_val(nnz_Rh) = beta * mvscale(jdx)
            beta_sum = beta_sum + beta
          end if
        end do
        nnz_Rh = nnz_Rh + 1
        Rh_indx(nnz_Rh) = idx
        Rh_val(nnz_Rh) = -beta_sum * mvscale(idx)
        Rh_rowptr(n_smooth+1) = nnz_Rh + 1

        beta_sum = 0.0_dp
        do jz = 2, nz_vg
          if (jz == iz) cycle
          jdx = vg%inv_idx(ix, jz)
          if (jdx <= 0) cycle
          dz = abs(vg%z(jz) - vg%z(iz))
          if (dz <= corr_v) then
            beta = exp(-dz*dz / (corr_v*corr_v))
            nnz_Rv = nnz_Rv + 1
            Rv_indx(nnz_Rv) = jdx
            Rv_val(nnz_Rv) = beta * mvscale(jdx)
            beta_sum = beta_sum + beta
          end if
        end do
        nnz_Rv = nnz_Rv + 1
        Rv_indx(nnz_Rv) = idx
        Rv_val(nnz_Rv) = -beta_sum * mvscale(idx)
        Rv_rowptr(n_smooth+1) = nnz_Rv + 1

      end do
    end do
  end do

  write(*,'(A,I6,A,I8,A,I8)') '  Gaussian smoothing: n_smooth=', n_smooth, &
    ' nnz_Rh=', nnz_Rh, ' nnz_Rv=', nnz_Rv

end subroutine build_gaussian_smoothing


subroutine lsqr_gaussian_solver(ndata, nparm, nnz, gindx, gval, &
  row_ptr, b, wsh, wsv, wdv, damp_area, &
  n_smooth, nnz_Rh, Rh_indx, Rh_val, Rh_rowptr, &
  nnz_Rv, Rv_indx, Rv_val, Rv_rowptr, &
  smooth_to_parm, m_sum, x, max_iter)

  implicit none
  integer, intent(in) :: ndata, nparm, nnz, max_iter
  integer, intent(in) :: gindx(nnz), row_ptr(ndata+1)
  real(kind=dp), intent(in) :: gval(nnz), b(ndata)
  real(kind=dp), intent(in) :: wsh, wsv, wdv
  real(kind=dp), intent(in) :: damp_area(nparm)
  integer, intent(in) :: n_smooth, nnz_Rh, nnz_Rv
  integer, intent(in) :: Rh_indx(nnz_Rh), Rh_rowptr(n_smooth+1)
  real(kind=dp), intent(in) :: Rh_val(nnz_Rh)
  integer, intent(in) :: Rv_indx(nnz_Rv), Rv_rowptr(n_smooth+1)
  real(kind=dp), intent(in) :: Rv_val(nnz_Rv)
  integer, intent(in) :: smooth_to_parm(n_smooth)
  real(kind=dp), intent(in) :: m_sum(nparm)
  real(kind=dp), intent(inout) :: x(nparm)

  integer :: ndata_aug, nnz_aug, irow, innz, i, k, k1, k2, idx
  integer, allocatable :: gindx_aug(:), row_ptr_aug(:)
  real(kind=dp), allocatable :: gval_aug(:), b_aug(:)
  real(kind=dp), allocatable :: Rh_msum(:), Rv_msum(:)
  real(kind=dp), allocatable :: u(:), v(:), w(:), Av(:), Atu(:), z(:)
  real(kind=dp) :: beta, alpha, c_bar, s_bar, phi_bar, rho_bar
  real(kind=dp) :: rho, theta, c, s, phi, rnorm, bnorm
  real(kind=dp) :: Anorm2, Anorm, Arnorm, test2, lsqr_atol
  integer :: iter

  allocate(Rh_msum(n_smooth), Rv_msum(n_smooth))
  call sparse_Ax(n_smooth, nparm, nnz_Rh, Rh_indx, Rh_val, Rh_rowptr, m_sum, Rh_msum)
  call sparse_Ax(n_smooth, nparm, nnz_Rv, Rv_indx, Rv_val, Rv_rowptr, m_sum, Rv_msum)

  ndata_aug = ndata + n_smooth + n_smooth + n_smooth
  nnz_aug = nnz + nnz_Rh + nnz_Rv + n_smooth

  allocate(gindx_aug(nnz_aug), gval_aug(nnz_aug))
  allocate(b_aug(ndata_aug), row_ptr_aug(ndata_aug+1))

  gindx_aug(1:nnz) = gindx(1:nnz)
  gval_aug(1:nnz) = gval(1:nnz)
  b_aug(1:ndata) = b(1:ndata)
  row_ptr_aug(1:ndata+1) = row_ptr(1:ndata+1)

  innz = nnz
  irow = ndata

  do i = 1, n_smooth
    irow = irow + 1
    b_aug(irow) = -wsh * Rh_msum(i)
    row_ptr_aug(irow) = innz + 1
    k1 = Rh_rowptr(i)
    k2 = Rh_rowptr(i+1) - 1
    do k = k1, k2
      innz = innz + 1
      gindx_aug(innz) = Rh_indx(k)
      gval_aug(innz) = wsh * Rh_val(k)
    end do
  end do

  do i = 1, n_smooth
    irow = irow + 1
    b_aug(irow) = -wsv * Rv_msum(i)
    row_ptr_aug(irow) = innz + 1
    k1 = Rv_rowptr(i)
    k2 = Rv_rowptr(i+1) - 1
    do k = k1, k2
      innz = innz + 1
      gindx_aug(innz) = Rv_indx(k)
      gval_aug(innz) = wsv * Rv_val(k)
    end do
  end do

  do i = 1, n_smooth
    irow = irow + 1
    idx = smooth_to_parm(i)
    b_aug(irow) = -wdv * damp_area(idx) * m_sum(idx)
    row_ptr_aug(irow) = innz + 1
    innz = innz + 1
    gindx_aug(innz) = idx
    gval_aug(innz) = wdv * damp_area(idx)
  end do

  row_ptr_aug(irow+1) = innz + 1

  write(*,'(A,I6,A,I10,A,I6)') '    augmented system: ndata_aug=', irow, &
    ' nnz_aug=', innz, ' n_smooth=', n_smooth

  deallocate(Rh_msum, Rv_msum)

  bnorm = sqrt(dot_product(b_aug(1:irow), b_aug(1:irow)))
  if (bnorm < 1.0e-30_dp) then
    x = 0.0_dp
    deallocate(gindx_aug, gval_aug, b_aug, row_ptr_aug)
    return
  end if

  allocate(u(irow), v(nparm), w(nparm), Av(irow), Atu(nparm), z(nparm))

  x = 0.0_dp

  u(1:irow) = b_aug(1:irow)
  beta = sqrt(dot_product(u, u))
  if (beta > 0.0_dp) then
    u = u / beta
  else
    x = 0.0_dp
    deallocate(u, v, w, Av, Atu, z, gindx_aug, gval_aug, b_aug, row_ptr_aug)
    return
  end if

  call sparse_ATx(irow, nparm, innz, gindx_aug, gval_aug, row_ptr_aug, u, v)
  alpha = sqrt(dot_product(v, v))
  if (alpha > 0.0_dp) then
    v = v / alpha
  else
    x = 0.0_dp
    deallocate(u, v, w, Av, Atu, z, gindx_aug, gval_aug, b_aug, row_ptr_aug)
    return
  end if

  w = v
  c_bar = 1.0_dp
  s_bar = 0.0_dp
  phi_bar = beta
  rho_bar = alpha
  Anorm2 = alpha**2 + beta**2
  lsqr_atol = 1.0e-8_dp

  do iter = 1, max_iter
    call sparse_Ax(irow, nparm, innz, gindx_aug, gval_aug, row_ptr_aug, v, Av)
    u = Av - alpha * u
    beta = sqrt(dot_product(u, u))
    if (beta > 0.0_dp) u = u / beta

    call sparse_ATx(irow, nparm, innz, gindx_aug, gval_aug, row_ptr_aug, u, Atu)
    v = Atu - beta * v
    alpha = sqrt(dot_product(v, v))
    if (alpha > 0.0_dp) v = v / alpha

    Anorm2 = Anorm2 + alpha**2 + beta**2

    rho = sqrt(rho_bar**2 + beta**2)
    c = rho_bar / rho
    s = beta / rho
    theta = s * alpha
    rho_bar = -c * alpha
    phi = c * phi_bar
    phi_bar = s * phi_bar

    z = w
    w = v - (theta / rho) * w
    x = x + (phi / rho) * z

    rnorm = abs(phi_bar)
    Arnorm = abs(phi_bar * alpha * c)
    Anorm = sqrt(Anorm2)

    if (Anorm > 0.0_dp .and. rnorm > 0.0_dp) then
      test2 = Arnorm / (Anorm * rnorm)
    else
      test2 = 0.0_dp
    end if

    if (mod(iter, 20) == 0) then
      write(*,'(A,I4,A,F12.6,A,F12.6,A,F12.6)') '    LSQR iter=', iter, &
        ' |r|/|b|=', rnorm/bnorm, ' test2=', test2, ' |x|=', sqrt(dot_product(x,x))
    end if

    if (test2 < lsqr_atol) exit
  end do

  write(*,'(A,I4,A,F12.6,A,F12.6,A,F12.6)') '    LSQR: ', iter, &
    ' iterations, |r|/|b|=', rnorm/bnorm, ' test2=', test2, ' |x|=', sqrt(dot_product(x,x))

  deallocate(u, v, w, Av, Atu, z)
  deallocate(gindx_aug, gval_aug, b_aug, row_ptr_aug)

end subroutine lsqr_gaussian_solver


subroutine cgnr_smooth_solver(ndata, nparm, nnz, gindx, gval, &
  row_ptr, b, wsv, wdv, &
  n_smooth, nnz_Rh, Rh_indx, Rh_val, Rh_rowptr, &
  nnz_Rv, Rv_indx, Rv_val, Rv_rowptr, &
  smooth_to_parm, m_sum, x, max_iter)

  implicit none
  integer, intent(in) :: ndata, nparm, nnz, max_iter
  integer, intent(in) :: gindx(nnz), row_ptr(ndata+1)
  real(kind=dp), intent(in) :: gval(nnz), b(ndata)
  real(kind=dp), intent(in) :: wsv, wdv
  integer, intent(in) :: n_smooth, nnz_Rh, nnz_Rv
  integer, intent(in) :: Rh_indx(nnz_Rh), Rh_rowptr(n_smooth+1)
  real(kind=dp), intent(in) :: Rh_val(nnz_Rh)
  integer, intent(in) :: Rv_indx(nnz_Rv), Rv_rowptr(n_smooth+1)
  real(kind=dp), intent(in) :: Rv_val(nnz_Rv)
  integer, intent(in) :: smooth_to_parm(n_smooth)
  real(kind=dp), intent(in) :: m_sum(nparm)
  real(kind=dp), intent(inout) :: x(nparm)

  real(kind=dp), allocatable :: rhs(:), r(:), p(:), Mp(:)
  real(kind=dp), allocatable :: Ax(:), ATb(:), Rh_p(:), RhT_Rh_p(:)
  real(kind=dp), allocatable :: Rv_p(:), RvT_Rv_p(:), tmp(:)
  real(kind=dp) :: rz_old, rz_new, pMp, alpha_cg, beta_cg
  real(kind=dp) :: rnorm, rhsnorm
  integer :: iter

  allocate(rhs(nparm), r(nparm), p(nparm), Mp(nparm))
  allocate(Ax(ndata), ATb(nparm))
  allocate(Rh_p(n_smooth), RhT_Rh_p(nparm))
  allocate(Rv_p(n_smooth), RvT_Rv_p(nparm))

  call sparse_ATx(ndata, nparm, nnz, gindx, gval, row_ptr, b, ATb)
  rhs = ATb

  if (wsv > 0.0_dp) then
    allocate(tmp(nparm))
    call sparse_Ax(n_smooth, nparm, nnz_Rh, Rh_indx, Rh_val, Rh_rowptr, m_sum, Rh_p)
    call sparse_ATx(n_smooth, nparm, nnz_Rh, Rh_indx, Rh_val, Rh_rowptr, Rh_p, tmp)
    rhs = rhs - wsv * wsv * tmp
    call sparse_Ax(n_smooth, nparm, nnz_Rv, Rv_indx, Rv_val, Rv_rowptr, m_sum, Rv_p)
    call sparse_ATx(n_smooth, nparm, nnz_Rv, Rv_indx, Rv_val, Rv_rowptr, Rv_p, tmp)
    rhs = rhs - wsv * wsv * tmp
    deallocate(tmp)
  end if

  if (wdv > 0.0_dp) then
    rhs = rhs - wdv * wdv * m_sum
  end if

  x = 0.0_dp
  r = rhs

  call compute_Mp(x)
  r = r - Mp

  p = r
  rz_old = dot_product(r, r)
  rhsnorm = sqrt(dot_product(rhs, rhs))

  if (rhsnorm < 1.0e-30_dp) then
    x = 0.0_dp
    deallocate(rhs, r, p, Mp, Ax, ATb, Rh_p, RhT_Rh_p, Rv_p, RvT_Rv_p)
    return
  end if

  do iter = 1, max_iter
    call compute_Mp(p)
    pMp = dot_product(p, Mp)
    if (abs(pMp) < 1.0e-30_dp) exit

    alpha_cg = rz_old / pMp
    x = x + alpha_cg * p
    r = r - alpha_cg * Mp

    rz_new = dot_product(r, r)
    rnorm = sqrt(rz_new)

    if (mod(iter, 50) == 0) then
      write(*,'(A,I5,A,F12.6,A,F12.6)') '    CG iter=', iter, &
        ' |r|=', rnorm, ' |r|/|rhs|=', rnorm/rhsnorm
    end if

    if (rnorm / rhsnorm < 1.0e-4_dp) exit

    beta_cg = rz_new / rz_old
    p = r + beta_cg * p
    rz_old = rz_new
  end do

  write(*,'(A,I5,A,F12.6,A,F12.6)') '    CGNR_smooth: ', iter, &
    ' iterations, |r|/|rhs|=', rnorm/rhsnorm, ' |x|=', sqrt(dot_product(x,x))

  deallocate(rhs, r, p, Mp, Ax, ATb, Rh_p, RhT_Rh_p, Rv_p, RvT_Rv_p)

contains

  subroutine compute_Mp(pvec)
    real(kind=dp), intent(in) :: pvec(nparm)
    real(kind=dp), allocatable :: tmp2(:)

    call sparse_Ax(ndata, nparm, nnz, gindx, gval, row_ptr, pvec, Ax)
    call sparse_ATx(ndata, nparm, nnz, gindx, gval, row_ptr, Ax, Mp)

    if (wsv > 0.0_dp) then
      call sparse_Ax(n_smooth, nparm, nnz_Rh, Rh_indx, Rh_val, Rh_rowptr, pvec, Rh_p)
      call sparse_ATx(n_smooth, nparm, nnz_Rh, Rh_indx, Rh_val, Rh_rowptr, Rh_p, RhT_Rh_p)
      Mp = Mp + wsv * wsv * RhT_Rh_p

      call sparse_Ax(n_smooth, nparm, nnz_Rv, Rv_indx, Rv_val, Rv_rowptr, pvec, Rv_p)
      call sparse_ATx(n_smooth, nparm, nnz_Rv, Rv_indx, Rv_val, Rv_rowptr, Rv_p, RvT_Rv_p)
      Mp = Mp + wsv * wsv * RvT_Rv_p
    end if

    if (wdv > 0.0_dp) then
      Mp = Mp + wdv * wdv * pvec
    end if

  end subroutine compute_Mp

end subroutine cgnr_smooth_solver


subroutine update_model_jumping(m_sum, mvscale)
  implicit none
  real(kind=dp), intent(in) :: m_sum(n_inv_parms)
  real(kind=dp), intent(in) :: mvscale(n_inv_parms)

  integer :: ig, ix, iz, idx
  integer :: nx_vg, nz_vg, vtype
  integer :: n_updated, n_skipped
  real(kind=dp) :: v_init, new_v, max_change
  type(Tvelocity_grid2d), pointer :: vg

  n_updated = 0
  n_skipped = 0
  max_change = 0.0_dp

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz

    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0) then
          v_init = vg%velocity_init(vtype, ix, iz)
          new_v = v_init + mvscale(idx) * m_sum(idx)
          if (new_v < 2.0_dp) new_v = 2.0_dp
          if (new_v > 10.0_dp) new_v = 10.0_dp
          vg%velocity(vtype, ix, iz) = new_v
          if (abs(new_v - v_init) > max_change) max_change = abs(new_v - v_init)
          if (abs(new_v - v_init) > 1.0e-6_dp) n_updated = n_updated + 1
        end if
      end do
    end do
  end do

  write(*,'(A,I6,A,I6,A,F12.6)') '  update_model_jumping: updated=', n_updated, &
    ' skipped_water=', n_skipped, ' max_dv=', max_change

end subroutine update_model_jumping


subroutine lsqr_augmented_solver(ndata_aug, nparm, nnz_aug, &
  aug_col, aug_val, aug_row, aug_rhs, x, max_iter)
  implicit none
  integer, intent(in) :: ndata_aug, nparm, nnz_aug, max_iter
  integer, intent(in) :: aug_col(nnz_aug), aug_row(nnz_aug)
  real(kind=dp), intent(in) :: aug_val(nnz_aug), aug_rhs(ndata_aug)
  real(kind=dp), intent(inout) :: x(nparm)

  integer, allocatable :: row_ptr(:), row_count(:), perm(:)
  integer, allocatable :: gindx_aug(:)
  real(kind=dp), allocatable :: gval_aug(:), b_aug(:)
  real(kind=dp), allocatable :: u(:), v(:), w(:), Av(:), Atu(:), z(:)
  real(kind=dp) :: beta, alpha, c_bar, s_bar, phi_bar, rho_bar
  real(kind=dp) :: rho, theta, c, s, phi, rnorm, bnorm
  real(kind=dp) :: arnorm_v, anorm_v, bnorm_v, dnorm_v, test2_v, atol_v
  integer :: iter, i, k, irow, innz, current_row

  allocate(row_ptr(ndata_aug+1), row_count(ndata_aug))
  allocate(gindx_aug(nnz_aug), gval_aug(nnz_aug), b_aug(ndata_aug))
  allocate(perm(nnz_aug))

  row_count = 0
  do i = 1, nnz_aug
    row_count(aug_row(i)) = row_count(aug_row(i)) + 1
  end do

  row_ptr(1) = 1
  do irow = 1, ndata_aug
    row_ptr(irow+1) = row_ptr(irow) + row_count(irow)
  end do

  row_count = 0
  do i = 1, nnz_aug
    irow = aug_row(i)
    innz = row_ptr(irow) + row_count(irow)
    gindx_aug(innz) = aug_col(i)
    gval_aug(innz) = aug_val(i)
    row_count(irow) = row_count(irow) + 1
  end do

  deallocate(row_count, perm)

  b_aug(1:ndata_aug) = aug_rhs(1:ndata_aug)

  bnorm = sqrt(dot_product(b_aug, b_aug))
  if (bnorm < 1.0e-30_dp) then
    x = 0.0_dp
    deallocate(row_ptr, gindx_aug, gval_aug, b_aug)
    return
  end if

  allocate(u(ndata_aug), v(nparm), w(nparm), Av(ndata_aug), Atu(nparm), z(nparm))

  x = 0.0_dp

  u = b_aug
  beta = sqrt(dot_product(u, u))
  if (beta > 0.0_dp) then
    u = u / beta
  else
    x = 0.0_dp
    deallocate(u, v, w, Av, Atu, z, row_ptr, gindx_aug, gval_aug, b_aug)
    return
  end if

  call sparse_ATx(ndata_aug, nparm, innz, gindx_aug, gval_aug, row_ptr, u, v)
  alpha = sqrt(dot_product(v, v))
  if (alpha > 0.0_dp) then
    v = v / alpha
  else
    x = 0.0_dp
    deallocate(u, v, w, Av, Atu, z, row_ptr, gindx_aug, gval_aug, b_aug)
    return
  end if

  w = v
  c_bar = 1.0_dp
  s_bar = 0.0_dp
  phi_bar = beta
  rho_bar = alpha
  bnorm_v = alpha*alpha + beta*beta
  dnorm_v = 1.0_dp / (rho_bar*rho_bar)
  atol_v = 1.0e-3_dp

  do iter = 1, max_iter
    call sparse_Ax(ndata_aug, nparm, innz, gindx_aug, gval_aug, row_ptr, v, Av)
    u = Av - alpha * u
    beta = sqrt(dot_product(u, u))
    if (beta > 0.0_dp) u = u / beta

    call sparse_ATx(ndata_aug, nparm, innz, gindx_aug, gval_aug, row_ptr, u, Atu)
    v = Atu - beta * v
    alpha = sqrt(dot_product(v, v))
    if (alpha > 0.0_dp) v = v / alpha

    bnorm_v = bnorm_v + alpha*alpha + beta*beta

    rho = sqrt(rho_bar**2 + beta**2)
    c = rho_bar / rho
    s = beta / rho
    theta = s * alpha
    rho_bar = -c * alpha
    phi = c * phi_bar
    phi_bar = s * phi_bar

    z = w
    w = v - (theta / rho) * w
    x = x + (phi / rho) * z

    dnorm_v = dnorm_v + dot_product(z, z) / (rho*rho)

    rnorm = abs(phi_bar)
    anorm_v = sqrt(bnorm_v)
    arnorm_v = phi_bar * alpha * abs(c)
    test2_v = arnorm_v / (anorm_v * rnorm + 1.0e-30_dp)

    if (mod(iter, 100) == 0) then
      write(*,'(A,I5,A,F12.6,A,ES10.3,A,ES10.3)') '    LSQR iter=', iter, &
        ' |r|/|b|=', rnorm/bnorm, ' test2=', test2_v, ' |x|=', sqrt(dot_product(x,x))
    end if

    if (test2_v < atol_v) exit
    if (1.0_dp + test2_v <= 1.0_dp) exit
  end do

  write(*,'(A,I5,A,F12.6,A,ES10.3,A,ES10.3)') '    LSQR: ', iter, &
    ' iterations, |r|/|b|=', rnorm/bnorm, ' test2=', test2_v, ' |x|=', sqrt(dot_product(x,x))

  deallocate(u, v, w, Av, Atu, z)
  deallocate(row_ptr, gindx_aug, gval_aug, b_aug)

end subroutine lsqr_augmented_solver


subroutine lsqr_damped_solver(ndata, nparm, nnz, gindx, gval, row_ptr, b, damping, x, max_iter)
  implicit none
  integer, intent(in) :: ndata, nparm, nnz, max_iter
  integer, intent(in) :: gindx(nnz), row_ptr(ndata+1)
  real(kind=dp), intent(in) :: gval(nnz), b(ndata), damping
  real(kind=dp), intent(inout) :: x(nparm)

  integer :: ndata_aug, nnz_aug, irow, innz, i, k, k1, k2, j
  integer, allocatable :: gindx_aug(:), row_ptr_aug(:)
  real(kind=dp), allocatable :: gval_aug(:), b_aug(:)
  real(kind=dp), allocatable :: u(:), v(:), w(:), Av(:), Atu(:), z(:)
  real(kind=dp) :: beta, alpha, c_bar, s_bar, phi_bar, rho_bar
  real(kind=dp) :: rho, theta, c, s, phi, rnorm, bnorm
  integer :: iter

  ndata_aug = ndata + nparm
  nnz_aug = nnz + nparm

  allocate(gindx_aug(nnz_aug), gval_aug(nnz_aug))
  allocate(b_aug(ndata_aug), row_ptr_aug(ndata_aug+1))

  gindx_aug(1:nnz) = gindx(1:nnz)
  gval_aug(1:nnz) = gval(1:nnz)
  b_aug(1:ndata) = b(1:ndata)
  row_ptr_aug(1:ndata+1) = row_ptr(1:ndata+1)

  innz = nnz
  do i = 1, nparm
    irow = ndata + i
    b_aug(irow) = 0.0_dp
    row_ptr_aug(irow) = innz + 1
    innz = innz + 1
    gindx_aug(innz) = i
    gval_aug(innz) = damping
  end do
  row_ptr_aug(ndata_aug + 1) = innz + 1

  bnorm = sqrt(dot_product(b_aug, b_aug))
  if (bnorm < 1.0e-30_dp) then
    x = 0.0_dp
    deallocate(gindx_aug, gval_aug, b_aug, row_ptr_aug)
    return
  end if

  allocate(u(ndata_aug), v(nparm), w(nparm), Av(ndata_aug), Atu(nparm), z(nparm))

  x = 0.0_dp

  u = b_aug
  beta = sqrt(dot_product(u, u))
  if (beta > 0.0_dp) then
    u = u / beta
  else
    x = 0.0_dp
    deallocate(u, v, w, Av, Atu, z, gindx_aug, gval_aug, b_aug, row_ptr_aug)
    return
  end if

  call sparse_ATx(ndata_aug, nparm, nnz_aug, gindx_aug, gval_aug, row_ptr_aug, u, v)
  alpha = sqrt(dot_product(v, v))
  if (alpha > 0.0_dp) then
    v = v / alpha
  else
    x = 0.0_dp
    deallocate(u, v, w, Av, Atu, z, gindx_aug, gval_aug, b_aug, row_ptr_aug)
    return
  end if

  w = v
  c_bar = 1.0_dp
  s_bar = 0.0_dp
  phi_bar = beta
  rho_bar = alpha

  do iter = 1, max_iter
    call sparse_Ax(ndata_aug, nparm, nnz_aug, gindx_aug, gval_aug, row_ptr_aug, v, Av)
    u = Av - alpha * u
    beta = sqrt(dot_product(u, u))
    if (beta > 0.0_dp) u = u / beta

    call sparse_ATx(ndata_aug, nparm, nnz_aug, gindx_aug, gval_aug, row_ptr_aug, u, Atu)
    v = Atu - beta * v
    alpha = sqrt(dot_product(v, v))
    if (alpha > 0.0_dp) v = v / alpha

    rho = sqrt(rho_bar**2 + beta**2)
    c = rho_bar / rho
    s = beta / rho
    theta = s * alpha
    rho_bar = -c * alpha
    phi = c * phi_bar
    phi_bar = s * phi_bar

    z = w
    w = v - (theta / rho) * w
    x = x + (phi / rho) * z

    rnorm = abs(phi_bar)

    if (mod(iter, 50) == 0) then
      write(*,'(A,I4,A,F12.6)') '    LSQR iter=', iter, ' |r|=', rnorm
    end if

    if (rnorm / bnorm < 1.0e-3_dp) exit
  end do

  write(*,'(A,I4,A,F12.6)') '    LSQR: ', iter, ' iterations, |r|=', rnorm

  deallocate(u, v, w, Av, Atu, z)
  deallocate(gindx_aug, gval_aug, b_aug, row_ptr_aug)

end subroutine lsqr_damped_solver


subroutine sparse_Ax(ndata, nparm, nnz, gindx, gval, row_ptr, x, y)
  implicit none
  integer, intent(in) :: ndata, nparm, nnz
  integer, intent(in) :: gindx(nnz), row_ptr(ndata+1)
  real(kind=dp), intent(in) :: gval(nnz), x(nparm)
  real(kind=dp), intent(out) :: y(ndata)

  integer :: i, k, k1, k2, j

  do i = 1, ndata
    y(i) = 0.0_dp
    k1 = row_ptr(i)
    k2 = row_ptr(i+1) - 1
    do k = k1, k2
      j = gindx(k)
      if (j >= 1 .and. j <= nparm) then
        y(i) = y(i) + gval(k) * x(j)
      end if
    end do
  end do

end subroutine sparse_Ax


subroutine sparse_ATx(ndata, nparm, nnz, gindx, gval, row_ptr, y, x)
  implicit none
  integer, intent(in) :: ndata, nparm, nnz
  integer, intent(in) :: gindx(nnz), row_ptr(ndata+1)
  real(kind=dp), intent(in) :: gval(nnz), y(ndata)
  real(kind=dp), intent(out) :: x(nparm)

  integer :: i, k, k1, k2, j

  x = 0.0_dp
  do i = 1, ndata
    k1 = row_ptr(i)
    k2 = row_ptr(i+1) - 1
    do k = k1, k2
      j = gindx(k)
      if (j >= 1 .and. j <= nparm) then
        x(j) = x(j) + gval(k) * y(i)
      end if
    end do
  end do

end subroutine sparse_ATx


subroutine cgnr_solver(ndata, nparm, nnz, gindx, gval, row_ptr, atb, damping, x, max_iter)
  implicit none
  integer, intent(in) :: ndata, nparm, nnz, max_iter
  integer, intent(in) :: gindx(nnz), row_ptr(ndata+1)
  real(kind=dp), intent(in) :: gval(nnz), atb(nparm), damping
  real(kind=dp), intent(inout) :: x(nparm)

  call lsqr_damped_solver(ndata, nparm, nnz, gindx, gval, &
                           row_ptr, atb, damping, x, max_iter)

end subroutine cgnr_solver


subroutine smooth_delta_m2d(delta_m, alpha)
  implicit none
  real(kind=dp), intent(inout) :: delta_m(n_inv_parms)
  real(kind=dp), intent(in) :: alpha

  integer :: ig, ix, iz, idx, jdx, vtype, nx_vg, nz_vg
  real(kind=dp), allocatable :: tmp(:)
  real(kind=dp) :: avg
  integer :: n_nb
  type(Tvelocity_grid2d), pointer :: vg

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    nx_vg = vg%nx
    nz_vg = vg%nz
    vtype = vgrids_to_be_inv(ig)%vtype

    allocate(tmp(nx_vg * nz_vg))
    tmp = 0.0_dp

    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0) then
          tmp((iz-1)*nx_vg + ix) = delta_m(idx)
        end if
      end do
    end do

    do iz = 2, nz_vg - 1
      do ix = 2, nx_vg - 1
        idx = vg%inv_idx(ix, iz)
        if (idx <= 0) cycle
        avg = 0.0_dp
        n_nb = 0
        jdx = vg%inv_idx(ix-1, iz)
        if (jdx > 0) then; avg = avg + delta_m(jdx); n_nb = n_nb + 1; end if
        jdx = vg%inv_idx(ix+1, iz)
        if (jdx > 0) then; avg = avg + delta_m(jdx); n_nb = n_nb + 1; end if
        jdx = vg%inv_idx(ix, iz-1)
        if (jdx > 0) then; avg = avg + delta_m(jdx); n_nb = n_nb + 1; end if
        jdx = vg%inv_idx(ix, iz+1)
        if (jdx > 0) then; avg = avg + delta_m(jdx); n_nb = n_nb + 1; end if
        if (n_nb > 0) then
          tmp((iz-1)*nx_vg + ix) = (1.0_dp - alpha) * delta_m(idx) + alpha * avg / dble(n_nb)
        end if
      end do
    end do

    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0) then
          delta_m(idx) = tmp((iz-1)*nx_vg + ix)
        end if
      end do
    end do

    deallocate(tmp)
  end do

end subroutine smooth_delta_m2d


subroutine update_model2d(delta_m, step_length, mvscale)
  implicit none
  real(kind=dp), intent(in) :: delta_m(n_inv_parms)
  real(kind=dp), intent(in) :: step_length
  real(kind=dp), intent(in) :: mvscale(n_inv_parms)

  integer :: ig, ii, iiface, ix, iz, idx
  integer :: nx_vg, nz_vg, vtype
  integer :: n_updated, n_skipped
  real(kind=dp) :: old_v, new_v, old_s, new_s, max_change
  type(Tvelocity_grid2d), pointer :: vg
  type(Tinterface2d), pointer :: ifc

  n_updated = 0
  n_skipped = 0
  max_change = 0.0_dp

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz

    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0) then
          old_v = vg%velocity(vtype, ix, iz)
          old_s = 1.0_dp / old_v
          new_s = old_s + step_length * mvscale(idx) * delta_m(idx)
          if (new_s > 0.0_dp) then
            new_v = 1.0_dp / new_s
          else
            new_v = old_v
          end if
          if (new_v < 2.0_dp) new_v = 2.0_dp
          if (new_v > 10.0_dp) new_v = 10.0_dp
          vg%velocity(vtype, ix, iz) = new_v
          if (abs(new_v - old_v) > max_change) max_change = abs(new_v - old_v)
          if (abs(new_v - old_v) > 1.0e-6_dp) n_updated = n_updated + 1
        end if
      end do
    end do
  end do

  write(*,'(A,I6,A,I6,A,F12.6)') '  update_model2d: updated=', n_updated, &
    ' skipped_water=', n_skipped, ' max_dv=', max_change

  ! Enforce velocity monotonicity near seafloor: first active node should not
  ! be slower than the node below it (avoids staircase artifact at convex seafloor)
  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    do ix = 1, nx_vg
      do iz = 1, nz_vg - 1
        if (vg%active(vtype, ix, iz)) then
          ! This is the first active node - check if it's slower than the one below
          if (vg%velocity(vtype, ix, iz) < vg%velocity(vtype, ix, iz+1)) then
            vg%velocity(vtype, ix, iz) = vg%velocity(vtype, ix, iz+1)
          end if
          exit
        end if
      end do
    end do
  end do

  do ii = 1, n_inv_iface
    iiface = ifaces_to_be_inv(ii)
    ifc => intrface(iiface)

    do ix = 1, ifc%nx
      idx = ifc%start_index + ix
      if (idx >= 1 .and. idx <= n_inv_parms) then
        ifc%z_nodes(ix) = ifc%z_nodes(ix) + step_length * mvscale(idx) * delta_m(idx)
      end if
    end do
  end do

  if (n_inv_iface > 0) then
    call rebuild_geometry2d()
  end if

end subroutine update_model2d


subroutine save_model2d(vel_save, z_save)
  implicit none
  real(kind=dp), allocatable, intent(out) :: vel_save(:)
  real(kind=dp), allocatable, intent(out) :: z_save(:)

  integer :: ig, ix, iz, idx, vtype, nx_vg, nz_vg, ii, iiface
  integer :: n_vel, n_z, ip
  type(Tvelocity_grid2d), pointer :: vg
  type(Tinterface2d), pointer :: ifc

  n_vel = 0
  n_z = 0
  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    n_vel = n_vel + vg%nx * vg%nz
  end do
  do ii = 1, n_inv_iface
    iiface = ifaces_to_be_inv(ii)
    ifc => intrface(iiface)
    n_z = n_z + ifc%nx
  end do

  allocate(vel_save(n_vel), z_save(n_z))

  ip = 0
  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx; nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        ip = ip + 1
        vel_save(ip) = vg%velocity(vtype, ix, iz)
      end do
    end do
  end do

  ip = 0
  do ii = 1, n_inv_iface
    iiface = ifaces_to_be_inv(ii)
    ifc => intrface(iiface)
    do ix = 1, ifc%nx
      ip = ip + 1
      z_save(ip) = ifc%z_nodes(ix)
    end do
  end do

end subroutine save_model2d


subroutine restore_model2d(vel_save, z_save)
  implicit none
  real(kind=dp), intent(in) :: vel_save(:), z_save(:)

  integer :: ig, ix, iz, idx, vtype, nx_vg, nz_vg, ii, iiface
  integer :: ip
  type(Tvelocity_grid2d), pointer :: vg
  type(Tinterface2d), pointer :: ifc

  ip = 0
  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx; nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        ip = ip + 1
        vg%velocity(vtype, ix, iz) = vel_save(ip)
      end do
    end do
  end do

  ip = 0
  do ii = 1, n_inv_iface
    iiface = ifaces_to_be_inv(ii)
    ifc => intrface(iiface)
    do ix = 1, ifc%nx
      ip = ip + 1
      ifc%z_nodes(ix) = z_save(ip)
    end do
  end do

  if (n_inv_iface > 0) then
    call rebuild_geometry2d()
  end if

end subroutine restore_model2d


subroutine compute_rms_after_forward(rms_out)
  implicit none
  real(kind=dp), intent(out) :: rms_out

  integer :: ircv, iray, ndata
  real(kind=dp) :: res_check, sum_sq
  integer :: t0, t1, rate
  real(kind=dp) :: t_fwd

  ! Line search only needs receiver travel times, not ray paths: skip the
  ! Frechet derivative pass to save one full sweep over all ray paths.
  skip_frechet_in_forward = .true.
  call system_clock(t0, rate)
  call run_forward2d()
  call system_clock(t1)
  skip_frechet_in_forward = .false.
  t_fwd = real(t1 - t0, kind=dp) / real(rate, kind=dp)
  write(*,'(A,F9.2,A)') '    forward (FMM+rays, no frechet): ', t_fwd, ' s'

  ndata = 0
  sum_sq = 0.0_dp
  do ircv = 1, n_receivers
    do iray = 1, receiver(ircv)%n_rays
      if (receiver(ircv)%ray(iray)%valid .and. &
          receiver(ircv)%ray(iray)%obs_time >= 0.0_dp) then
        res_check = receiver(ircv)%ray(iray)%obs_time - &
                    receiver(ircv)%ray(iray)%receiver_time
        if (abs(res_check) <= inv_crit_chi) then
          ndata = ndata + 1
          sum_sq = sum_sq + res_check**2
        end if
      end if
    end do
  end do

  if (ndata > 0) then
    rms_out = sqrt(sum_sq / ndata)
  else
    rms_out = huge_time
  end if

end subroutine compute_rms_after_forward


subroutine write_residuals2d(filename)
  implicit none
  character(len=*), intent(in) :: filename

  integer :: ircv, iray, iunit, ndata
  real(kind=dp) :: res_check, obs_t, calc_t, rcv_x, rcv_z
  integer :: src_id, phase_id

  iunit = 30
  open(unit=iunit, file=trim(filename), status='replace', action='write')
  write(iunit,'(A)') '# src_id  rcv_x  rcv_z  phase_id  obs_time  calc_time  residual'

  ndata = 0
  do ircv = 1, n_receivers
    do iray = 1, receiver(ircv)%n_rays
      if (.not. receiver(ircv)%ray(iray)%valid) cycle
      if (receiver(ircv)%ray(iray)%obs_time < 0.0_dp) cycle

      src_id = receiver(ircv)%ray(iray)%source_id
      obs_t = receiver(ircv)%ray(iray)%obs_time
      calc_t = receiver(ircv)%ray(iray)%receiver_time
      res_check = obs_t - calc_t
      phase_id = receiver(ircv)%ray(iray)%raypath_id
      rcv_x = receiver(ircv)%x
      rcv_z = receiver(ircv)%z

      write(iunit,'(I6,F12.3,F10.3,I4,F12.6,F12.6,ES14.6)') &
        src_id, rcv_x, rcv_z, phase_id, obs_t, calc_t, res_check
      ndata = ndata + 1
    end do
  end do

  close(iunit)
  write(*,'(A,A,A,I6)') '  Written residuals to ', trim(filename), ' ndata=', ndata

end subroutine write_residuals2d


subroutine write_ray_density2d(filename)
  implicit none
  character(len=*), intent(in) :: filename

  integer :: ircv, iray, isec, ip, ix, iz
  integer :: ivx, ivz, vtype, iunit
  real(kind=dp) :: px, pz, wx, wz
  real(kind=dp), allocatable :: ray_count(:,:)
  type(Tvelocity_grid2d), pointer :: vg
  type(Tray_section2d), pointer :: raysec

  vg => vgrid(vgrids_to_be_inv(1)%igrid, vgrids_to_be_inv(1)%vtype)
  vtype = vgrids_to_be_inv(1)%vtype

  allocate(ray_count(vg%nx, vg%nz))
  ray_count = 0.0_dp

  do ircv = 1, n_receivers
    do iray = 1, receiver(ircv)%n_rays
      if (.not. receiver(ircv)%ray(iray)%valid) cycle
      if (receiver(ircv)%ray(iray)%obs_time < 0.0_dp) cycle

      do isec = 1, receiver(ircv)%ray(iray)%nsections
        raysec => receiver(ircv)%ray(iray)%section(isec)
        do ip = 2, raysec%npoints
          px = 0.5_dp * (raysec%point(1, ip) + raysec%point(1, ip-1))
          pz = 0.5_dp * (raysec%point(2, ip) + raysec%point(2, ip-1))

          ivx = floor((px - vg%x(1)) / vg%dx) + 1
          if (ivx < 1) ivx = 1
          if (ivx >= vg%nx) ivx = vg%nx - 1
          wx = (px - vg%x(ivx)) / vg%dx
          if (wx < 0.0_dp) wx = 0.0_dp
          if (wx > 1.0_dp) wx = 1.0_dp

          if (pz < vg%z(1) .or. pz >= vg%z(vg%nz)) cycle
          ivz = floor((pz - vg%z(1)) / vg%dz) + 1
          if (ivz < 1) ivz = 1
          if (ivz >= vg%nz) ivz = vg%nz - 1
          wz = (pz - vg%z(ivz)) / vg%dz
          if (wz < 0.0_dp) wz = 0.0_dp
          if (wz > 1.0_dp) wz = 1.0_dp

          ray_count(ivx, ivz) = ray_count(ivx, ivz) + (1.0_dp-wx)*(1.0_dp-wz)
          ray_count(ivx+1, ivz) = ray_count(ivx+1, ivz) + wx*(1.0_dp-wz)
          ray_count(ivx, ivz+1) = ray_count(ivx, ivz+1) + (1.0_dp-wx)*wz
          ray_count(ivx+1, ivz+1) = ray_count(ivx+1, ivz+1) + wx*wz
        end do
      end do
    end do
  end do

  iunit = 31
  open(unit=iunit, file=trim(filename), status='replace', action='write')
  write(iunit,'(I6,I6)') vg%nx, vg%nz
  ! Write x coordinates
  do ix = 1, vg%nx
    write(iunit,'(ES14.6)') vg%x(ix)
  end do
  ! Write z coordinates
  do iz = 1, vg%nz
    write(iunit,'(ES14.6)') vg%z(iz)
  end do
  ! Write ray count (sparse: only non-zero)
  do ix = 1, vg%nx
    do iz = 1, vg%nz
      if (ray_count(ix, iz) > 0.0_dp) then
        write(iunit,'(I6,I6,ES14.6)') ix, iz, ray_count(ix, iz)
      end if
    end do
  end do
  close(iunit)

  write(*,'(A,A,A,F10.1)') '  Written ray density to ', trim(filename), &
    ' max_count=', maxval(ray_count)

  deallocate(ray_count)

end subroutine write_ray_density2d


subroutine update_model_creeping2d(dmodel, mvscale)
  implicit none
  real(kind=dp), intent(in) :: dmodel(n_inv_parms)
  real(kind=dp), intent(in) :: mvscale(n_inv_parms)

  integer :: ig, ii, iiface, ix, iz, idx
  integer :: nx_vg, nz_vg, vtype
  integer :: n_updated
  real(kind=dp) :: cur_v, cur_s, new_v, new_s, ds, max_change
  type(Tvelocity_grid2d), pointer :: vg
  type(Tinterface2d), pointer :: ifc

  n_updated = 0
  max_change = 0.0_dp

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz

    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0) then
          cur_v = vg%velocity(vtype, ix, iz)
          cur_s = 1.0_dp / cur_v
          ds = mvscale(idx) * dmodel(idx)
          new_s = cur_s + ds
          if (new_s > 0.0_dp) then
            new_v = 1.0_dp / new_s
          else
            new_v = cur_v
          end if
          if (new_v < 2.0_dp) new_v = 2.0_dp
          if (new_v > 10.0_dp) new_v = 10.0_dp
          vg%velocity(vtype, ix, iz) = new_v
          if (abs(new_v - cur_v) > max_change) max_change = abs(new_v - cur_v)
          if (abs(new_v - cur_v) > 1.0e-6_dp) n_updated = n_updated + 1
        end if
      end do
    end do
  end do

  write(*,'(A,I6,A,F12.6)') '  update_model: updated=', n_updated, &
    ' max_dv=', max_change

  do ii = 1, n_inv_iface
    iiface = ifaces_to_be_inv(ii)
    ifc => intrface(iiface)

    do ix = 1, ifc%nx
      idx = ifc%start_index + ix
      if (idx >= 1 .and. idx <= n_inv_parms) then
        ifc%z_nodes(ix) = ifc%z_nodes(ix) + mvscale(idx) * dmodel(idx)
      end if
    end do
  end do

  if (n_inv_iface > 0) then
    call rebuild_geometry2d()
  end if

end subroutine update_model_creeping2d


subroutine update_model_from_init2d(dmodel_total, dmodel_total_sum, mvscale)
  implicit none
  real(kind=dp), intent(in) :: dmodel_total(n_inv_parms)
  real(kind=dp), intent(inout) :: dmodel_total_sum(n_inv_parms)
  real(kind=dp), intent(in) :: mvscale(n_inv_parms)

  integer :: ig, ii, iiface, ix, iz, idx
  integer :: nx_vg, nz_vg, vtype
  integer :: n_updated
  real(kind=dp) :: cur_v, init_v, init_s, new_v, new_s, max_change
  type(Tvelocity_grid2d), pointer :: vg
  type(Tinterface2d), pointer :: ifc

  n_updated = 0
  max_change = 0.0_dp

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz

    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0) then
          cur_v = vg%velocity(vtype, ix, iz)
          init_v = vg%velocity_init(vtype, ix, iz)
          init_s = 1.0_dp / init_v
          new_s = init_s + mvscale(idx) * dmodel_total_sum(idx)
          if (new_s > 0.0_dp) then
            new_v = 1.0_dp / new_s
          else
            new_v = init_v
          end if
          if (new_v < 2.0_dp) new_v = 2.0_dp
          if (new_v > 10.0_dp) new_v = 10.0_dp
          vg%velocity(vtype, ix, iz) = new_v
          if (abs(new_v - cur_v) > 1.0e-6_dp) then
            dmodel_total_sum(idx) = (1.0_dp/new_v - init_s) / mvscale(idx)
          end if
          if (abs(new_v - cur_v) > max_change) max_change = abs(new_v - cur_v)
          if (abs(new_v - cur_v) > 1.0e-6_dp) n_updated = n_updated + 1
        end if
      end do
    end do
  end do

  write(*,'(A,I6,A,F12.6)') '  update_model: updated=', n_updated, &
    ' max_dv=', max_change

  do ii = 1, n_inv_iface
    iiface = ifaces_to_be_inv(ii)
    ifc => intrface(iiface)

    do ix = 1, ifc%nx
      idx = ifc%start_index + ix
      if (idx >= 1 .and. idx <= n_inv_parms) then
        ifc%z_nodes(ix) = ifc%z_nodes_init(ix) + mvscale(idx) * dmodel_total_sum(idx)
      end if
    end do
  end do

  if (n_inv_iface > 0) then
    call rebuild_geometry2d()
  end if

end subroutine update_model_from_init2d

subroutine update_model_from_increment2d(dmodel_total, mvscale)
  implicit none
  real(kind=dp), intent(in) :: dmodel_total(n_inv_parms)
  real(kind=dp), intent(in) :: mvscale(n_inv_parms)

  integer :: ig, ii, iiface, ix, iz, idx
  integer :: nx_vg, nz_vg, vtype
  integer :: n_updated
  real(kind=dp) :: cur_v, cur_s, new_v, new_s, max_change
  type(Tvelocity_grid2d), pointer :: vg
  type(Tinterface2d), pointer :: ifc

  n_updated = 0
  max_change = 0.0_dp

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz

    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0) then
          cur_v = vg%velocity(vtype, ix, iz)
          cur_s = 1.0_dp / cur_v
          new_s = cur_s + mvscale(idx) * dmodel_total(idx)
          if (new_s > 0.0_dp) then
            new_v = 1.0_dp / new_s
          else
            new_v = cur_v
          end if
          if (new_v < 2.0_dp) new_v = 2.0_dp
          if (new_v > 10.0_dp) new_v = 10.0_dp
          vg%velocity(vtype, ix, iz) = new_v
          if (abs(new_v - cur_v) > max_change) max_change = abs(new_v - cur_v)
          if (abs(new_v - cur_v) > 1.0e-6_dp) n_updated = n_updated + 1
        end if
      end do
    end do
  end do

  write(*,'(A,I6,A,F12.6)') '  update_model_increment: updated=', n_updated, &
    ' max_dv=', max_change

  do ii = 1, n_inv_iface
    iiface = ifaces_to_be_inv(ii)
    ifc => intrface(iiface)

    do ix = 1, ifc%nx
      idx = ifc%start_index + ix
      if (idx >= 1 .and. idx <= n_inv_parms) then
        ifc%z_nodes(ix) = ifc%z_nodes(ix) + mvscale(idx) * dmodel_total(idx)
      end if
    end do
  end do

  if (n_inv_iface > 0) then
    call rebuild_geometry2d()
  end if

end subroutine update_model_from_increment2d


subroutine rebuild_geometry2d()
  implicit none

  integer :: n, ireg, ix, iz

  do n = 1, n_intersections
    if (associated(intersection(n)%x)) then
      deallocate(intersection(n)%x)
      intersection(n)%x => null()
    end if
    if (associated(intersection(n)%z)) then
      deallocate(intersection(n)%z)
      intersection(n)%z => null()
    end if
    if (associated(intersection(n)%normal)) then
      deallocate(intersection(n)%normal)
      intersection(n)%normal => null()
    end if
    if (associated(intersection(n)%ccells)) then
      deallocate(intersection(n)%ccells)
      intersection(n)%ccells => null()
    end if
    if (associated(intersection(n)%n_inodes)) then
      deallocate(intersection(n)%n_inodes)
      intersection(n)%n_inodes => null()
    end if
    if (associated(intersection(n)%inodes)) then
      deallocate(intersection(n)%inodes)
      intersection(n)%inodes => null()
    end if
    if (associated(intersection(n)%ccell_from_inode)) then
      deallocate(intersection(n)%ccell_from_inode)
      intersection(n)%ccell_from_inode => null()
    end if
    if (associated(intersection(n)%rright_node_id)) then
      deallocate(intersection(n)%rright_node_id)
      intersection(n)%rright_node_id => null()
    end if
    if (associated(intersection(n)%rleft_node_id)) then
      deallocate(intersection(n)%rleft_node_id)
      intersection(n)%rleft_node_id => null()
    end if
    if (associated(intersection(n)%arrivaltime)) then
      deallocate(intersection(n)%arrivaltime)
      intersection(n)%arrivaltime => null()
    end if
    if (associated(intersection(n)%time_gradient)) then
      deallocate(intersection(n)%time_gradient)
      intersection(n)%time_gradient => null()
    end if
    if (associated(intersection(n)%vel_top)) then
      deallocate(intersection(n)%vel_top)
      intersection(n)%vel_top => null()
    end if
    if (associated(intersection(n)%vel_bot)) then
      deallocate(intersection(n)%vel_bot)
      intersection(n)%vel_bot => null()
    end if
  end do

  do ireg = 1, n_regions
    if (associated(region(ireg)%node)) then
      deallocate(region(ireg)%node)
      region(ireg)%node => null()
    end if
    if (associated(region(ireg)%node_status)) then
      deallocate(region(ireg)%node_status)
      region(ireg)%node_status => null()
    end if
    if (associated(region(ireg)%arrivaltime)) then
      deallocate(region(ireg)%arrivaltime)
      region(ireg)%arrivaltime => null()
    end if
    if (associated(region(ireg)%time_gradient)) then
      deallocate(region(ireg)%time_gradient)
      region(ireg)%time_gradient => null()
    end if
    if (associated(region(ireg)%velocity)) then
      deallocate(region(ireg)%velocity)
      region(ireg)%velocity => null()
    end if
    if (associated(region(ireg)%x)) then
      deallocate(region(ireg)%x)
      region(ireg)%x => null()
    end if
    if (associated(region(ireg)%z)) then
      deallocate(region(ireg)%z)
      region(ireg)%z => null()
    end if
    if (associated(region(ireg)%region_rnode_id)) then
      deallocate(region(ireg)%region_rnode_id)
      region(ireg)%region_rnode_id => null()
    end if
    if (associated(region(ireg)%itop_rnode_id)) then
      deallocate(region(ireg)%itop_rnode_id)
      region(ireg)%itop_rnode_id => null()
    end if
    if (associated(region(ireg)%ibot_rnode_id)) then
      deallocate(region(ireg)%ibot_rnode_id)
      region(ireg)%ibot_rnode_id => null()
    end if
    region(ireg)%nnode = 0
  end do

  if (associated(pgrid%rnode_id)) then
    pgrid%rnode_id = 0
  end if
  if (associated(pgrid%node_region)) then
    pgrid%node_region = 0
  end if

  ! Cut-cell lookup must be cleared before rebuilding the intersections.
  ! find_intersection2d only fills the cells cut by the *current* interface
  ! position and otherwise reuses whatever pointer is already there, so a cell
  ! that was cut by an earlier (different) interface geometry keeps a stale
  ! cell index. That index refers to a reallocated ccells array and corrupts
  ! reflected-phase tracing. Releasing the pointers forces a fresh, zeroed
  ! allocation on the next find_intersection2d call.
  if (associated(pgrid%ccind_from_2dc)) then
    do ix = 1, pgrid%nx
      do iz = 1, pgrid%nz
        if (associated(pgrid%ccind_from_2dc(ix,iz)%p)) then
          deallocate(pgrid%ccind_from_2dc(ix,iz)%p)
          nullify(pgrid%ccind_from_2dc(ix,iz)%p)
        end if
      end do
    end do
  end if

  do n = 1, n_intersections
    call find_intersection2d(intersection(n), intrface(n), pgrid)
  end do

  do ireg = 1, n_regions
    call define_region2d(region(ireg), region(ireg)%itop, region(ireg)%ibot, pgrid)
  end do

  call tag_regular_nodes2d(pgrid)
  call velocities_to_grid2d(pgrid)

  do n = 1, n_intersections
    call velocities_to_intersection2d(intersection(n))
  end do
  do ireg = 1, n_regions
    if (region(ireg)%nnode > 0) then
      call velocities_to_region2d(region(ireg), pgrid)
    end if
  end do

  do ireg = 2, n_regions
    if (region(ireg)%nnode > 0) then
      call tag_active_vgrid_nodes2d(region(ireg))
    end if
  end do

end subroutine rebuild_geometry2d


function compute_rms2d(ndata, residual) result(rms)
  implicit none
  integer, intent(in) :: ndata
  real(kind=dp), intent(in) :: residual(ndata)
  real(kind=dp) :: rms

  integer :: i
  real(kind=dp) :: sumsq

  sumsq = 0.0_dp
  do i = 1, ndata
    sumsq = sumsq + residual(i) * residual(i)
  end do

  if (ndata > 0) then
    rms = sqrt(sumsq / real(ndata, kind=dp))
  else
    rms = 0.0_dp
  end if

end function compute_rms2d


subroutine write_model2d(iter)
  implicit none
  integer, intent(in) :: iter

  integer :: ig, iiface, ii, ix, iz, iunit
  integer :: nx_vg, nz_vg, vtype
  integer :: iz_out, iz_lo, iz_hi
  real(kind=dp) :: vwater, vair, wz, z_abs, v_interp
  character(len=256) :: filename
  type(Tvelocity_grid2d), pointer :: vg
  type(Tinterface2d), pointer :: ifc

  iunit = 20

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz

    write(filename,'(A,I0,A,I0,A)') 'vgrid_', vgrids_to_be_inv(ig)%igrid, &
      '_vtype', vtype, '.smesh'

    open(unit=iunit, file=trim(filename), status='replace', action='write')

    ! Water and air velocities are taken from the input smesh header.
    ! Hardcoding them made the written model inconsistent with the model the
    ! inversion actually used: any later forward run on this file then has a
    ! systematic ~40 ms error in the water leg.
    vwater = vwater_global
    vair = vair_global

    if (allocated(z_orig) .and. nz_orig > 0) then
      ! Output in original format: z relative to seafloor, no water layer
      write(iunit,'(I6,I6,F10.3,F10.3)') nx_vg, nz_orig, vwater, vair
      write(iunit,'(10F12.4)') (vg%x(ix), ix=1,nx_vg)
      write(iunit,'(10F12.4)') (intrface(2)%z_nodes(ix), ix=1,nx_vg)
      write(iunit,'(10F12.4)') (z_orig(iz), iz=1,nz_orig)
      do ix = 1, nx_vg
        do iz = 1, nz_orig
          z_abs = intrface(2)%z_nodes(ix) + z_orig(iz)
          ! Find bracketing indices in resampled grid
          iz_lo = 1
          iz_hi = nz_vg
          do iz_out = 1, nz_vg - 1
            if (vg%z(iz_out) <= z_abs .and. vg%z(iz_out+1) >= z_abs) then
              iz_lo = iz_out
              iz_hi = iz_out + 1
              exit
            end if
          end do
          if (z_abs <= vg%z(1)) then
            v_interp = vg%velocity(vtype, ix, 1)
          else if (z_abs >= vg%z(nz_vg)) then
            v_interp = vg%velocity(vtype, ix, nz_vg)
          else
            wz = (z_abs - vg%z(iz_lo)) / (vg%z(iz_hi) - vg%z(iz_lo))
            ! Avoid interpolating across water-crust boundary:
            ! if iz_lo is in water (inactive or velocity < 2.0), use iz_hi velocity
            if (vg%active(vtype, ix, iz_lo) .and. &
                vg%velocity(vtype, ix, iz_lo) >= 2.0_dp) then
              v_interp = (1.0_dp - wz) * vg%velocity(vtype, ix, iz_lo) + &
                          wz * vg%velocity(vtype, ix, iz_hi)
            else
              v_interp = vg%velocity(vtype, ix, iz_hi)
            end if
          end if
          write(iunit,'(F10.4)', advance='no') v_interp
        end do
        write(iunit,*)
      end do
    else
      ! Fallback: output resampled grid as-is
      write(iunit,'(I6,I6,F10.3,F10.3)') nx_vg, nz_vg, vwater, vair
      write(iunit,'(10F12.4)') (vg%x(ix), ix=1,nx_vg)
      write(iunit,'(10F12.4)') (intrface(2)%z_nodes(ix), ix=1,nx_vg)
      write(iunit,'(10F12.4)') (vg%z(iz), iz=1,nz_vg)
      do ix = 1, nx_vg
        write(iunit,'(10F10.4)') (vg%velocity(vtype,ix,iz), iz=1,nz_vg)
      end do
    end if

    close(iunit)
  end do

  do ii = 1, n_inv_iface
    iiface = ifaces_to_be_inv(ii)
    ifc => intrface(iiface)

    write(filename,'(A,I0,A)') 'iface_', iiface, '.refl'

    open(unit=iunit, file=trim(filename), status='replace', action='write')
    do ix = 1, ifc%nx
      write(iunit,'(F12.3,F12.3)') ifc%x(ix), ifc%z_nodes(ix)
    end do
    close(iunit)
  end do

  write(filename,'(A,I0,A)') 'model_iter', iter, '.dat'
  open(unit=iunit, file=trim(filename), status='replace', action='write')
  write(iunit,'(A,I0)') 'Iteration ', iter
  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    write(iunit,'(A,I0,A,I0)') 'vgrid ', vgrids_to_be_inv(ig)%igrid, &
      ' vtype ', vtype
    do ix = 1, nx_vg
      write(iunit,'(10F10.4)') (vg%velocity(vtype,ix,iz), iz=1,nz_vg)
    end do
  end do
  do ii = 1, n_inv_iface
    iiface = ifaces_to_be_inv(ii)
    ifc => intrface(iiface)
    write(iunit,'(A,I0)') 'iface ', iiface
    write(iunit,'(10F12.4)') (ifc%z_nodes(ix), ix=1,ifc%nx)
  end do
  close(iunit)

end subroutine write_model2d


subroutine read_frechet_data2d(ndata_total, nnz_total, gindx, gval, residual, row_ptr)
  implicit none
  integer, intent(out) :: ndata_total, nnz_total
  integer, intent(out) :: gindx(:), row_ptr(:)
  real(kind=dp), intent(out) :: gval(:), residual(:)

  integer :: ircv, iray, idata, k, kk
  integer :: ipath, max_path
  integer, allocatable :: n_per_path(:)
  real(kind=dp), allocatable :: ssq_per_path(:)
  real(kind=dp) :: res_check

  max_path = 0
  do ircv = 1, n_receivers
    do iray = 1, receiver(ircv)%n_rays
      if (receiver(ircv)%ray(iray)%raypath_id > max_path) then
        max_path = receiver(ircv)%ray(iray)%raypath_id
      end if
    end do
  end do
  if (max_path < 3) max_path = 3

  allocate(n_per_path(max_path), ssq_per_path(max_path))
  n_per_path = 0
  ssq_per_path = 0.0_dp

  ndata_total = 0
  nnz_total = 0
  idata = 0
  kk = 0
  row_ptr(1) = 1

  do ircv = 1, n_receivers
    do iray = 1, receiver(ircv)%n_rays
      if (.not. receiver(ircv)%ray(iray)%valid) cycle
      if (receiver(ircv)%ray(iray)%obs_time < 0.0_dp) cycle

      res_check = receiver(ircv)%ray(iray)%obs_time - receiver(ircv)%ray(iray)%receiver_time
      if (abs(res_check) > inv_crit_chi) then
        if (abs(res_check) > 0.2_dp) then
          write(*,'(A,I0,A,I0,A,I0,A,F12.6,A,F12.6,A,F12.6)') '    OUTLIER: rec=', ircv, &
            ' ray=', iray, ' path=', receiver(ircv)%ray(iray)%raypath_id, &
            ' obs=', receiver(ircv)%ray(iray)%obs_time, &
            ' calc=', receiver(ircv)%ray(iray)%receiver_time, &
            ' res=', res_check
        end if
        cycle
      end if

      idata = idata + 1
      ndata_total = ndata_total + 1

      residual(idata) = res_check

      ipath = receiver(ircv)%ray(iray)%raypath_id
      if (ipath >= 1 .and. ipath <= max_path) then
        n_per_path(ipath) = n_per_path(ipath) + 1
        ssq_per_path(ipath) = ssq_per_path(ipath) + residual(idata)**2
      end if

      do k = 1, receiver(ircv)%ray(iray)%n_pdev
        kk = kk + 1
        nnz_total = nnz_total + 1
        if (kk <= size(gindx)) then
          gindx(kk) = receiver(ircv)%ray(iray)%pdev_indx(k)
          gval(kk) = receiver(ircv)%ray(iray)%pdev(k)
        end if
      end do

      if (idata+1 <= size(row_ptr)) then
        row_ptr(idata+1) = kk + 1
      end if

    end do
  end do

  write(*,'(A)') '  RMS by phase:'
  do ipath = 1, max_path
    if (n_per_path(ipath) > 0) then
      write(*,'(A,I1,A,I6,A,F12.6)') '    path=', ipath, &
        ' n=', n_per_path(ipath), &
        ' RMS=', sqrt(ssq_per_path(ipath) / n_per_path(ipath))
    end if
  end do

  deallocate(n_per_path, ssq_per_path)

end subroutine read_frechet_data2d


subroutine check_vg_pgrid_consistency()
  implicit none

  integer :: ig, vtype, ix, iz, ivx, ivz, nx_vg, nz_vg
  real(kind=dp) :: px, pz, pz_rel, z_sf, wx_z, wx, wz
  real(kind=dp) :: v_vg, v_pgrid, max_diff, rel_diff
  integer :: ivx_z, n_compare, max_diff_ix, max_diff_iz
  type(Tvelocity_grid2d), pointer :: vg
  type(Tpropagation_grid2d), pointer :: pgrid

  ig = 1
  vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
  vtype = vgrids_to_be_inv(ig)%vtype
  nx_vg = vg%nx
  nz_vg = vg%nz

  pgrid => region(2)%grid

  max_diff = 0.0_dp
  n_compare = 0
  max_diff_ix = 0
  max_diff_iz = 0

  do ix = 1, pgrid%nx
    px = pgrid%x(ix)
    do iz = 1, pgrid%nz
      pz = pgrid%z(iz)

      ivx_z = floor((px - intrface(2)%x(1)) / intrface(2)%dx) + 1
      if (ivx_z < 1) ivx_z = 1
      if (ivx_z >= intrface(2)%nx) ivx_z = intrface(2)%nx - 1
      wx_z = (px - intrface(2)%x(ivx_z)) / intrface(2)%dx
      if (wx_z < 0.0_dp) wx_z = 0.0_dp
      if (wx_z > 1.0_dp) wx_z = 1.0_dp
      z_sf = (1.0_dp - wx_z) * intrface(2)%z_nodes(ivx_z) + &
             wx_z * intrface(2)%z_nodes(ivx_z + 1)
      pz_rel = pz - z_sf

      if (pz_rel < vg%z(1) .or. pz_rel > vg%z(nz_vg)) cycle

      v_pgrid = pgrid%velocity(vtype, ix, iz)

      ivx = floor((px - vg%x(1)) / vg%dx) + 1
      if (ivx < 1) ivx = 1
      if (ivx >= nx_vg) ivx = nx_vg - 1
      if (pz_rel >= vg%z(nz_vg)) then
        ivz = nz_vg - 1
        wz = 1.0_dp
      else
        ivz = floor((pz_rel - vg%z(1)) / vg%dz) + 1
        if (ivz < 1) ivz = 1
        if (ivz >= nz_vg) ivz = nz_vg - 1
        wz = (pz_rel - vg%z(ivz)) / vg%dz
      end if
      wx = (px - vg%x(ivx)) / vg%dx
      if (wx < 0.0_dp) wx = 0.0_dp
      if (wx > 1.0_dp) wx = 1.0_dp
      if (wz < 0.0_dp) wz = 0.0_dp
      if (wz > 1.0_dp) wz = 1.0_dp

      v_vg = (1.0_dp-wx)*(1.0_dp-wz)*vg%velocity(vtype,ivx,ivz) + &
             wx*(1.0_dp-wz)*vg%velocity(vtype,ivx+1,ivz) + &
             (1.0_dp-wx)*wz*vg%velocity(vtype,ivx,ivz+1) + &
             wx*wz*vg%velocity(vtype,ivx+1,ivz+1)

      n_compare = n_compare + 1
      if (abs(v_pgrid - v_vg) > max_diff) then
        max_diff = abs(v_pgrid - v_vg)
        rel_diff = max_diff / v_vg * 100.0_dp
        max_diff_ix = ix
        max_diff_iz = iz
      end if
    end do
  end do

  write(*,'(A,I6,A,I6)') '    vg-pgrid comparison: n_compare=', n_compare, &
    ' nx_vg=', nx_vg, ' nz_vg=', nz_vg
  write(*,'(A,ES12.4,A,F8.3,A,I6,A,I6)') '    max |v_pgrid - v_vg|=', max_diff, &
    ' rel=', rel_diff, '% at pgrid(', max_diff_ix, ',', max_diff_iz, ')'

  if (max_diff_ix > 0 .and. max_diff_iz > 0) then
    px = pgrid%x(max_diff_ix)
    pz = pgrid%z(max_diff_iz)
    write(*,'(A,F10.3,A,F10.3,A,F10.3,A,F10.3)') '    px=', px, ' pz=', pz, &
      ' v_pgrid=', pgrid%velocity(vtype,max_diff_ix,max_diff_iz), ' v_vg=', v_vg
  end if

end subroutine check_vg_pgrid_consistency


subroutine run_forward2d()
  implicit none

  integer :: isrc, ircv, iray, ipath, itf_idx, itf, ireg
  integer :: ntf, vtype, n, m, mm, nn, n_updated
  integer :: i1, i2, i3, icell, ixx, izz, ii, jj, ix, m_top, m_below
  real(kind=dp) :: dtdx, dtdz, dtdz_sq, dtdz_orig, v_loc, dx_if
  real(kind=dp), allocatable :: t_hw(:), v_hw(:)
  logical :: found
  type(Tsource2d), pointer :: src
  type(Tregion2d), pointer :: reg
  type(Ttime_field2d), pointer :: tf, tf_prev
  type(Tintersection2d), pointer :: isec

  integer :: t_c0, t_c1, t_c2, t_c3, t_rate
  real(kind=dp) :: t_fmm_acc, t_trace_acc, t_fr_acc

  t_fmm_acc = 0.0_dp
  t_trace_acc = 0.0_dp
  t_fr_acc = 0.0_dp
  call system_clock(t_c0, t_rate)

  call velocities_to_grid2d(pgrid)

  write(*,'(A,F10.4,A,F10.4)') '  pgrid v range:', &
    minval(pgrid%velocity(1, pgrid%nx/4+1:3*pgrid%nx/4, :)), &
    ' to', maxval(pgrid%velocity(1,:,:))

  do itf = 1, n_intersections
    call velocities_to_intersection2d(intersection(itf))
  end do

  do ireg = 1, n_regions
    call velocities_to_region2d(region(ireg), pgrid)
    write(*,'(A,I1,A,F10.4,A,F10.4)') '  region ', ireg, ' v range:', &
      minval(region(ireg)%velocity(VTYPE_P, 1:region(ireg)%nnode)), &
      ' to', maxval(region(ireg)%velocity(VTYPE_P, 1:region(ireg)%nnode))
  end do

  do isrc = 1, n_sources
    src => source(isrc)

    do itf = 1, src%n_time_fields
      tf => src%time_field(itf)
      if (.not. associated(tf%reg)) cycle
      reg => tf%reg
      if (associated(tf%arrivaltime)) then
        if (size(tf%arrivaltime) /= reg%nnode) then
          deallocate(tf%arrivaltime)
          tf%arrivaltime => null()
        end if
      end if
      if (.not. associated(tf%arrivaltime)) then
        allocate(tf%arrivaltime(reg%nnode))
        tf%arrivaltime = huge_time
      end if
      if (associated(tf%time_gradient)) then
        if (size(tf%time_gradient, 2) /= reg%nnode) then
          deallocate(tf%time_gradient)
          tf%time_gradient => null()
        end if
      end if
      if (.not. associated(tf%time_gradient)) then
        allocate(tf%time_gradient(2, reg%nnode))
        tf%time_gradient = 0.0_dp
      end if
    end do

    call system_clock(t_c0)
    do ipath = 1, src%n_paths
      if (.not. src%path(ipath)%valid) cycle
      if (.not. src%path(ipath)%used) cycle

      ntf = src%path(ipath)%n_tf

      do itf_idx = 1, ntf
        itf = src%path(ipath)%tf_sequence(itf_idx)
        if (itf < 1 .or. itf > src%n_time_fields) cycle

        tf => src%time_field(itf)
        reg => tf%reg
        vtype = tf%vtype

        reg%node_status = NODE_FAR
        reg%arrivaltime = huge_time
        reg%time_gradient = 0.0_dp

        if (itf_idx == 1) then
          call initialize_region_from_source2d(reg, src%x, src%z, 0.0_dp, vtype)
        else
          if (associated(tf%istart)) then
            tf_prev => src%time_field(src%path(ipath)%tf_sequence(itf_idx-1))
            if (associated(tf_prev%reg)) then
              tf_prev%reg%arrivaltime = tf_prev%arrivaltime
              tf_prev%reg%time_gradient = tf_prev%time_gradient
            end if
            call transfer_times_to_intersection(tf%istart, tf_prev)
            if ((ipath == 2 .and. itf_idx == 2) .or. &
                  (ipath == 3 .and. itf_idx == 2) .or. &
                  (ipath == 4 .and. itf_idx == 2) .or. &
                  (ipath == 4 .and. itf_idx == 3) .or. &
                  (ipath == 4 .and. itf_idx == 4) .or. &
                  (ipath == 5 .and. itf_idx == 2) .or. &
                  (ipath == 5 .and. itf_idx == 3) .or. &
                  (ipath == 5 .and. itf_idx == 4)) then
              call reflect_gradient_at_interface(tf%istart, reg, vtype)
            end if
            reg%node_status = NODE_FAR
            reg%arrivaltime = huge_time
            reg%time_gradient = 0.0_dp
            call sweep_region_from_interface2d(reg, tf%istart, vtype, src)
          end if
        end if

        call propagate2d(reg, vtype)

        ! Fix boundary gradient for Pg_ss/Pg_ss2 water-layer time fields
        ! Pg_ss:  tf6(itf_idx=2,up), tf7(itf_idx=3,down), tf8(itf_idx=4,up)
        ! Pg_ss2: tf9(itf_idx=1,up), tf10(itf_idx=2,down), tf12(itf_idx=4,up)
        if (reg%id == 1) then
          if ((ipath == 4 .and. (itf_idx == 2 .or. itf_idx == 4)) .or. &
              (ipath == 5 .and. (itf_idx == 1 .or. itf_idx == 4))) then
            ! Upward FMM: dtdz must be < 0 at top boundary
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
          else if ((ipath == 4 .and. itf_idx == 3) .or. &
                   (ipath == 5 .and. itf_idx == 2)) then
            ! Downward FMM: dtdz must be > 0 at top boundary
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
          end if
        end if

        if (ipath == 3 .and. itf_idx == 2) then
          isec => tf%istart
          if (isec%nnode > 0) then
            allocate(t_hw(isec%nnode), v_hw(isec%nnode))
            do n = 1, isec%nnode
              m = 0
              if (associated(isec%rright_node_id)) m = isec%rright_node_id(n)
              if (m == 0 .and. associated(isec%rleft_node_id)) m = isec%rleft_node_id(n)
              if (m > 0 .and. m <= reg%nnode) then
                t_hw(n) = reg%arrivaltime(m)
              else
                t_hw(n) = huge_time
              endif
              v_hw(n) = 0.0_dp
              if (m > 0 .and. m <= reg%nnode) then
                if (reg%node(m)%i1 == 0) then
                  i1 = 0
                  i2 = reg%node(m)%i2
                  i3 = reg%node(m)%i3
                  do nn = 1, 4
                    icell = isec%ccell_from_inode(nn, i3)
                    if (icell > 0) then
                      ixx = isec%ccells(icell)%ix
                      izz = isec%ccells(icell)%iz
                      do ii = 0, 1
                        do jj = 0, 1
                          ixx = isec%ccells(icell)%ix + ii
                          izz = isec%ccells(icell)%iz + jj
                          if (ixx >= 1 .and. ixx <= pgrid%nx .and. izz >= 1 .and. izz <= pgrid%nz) then
                            if (pgrid%node_region(ixx, izz) == reg%id) then
                              mm = pgrid%rnode_id(ixx, izz)
                              if (mm > 0 .and. mm <= reg%nnode) then
                                if (reg%node(mm)%i1 > 0) then
                                  if (reg%z(mm) > isec%z(n) + 0.5_dp * pgrid%dz) then
                                    v_hw(n) = max(v_hw(n), reg%velocity(VTYPE_P, mm))
                                  endif
                                endif
                              endif
                            endif
                          endif
                        end do
                      end do
                    endif
                  end do
                else
                  v_hw(n) = reg%velocity(VTYPE_P, m)
                endif
              endif
            end do

            dx_if = isec%x(2) - isec%x(1)
            do n = 2, isec%nnode
              v_loc = v_hw(n)
              if (v_loc > 0.0_dp) then
                t_hw(n) = min(t_hw(n), t_hw(n-1) + dx_if / v_loc)
              endif
            end do
            do n = isec%nnode-1, 1, -1
              v_loc = v_hw(n)
              if (v_loc > 0.0_dp) then
                t_hw(n) = min(t_hw(n), t_hw(n+1) + dx_if / v_loc)
              endif
            end do

            n_updated = 0
            do n = 1, isec%nnode
              m = 0
              if (associated(isec%rright_node_id)) m = isec%rright_node_id(n)
              if (m == 0 .and. associated(isec%rleft_node_id)) m = isec%rleft_node_id(n)
              if (m > 0 .and. m <= reg%nnode) then
                if (t_hw(n) < reg%arrivaltime(m) - 1.0e-10_dp) then
                  reg%arrivaltime(m) = t_hw(n)
                  n_updated = n_updated + 1
                endif
              endif
            end do

            do n = 1, isec%nnode
              m = 0
              if (associated(isec%rright_node_id)) m = isec%rright_node_id(n)
              if (m == 0 .and. associated(isec%rleft_node_id)) m = isec%rleft_node_id(n)
              if (m > 0 .and. m <= reg%nnode) then
                if (reg%node_status(m) == 0) then
                  if (reg%node(m)%i1 > 0) then
                    call fit_gradient_at_node2d(reg, m, reg%time_gradient(1,m), reg%time_gradient(2,m))
                  else
                    if (n > 1 .and. n < isec%nnode) then
                      dtdx = (t_hw(n+1) - t_hw(n-1)) / (2.0_dp * dx_if)
                    elseif (n == 1 .and. isec%nnode > 1) then
                      dtdx = (t_hw(2) - t_hw(1)) / dx_if
                    elseif (n == isec%nnode .and. isec%nnode > 1) then
                      dtdx = (t_hw(isec%nnode) - t_hw(isec%nnode-1)) / dx_if
                    else
                      dtdx = 0.0_dp
                    endif
                    v_loc = v_hw(n)
                    if (v_loc > 0.0_dp) then
                      dtdz_sq = 1.0_dp / (v_loc * v_loc) - dtdx * dtdx
                      if (dtdz_sq > 0.0_dp) then
                        dtdz = sqrt(dtdz_sq)
                      else
                        dtdz = 0.0_dp
                      endif
                    else
                      dtdz = 0.0_dp
                    endif
                    reg%time_gradient(1, m) = dtdx
                    reg%time_gradient(2, m) = dtdz
                  endif
                endif
              endif
            end do

            deallocate(t_hw, v_hw)
          endif
        endif

        if (associated(tf%arrivaltime)) then
          tf%arrivaltime(1:reg%nnode) = reg%arrivaltime(1:reg%nnode)
        end if
        if (associated(tf%time_gradient)) then
          tf%time_gradient(1:2,1:reg%nnode) = reg%time_gradient(1:2,1:reg%nnode)
        end if

        if (ipath == 3 .and. itf_idx == 2) then
          isec => tf%istart
          call transfer_times_to_intersection(isec, src%time_field(itf))
        end if

        if (ipath == 3 .and. itf_idx == 3) then
          isec => tf%istart
          if (associated(isec)) then
            do n = 1, isec%nnode
              m = 0
              if (associated(isec%rleft_node_id)) m = isec%rleft_node_id(n)
              if (m == 0 .and. associated(isec%rright_node_id)) m = isec%rright_node_id(n)
              if (m > 0 .and. m <= reg%nnode) then
                if (reg%node_status(m) > 0 .or. reg%node_status(m) == 0) then
                  dtdx = src%time_field(1)%time_gradient(1, m)
                  dtdz_orig = src%time_field(1)%time_gradient(2, m)
                  v_loc = isec%vel_top(VTYPE_P, n)
                  if (v_loc > 0.0_dp) then
                    dtdz_sq = 1.0_dp / (v_loc * v_loc) - dtdx * dtdx
                    if (dtdz_sq > 0.0_dp) then
                      reg%time_gradient(1, m) = dtdx
                      reg%time_gradient(2, m) = -sqrt(dtdz_sq)
                    else
                      reg%time_gradient(1, m) = dtdx
                      if (dtdz_orig < 0.0_dp) then
                        reg%time_gradient(2, m) = dtdz_orig
                      else
                        reg%time_gradient(2, m) = -abs(dtdz_orig)
                      endif
                    endif
                  endif
                endif
              endif
            end do
          endif

          call propagate2d(reg, vtype)

          do n = 1, isec%nnode
            m = 0
            if (associated(isec%rleft_node_id)) m = isec%rleft_node_id(n)
            if (m == 0 .and. associated(isec%rright_node_id)) m = isec%rright_node_id(n)
            if (m > 0 .and. m <= reg%nnode) then
              if (reg%node_status(m) == 0) then
                call fit_gradient_at_node2d(reg, m, reg%time_gradient(1,m), reg%time_gradient(2,m))
              endif
            endif
          end do

          if (associated(tf%arrivaltime)) then
            tf%arrivaltime(1:reg%nnode) = reg%arrivaltime(1:reg%nnode)
          end if
          if (associated(tf%time_gradient)) then
            tf%time_gradient(1:2,1:reg%nnode) = reg%time_gradient(1:2,1:reg%nnode)
          end if
        end if

      end do
    end do
    call system_clock(t_c1)
    t_fmm_acc = t_fmm_acc + real(t_c1 - t_c0, kind=dp) / real(t_rate, kind=dp)

    ! The FMM pass above writes the region time arrays directly, so whatever
    ! the ray tracer reads from now on has to be re-loaded. Clearing the
    ! bookkeeping forces the first load of each time field into its region.
    do ireg = 1, n_regions
      region(ireg)%loaded_tf => null()
    end do

    call system_clock(t_c0)
    do ircv = 1, n_receivers
      do iray = 1, receiver(ircv)%n_rays
        if (receiver(ircv)%ray(iray)%source_id /= isrc) cycle
        ipath = receiver(ircv)%ray(iray)%raypath_id
        itf = src%path(ipath)%tf_sequence(src%path(ipath)%n_tf)

        if (in_water(ircv)) then
          if (ipath == 4 .or. ipath == 5) then
            ! Pg_ss / Pg_ss2: path ends at sea surface (receiver position).
            ! Restore receiver to original water-column position so the
            ! final time field can be interpolated directly. No water_time needed.
            receiver(ircv)%x = rec_x_orig(ircv)
            receiver(ircv)%z = rec_z_orig(ircv)
            water_time(ircv) = 0.0_dp
          else
            if (use_entry_points_file) then
               call get_entry_point(ircv, src%id, ipath, &
                    receiver(ircv)%x, receiver(ircv)%z, water_time(ircv), found)
               if (.not. found) then
                  call find_optimal_entry_inv(src, ircv, ipath, itf)
               endif
            else
               call find_optimal_entry_inv(src, ircv, ipath, itf)
            endif
          endif
       endif
        call trace_ray_from_receiver2d(receiver(ircv), src, receiver(ircv)%ray(iray))

        if (in_water(ircv)) then
          ! Restore receiver to seafloor for other phases
          receiver(ircv)%x = rec_x_orig(ircv)
          call interpolate_iface_z_inv(intrface(2), rec_x_orig(ircv), &
               receiver(ircv)%z)
        endif

        if (receiver(ircv)%ray(iray)%valid) then
          receiver(ircv)%ray(iray)%receiver_time = receiver(ircv)%ray(iray)%receiver_time &
               + water_time(ircv)
        endif

        if (receiver(ircv)%ray(iray)%valid .and. .not. skip_frechet_in_forward) then
          call system_clock(t_c2)
          call compute_frechet_for_ray2d(receiver(ircv)%ray(iray))
          call system_clock(t_c3)
          t_fr_acc = t_fr_acc + real(t_c3 - t_c2, kind=dp) / real(t_rate, kind=dp)
        end if
      end do
    end do
    call system_clock(t_c1)
    t_trace_acc = t_trace_acc + real(t_c1 - t_c0, kind=dp) / real(t_rate, kind=dp)

  end do

  write(*,'(A,F9.2,A,F9.2,A,F9.2,A)') '    forward breakdown: FMM=', t_fmm_acc, &
    ' s, trace=', t_trace_acc - t_fr_acc, ' s, frechet=', t_fr_acc, ' s'
  flush(6)

end subroutine run_forward2d


subroutine compute_frechet_for_ray2d(ray)
  implicit none
  type(Tray2d) :: ray

  integer :: isec_idx, ipt, npts, ivx, ivz
  integer :: n, k, ig, idx
  integer :: nx_vg, nz_vg, vtype, reg_ivgrid
  real(kind=dp) :: px, pz, seg_len, vel
  real(kind=dp) :: wx, wz
  integer :: max_pdev, n_pdev
  integer, allocatable :: tmp_indx(:)
  real(kind=dp), allocatable :: tmp_pdev(:)
  real(kind=dp), allocatable :: cell_deriv(:)
  real(kind=dp) :: path_time
  type(Tvelocity_grid2d), pointer :: vg
  type(Tregion2d), pointer :: reg

  real(kind=dp) :: gradt_in(2), gradt_out(2), normal(2), h
  real(kind=dp) :: gradt_in_perp, gradt_out_perp, gradt_in_par
  real(kind=dp) :: vel_in, vel_out, geo_factor, det, norm_len
  type(Tintersection2d), pointer :: isec_if
  type(Tinterface2d), pointer :: iface
  type(Tray_section2d), pointer :: raysec, raysec_prev
  logical :: rec_region_done, do_interface_partials

  max_pdev = n_inv_parms

  allocate(tmp_indx(max_pdev), tmp_pdev(max_pdev))
  allocate(cell_deriv(n_inv_parms))
  cell_deriv = 0.0_dp
  n_pdev = 0
  path_time = 0.0_dp

  do isec_idx = 1, ray%nsections
    if (ray%section(isec_idx)%npoints <= 1) cycle
    reg => ray%section(isec_idx)%reg
    reg_ivgrid = reg%ivgrid

    ! Skip water-layer sections for Frechet derivative computation.
    ! Water velocity is fixed (not inverted), so water-layer ray paths
    ! should not contribute to crust velocity sensitivity.
    ! The bilinear interpolation leaks weight to crust nodes (inv_idx>0)
    ! at the seafloor boundary, giving incorrect Frechet derivatives.
    if (reg%id == 1) cycle  ! region 1 = water layer

    do ig = 1, n_inv_vgrid
      if (vgrids_to_be_inv(ig)%igrid /= reg_ivgrid) cycle

      vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
      vtype = vgrids_to_be_inv(ig)%vtype
      nx_vg = vg%nx
      nz_vg = vg%nz

      npts = ray%section(isec_idx)%npoints

      do ipt = 1, npts - 1
        px = 0.5_dp * (ray%section(isec_idx)%point(1,ipt) + ray%section(isec_idx)%point(1,ipt+1))
        pz = 0.5_dp * (ray%section(isec_idx)%point(2,ipt) + ray%section(isec_idx)%point(2,ipt+1))

        seg_len = sqrt((ray%section(isec_idx)%point(1,ipt+1) - ray%section(isec_idx)%point(1,ipt))**2 + &
                       (ray%section(isec_idx)%point(2,ipt+1) - ray%section(isec_idx)%point(2,ipt))**2)

        if (seg_len < 1.0e-10_dp) cycle

        ivx = floor((px - vg%x(1)) / vg%dx) + 1
        if (ivx < 1) ivx = 1
        if (ivx >= nx_vg) ivx = nx_vg - 1

        if (pz < vg%z(1)) then
           cycle
        else if (pz >= vg%z(nz_vg)) then
           ivz = nz_vg - 1
           wz = 1.0_dp
        else
           ivz = floor((pz - vg%z(1)) / vg%dz) + 1
           if (ivz < 1) ivz = 1
           if (ivz >= nz_vg) ivz = nz_vg - 1
           wz = (pz - vg%z(ivz)) / vg%dz
        endif

        if (ivx < 1) ivx = 1
        if (ivx >= nx_vg) ivx = nx_vg - 1
        if (ivz < 1) ivz = 1
        if (ivz >= nz_vg) ivz = nz_vg - 1

        wx = (px - vg%x(ivx)) / vg%dx
        if (wx < 0.0_dp) wx = 0.0_dp
        if (wx > 1.0_dp) wx = 1.0_dp
        if (wz < 0.0_dp) wz = 0.0_dp
        if (wz > 1.0_dp) wz = 1.0_dp

        vel = (1.0_dp-wx)*(1.0_dp-wz)*vg%velocity(vtype,ivx,ivz) + &
              wx*(1.0_dp-wz)*vg%velocity(vtype,ivx+1,ivz) + &
              (1.0_dp-wx)*wz*vg%velocity(vtype,ivx,ivz+1) + &
              wx*wz*vg%velocity(vtype,ivx+1,ivz+1)

        if (vel < 1.0e-6_dp) vel = 1.0e-6_dp

        path_time = path_time + seg_len / vel

        idx = vg%inv_idx(ivx, ivz)
        if (idx > 0) then
          cell_deriv(idx) = cell_deriv(idx) + (1.0_dp-wx)*(1.0_dp-wz) * seg_len
        end if

        idx = vg%inv_idx(ivx+1, ivz)
        if (idx > 0) then
          cell_deriv(idx) = cell_deriv(idx) + wx*(1.0_dp-wz) * seg_len
        end if

        idx = vg%inv_idx(ivx, ivz+1)
        if (idx > 0) then
          cell_deriv(idx) = cell_deriv(idx) + (1.0_dp-wx)*wz * seg_len
        end if

        idx = vg%inv_idx(ivx+1, ivz+1)
        if (idx > 0) then
          cell_deriv(idx) = cell_deriv(idx) + wx*wz * seg_len
        end if

      end do
    end do



  end do

  nullify(raysec_prev)
  rec_region_done = .false.
  do_interface_partials = .false.

  do isec_idx = ray%nsections, 1, -1

     raysec => ray%section(isec_idx)

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

     if (do_interface_partials .and. n_inv_iface > 0) then

        isec_if => raysec%istart
        if (.not. associated(isec_if)) then
           raysec_prev => raysec
           cycle
        endif

        iface => intrface(isec_if%iface_id)

        if (iface%to_be_inverted) then

           px = raysec%point(1,1)
           pz = raysec%point(2,1)

           call interpolate_time_gradient2d(raysec%reg, px, pz, gradt_in(1), gradt_in(2))

           vel_in = get_velocity_at2d(raysec%reg, px, pz, VTYPE_P)
           norm_len = sqrt(sum(gradt_in**2))
           if (norm_len > 0.0_dp) then
              gradt_in = gradt_in / (vel_in * norm_len)
           endif

           call interface_normal2d(px, isec_if, normal(1), normal(2), h)

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

           idx = iface%start_index + ivx - 1
           if (idx >= 1 .and. idx <= n_inv_parms) then
             cell_deriv(idx + 1) = cell_deriv(idx + 1) + (1.0_dp - wx)*geo_factor
           end if

           idx = iface%start_index + ivx
           if (idx >= 1 .and. idx <= n_inv_parms) then
             cell_deriv(idx + 1) = cell_deriv(idx + 1) + wx*geo_factor
           end if

        endif

     endif

     raysec_prev => raysec

  end do

  n_pdev = 0
  do n = 1, n_inv_parms
    if (abs(cell_deriv(n)) > 1.0e-30_dp) then
      n_pdev = n_pdev + 1
      tmp_indx(n_pdev) = n
      tmp_pdev(n_pdev) = cell_deriv(n)
    end if
  end do

  deallocate(cell_deriv)

  ray%n_pdev = n_pdev
  if (associated(ray%pdev_indx)) deallocate(ray%pdev_indx)
  if (associated(ray%pdev)) deallocate(ray%pdev)
  if (n_pdev > 0) then
    allocate(ray%pdev_indx(n_pdev), ray%pdev(n_pdev))
    ray%pdev_indx(1:n_pdev) = tmp_indx(1:n_pdev)
    ray%pdev(1:n_pdev) = tmp_pdev(1:n_pdev)
  else
    allocate(ray%pdev_indx(1), ray%pdev(1))
    ray%pdev_indx = 0
    ray%pdev = 0.0_dp
  end if

  deallocate(tmp_indx, tmp_pdev)

end subroutine compute_frechet_for_ray2d


subroutine gauss_solve2d(n, a, b, x, info)
  implicit none
  integer, intent(in) :: n
  real(kind=dp), intent(in) :: a(n, n), b(n)
  real(kind=dp), intent(out) :: x(n)
  integer, intent(out) :: info

  integer :: i, j, k, pivot_row
  real(kind=dp) :: max_val, factor, temp
  real(kind=dp), allocatable :: aug(:,:)

  info = 0

  if (n <= 0) then
    info = -1
    return
  end if

  allocate(aug(n, n+1))

  do i = 1, n
    do j = 1, n
      aug(i,j) = a(i,j)
    end do
    aug(i,n+1) = b(i)
  end do

  do k = 1, n

    max_val = abs(aug(k,k))
    pivot_row = k
    do i = k+1, n
      if (abs(aug(i,k)) > max_val) then
        max_val = abs(aug(i,k))
        pivot_row = i
      end if
    end do

    if (max_val < 1.0e-30_dp) then
      info = k
      deallocate(aug)
      return
    end if

    if (pivot_row /= k) then
      do j = k, n+1
        temp = aug(k,j)
        aug(k,j) = aug(pivot_row,j)
        aug(pivot_row,j) = temp
      end do
    end if

    do i = k+1, n
      factor = aug(i,k) / aug(k,k)
      do j = k, n+1
        aug(i,j) = aug(i,j) - factor * aug(k,j)
      end do
    end do

  end do

  x(n) = aug(n,n+1) / aug(n,n)
  do i = n-1, 1, -1
    x(i) = aug(i,n+1)
    do j = i+1, n
      x(i) = x(i) - aug(i,j) * x(j)
    end do
    x(i) = x(i) / aug(i,i)
  end do

  deallocate(aug)

end subroutine gauss_solve2d


subroutine read_observed_times2d(filename)
  implicit none
  character(len=*), intent(in) :: filename

  integer :: iunit, ios, i, nobs_file
  integer :: ircv, isrc, ipath, iray
  real(kind=dp) :: t_obs

  iunit = 30
  open(unit=iunit, file=trim(filename), status='old', action='read', iostat=ios)
  if (ios /= 0) then
    write(*,'(A,A)') 'Warning: cannot open observed times file: ', trim(filename)
    return
  end if

  read(iunit, *, iostat=ios) nobs_file
  if (ios /= 0 .or. nobs_file <= 0) then
    close(iunit)
    return
  end if

  do ircv = 1, n_receivers
    do iray = 1, receiver(ircv)%n_rays
      receiver(ircv)%ray(iray)%obs_time = -1.0_dp
    end do
  end do

  inv_nobs = 0
  do i = 1, nobs_file
    read(iunit, *, iostat=ios) ircv, isrc, ipath, t_obs
    if (ios /= 0) then
      exit
    end if
    if (ircv >= 1 .and. ircv <= n_receivers) then
      do iray = 1, receiver(ircv)%n_rays
        if (receiver(ircv)%ray(iray)%source_id == isrc .and. &
            receiver(ircv)%ray(iray)%raypath_id == ipath) then
          receiver(ircv)%ray(iray)%obs_time = t_obs
          inv_nobs = inv_nobs + 1
          exit
        end if
      end do
    end if
  end do

  close(iunit)

  write(*,'(A,I0,A,I0)') '  Read observed times: nobs_file=', nobs_file, ' matched=', inv_nobs

end subroutine read_observed_times2d


subroutine build_augmented_system(nparm, ndata, nnz, gindx, gval, residual, row_ptr, &
  lam_vel, smooth_scale, n_aug_rows, n_aug_nnz, &
  aug_row, aug_col, aug_val, aug_rhs)
  implicit none
  integer, intent(in) :: nparm, ndata, nnz
  integer, intent(in) :: gindx(nnz), row_ptr(ndata+1)
  real(kind=dp), intent(in) :: gval(nnz), residual(ndata)
  real(kind=dp), intent(in) :: lam_vel, smooth_scale
  integer, intent(out) :: n_aug_rows, n_aug_nnz
  integer, allocatable, intent(out) :: aug_row(:), aug_col(:)
  real(kind=dp), allocatable, intent(out) :: aug_val(:), aug_rhs(:)

  integer :: ig, ix, iz, idx, vtype, nx_vg, nz_vg
  integer :: n_damp, n_smooth, n_smooth_nnz, i, j, k, row_off, n_neighbors
  integer, allocatable :: smooth_row_tmp(:), smooth_col_tmp(:)
  real(kind=dp), allocatable :: smooth_val_tmp(:)
  type(Tvelocity_grid2d), pointer :: vg

  n_damp = 0
  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx <= 0) cycle
        n_damp = n_damp + 1
      end do
    end do
  end do

  n_smooth = 0
  n_smooth_nnz = 0
  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx <= 0) cycle
        n_neighbors = 0
        if (vg%inv_idx(ix-1, iz) > 0) n_neighbors = n_neighbors + 1
        if (vg%inv_idx(ix+1, iz) > 0) n_neighbors = n_neighbors + 1
        if (vg%inv_idx(ix, iz-1) > 0) n_neighbors = n_neighbors + 1
        if (vg%inv_idx(ix, iz+1) > 0) n_neighbors = n_neighbors + 1
        if (n_neighbors > 0) then
          n_smooth = n_smooth + 1
          n_smooth_nnz = n_smooth_nnz + 1 + n_neighbors
        end if
      end do
    end do
  end do

  n_aug_rows = ndata + n_damp + n_smooth
  n_aug_nnz = nnz + n_damp + n_smooth_nnz

  allocate(aug_row(n_aug_nnz), aug_col(n_aug_nnz), aug_val(n_aug_nnz))
  allocate(aug_rhs(n_aug_rows))

  j = 0
  do i = 1, ndata
    do k = row_ptr(i), row_ptr(i+1) - 1
      j = j + 1
      aug_row(j) = i
      aug_col(j) = gindx(k)
      aug_val(j) = gval(k)
    end do
  end do

  do i = 1, ndata
    aug_rhs(i) = residual(i)
  end do

  row_off = ndata

  j = nnz
  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx <= 0) cycle
        row_off = row_off + 1
        j = j + 1
        aug_row(j) = row_off
        aug_col(j) = idx
        aug_val(j) = lam_vel
        aug_rhs(row_off) = 0.0_dp
      end do
    end do
  end do

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx <= 0) cycle
        n_neighbors = 0
        if (vg%inv_idx(ix-1, iz) > 0) n_neighbors = n_neighbors + 1
        if (vg%inv_idx(ix+1, iz) > 0) n_neighbors = n_neighbors + 1
        if (vg%inv_idx(ix, iz-1) > 0) n_neighbors = n_neighbors + 1
        if (vg%inv_idx(ix, iz+1) > 0) n_neighbors = n_neighbors + 1
        if (n_neighbors > 0) then
          row_off = row_off + 1
          j = j + 1
          aug_row(j) = row_off
          aug_col(j) = idx
          aug_val(j) = -n_neighbors * smooth_scale
          aug_rhs(row_off) = 0.0_dp
          if (vg%inv_idx(ix-1, iz) > 0) then
            j = j + 1
            aug_row(j) = row_off
            aug_col(j) = vg%inv_idx(ix-1, iz)
            aug_val(j) = smooth_scale
          end if
          if (vg%inv_idx(ix+1, iz) > 0) then
            j = j + 1
            aug_row(j) = row_off
            aug_col(j) = vg%inv_idx(ix+1, iz)
            aug_val(j) = smooth_scale
          end if
          if (vg%inv_idx(ix, iz-1) > 0) then
            j = j + 1
            aug_row(j) = row_off
            aug_col(j) = vg%inv_idx(ix, iz-1)
            aug_val(j) = smooth_scale
          end if
          if (vg%inv_idx(ix, iz+1) > 0) then
            j = j + 1
            aug_row(j) = row_off
            aug_col(j) = vg%inv_idx(ix, iz+1)
            aug_val(j) = smooth_scale
          end if
        end if
      end do
    end do
  end do

  n_aug_nnz = j

  write(*,'(A,I8,A,I8,A,I8,A,I8)') '  Augmented system: ndata=', ndata, &
    ' n_damp=', n_damp, ' n_smooth=', n_smooth, ' nnz=', n_aug_nnz

end subroutine build_augmented_system


subroutine update_model_vel(delta_m, step)
  implicit none
  real(kind=dp), intent(in) :: delta_m(n_inv_parms)
  real(kind=dp), intent(in) :: step

  integer :: ig, ix, iz, idx, vtype, nx_vg, nz_vg
  integer :: n_updated, n_skipped
  real(kind=dp) :: old_v, new_v, max_change
  type(Tvelocity_grid2d), pointer :: vg

  n_updated = 0
  n_skipped = 0
  max_change = 0.0_dp

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz

    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0) then
          old_v = vg%velocity(vtype, ix, iz)
          vg%velocity(vtype, ix, iz) = vg%velocity(vtype, ix, iz) + step * delta_m(idx)
          if (vg%velocity(vtype, ix, iz) < 2.0_dp) vg%velocity(vtype, ix, iz) = 2.0_dp
          if (vg%velocity(vtype, ix, iz) > 10.0_dp) vg%velocity(vtype, ix, iz) = 10.0_dp
          new_v = vg%velocity(vtype, ix, iz)
          if (abs(new_v - old_v) > max_change) max_change = abs(new_v - old_v)
          if (abs(new_v - old_v) > 1.0e-6_dp) n_updated = n_updated + 1
        end if
      end do
    end do
  end do

  write(*,'(A,I6,A,I6,A,F12.6)') '  update_model_vel: updated=', n_updated, &
    ' skipped_water=', n_skipped, ' max_dv=', max_change

end subroutine update_model_vel


subroutine lsqr_damped_smooth_solver(ndata, nparm, nnz, gindx, gval, &
  row_ptr, b, lam_vel, lam_smooth, mvscale, dmodel_total, x, max_iter)

  implicit none
  integer, intent(in) :: ndata, nparm, nnz, max_iter
  integer, intent(in) :: gindx(nnz), row_ptr(ndata+1)
  real(kind=dp), intent(in) :: gval(nnz), b(ndata)
  real(kind=dp), intent(in) :: lam_vel, lam_smooth
  real(kind=dp), intent(in) :: mvscale(nparm)
  real(kind=dp), intent(in) :: dmodel_total(nparm)
  real(kind=dp), intent(inout) :: x(nparm)

  integer :: ndata_aug, nnz_aug, irow, innz, i, j, k, idx
  integer :: ig, ix, iz, vtype, nx_vg, nz_vg
  integer :: n_damp, n_smooth_rows, n_smooth_nnz, n_damp_nnz
  integer :: ii, kk, jdx, n_h_neighbors, n_v_neighbors, jdx1, jdx2, jdx3
  real(kind=dp) :: Lh, Lv, Lh2, Lv2, dx_node, dz_node, beta_wt, beta_sum
  real(kind=dp) :: cell_area
  integer, allocatable :: gindx_aug(:), row_ptr_aug(:)
  real(kind=dp), allocatable :: gval_aug(:), b_aug(:)
  real(kind=dp), allocatable :: u(:), v(:), w(:), Av(:), Atu(:), z(:)
  real(kind=dp), allocatable :: gtg_diag(:)
  real(kind=dp) :: beta, alpha, c_bar, s_bar, phi_bar, rho_bar
  real(kind=dp) :: rho, theta, c, s, phi, rnorm, bnorm
  real(kind=dp) :: Anorm2, Anorm, Arnorm, test2, lsqr_atol
  integer :: iter
  type(Tvelocity_grid2d), pointer :: vg

  Lh = 114.0_dp
  Lv = 3.0_dp
  Lh2 = Lh * Lh
  Lv2 = Lv * Lv

  allocate(gtg_diag(nparm))
  gtg_diag = 0.0_dp
  do i = 1, ndata
    do k = row_ptr(i), row_ptr(i+1) - 1
      j = gindx(k)
      if (j >= 1 .and. j <= nparm) then
        gtg_diag(j) = gtg_diag(j) + gval(k) * gval(k)
      end if
    end do
  end do

  if (lam_vel <= 0.0_dp .and. lam_smooth <= 0.0_dp) then
    ndata_aug = ndata
    nnz_aug = nnz
    allocate(gindx_aug(nnz_aug), gval_aug(nnz_aug))
    allocate(b_aug(ndata_aug), row_ptr_aug(ndata_aug+1))
    gindx_aug(1:nnz) = gindx(1:nnz)
    gval_aug(1:nnz) = gval(1:nnz)
    b_aug(1:ndata) = b(1:ndata)
    row_ptr_aug(1:ndata+1) = row_ptr(1:ndata+1)
    irow = ndata
    innz = nnz
    n_damp = 0
    n_smooth_rows = 0
  else

  n_damp = 0
  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx <= 0) cycle
        if (gtg_diag(idx) > 0.0_dp) then
          n_damp = n_damp + 1
        end if
      end do
    end do
  end do

  n_smooth_rows = 0
  n_smooth_nnz = 0
  if (lam_smooth > 0.0_dp) then
    do ig = 1, n_inv_vgrid
      vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
      vtype = vgrids_to_be_inv(ig)%vtype
      nx_vg = vg%nx
      nz_vg = vg%nz
      do iz = 1, nz_vg
        do ix = 1, nx_vg
          idx = vg%inv_idx(ix, iz)
          if (idx <= 0) cycle
          if (gtg_diag(idx) <= 0.0_dp) cycle

          n_h_neighbors = 0
          do ii = 1, nx_vg
            if (ii == ix) cycle
            jdx = vg%inv_idx(ii, iz)
            if (jdx <= 0) cycle
            dx_node = abs(dble(ii - ix)) * vg%dx
            if (dx_node <= Lh) n_h_neighbors = n_h_neighbors + 1
          end do
          if (n_h_neighbors > 0) then
            n_smooth_rows = n_smooth_rows + 1
            n_smooth_nnz = n_smooth_nnz + 1 + n_h_neighbors
          end if

          n_v_neighbors = 0
          do kk = 1, nz_vg
            if (kk == iz) cycle
            jdx = vg%inv_idx(ix, kk)
            if (jdx <= 0) cycle
            dz_node = abs(dble(kk - iz)) * vg%dz
            if (dz_node <= Lv) n_v_neighbors = n_v_neighbors + 1
          end do
          if (n_v_neighbors > 0) then
            n_smooth_rows = n_smooth_rows + 1
            n_smooth_nnz = n_smooth_nnz + 1 + n_v_neighbors
          end if
        end do
      end do
    end do
  end if

  ndata_aug = ndata + n_damp + n_smooth_rows
  nnz_aug = nnz + n_damp + n_smooth_nnz

  allocate(gindx_aug(nnz_aug), gval_aug(nnz_aug))
  allocate(b_aug(ndata_aug), row_ptr_aug(ndata_aug+1))

  gindx_aug(1:nnz) = gindx(1:nnz)
  gval_aug(1:nnz) = gval(1:nnz)
  b_aug(1:ndata) = b(1:ndata)
  row_ptr_aug(1:ndata+1) = row_ptr(1:ndata+1)

  innz = nnz
  irow = ndata

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx <= 0) cycle
        if (gtg_diag(idx) <= 0.0_dp) cycle
        irow = irow + 1
        innz = innz + 1
        gindx_aug(innz) = idx
        gval_aug(innz) = lam_vel
        b_aug(irow) = 0.0_dp
        row_ptr_aug(irow) = innz
      end do
    end do
  end do

  if (lam_smooth > 0.0_dp) then
    do ig = 1, n_inv_vgrid
      vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
      vtype = vgrids_to_be_inv(ig)%vtype
      nx_vg = vg%nx
      nz_vg = vg%nz
      do iz = 1, nz_vg
        do ix = 1, nx_vg
          idx = vg%inv_idx(ix, iz)
          if (idx <= 0) cycle
          if (gtg_diag(idx) <= 0.0_dp) cycle

          n_h_neighbors = 0
          beta_sum = 0.0_dp
          do ii = 1, nx_vg
            if (ii == ix) cycle
            jdx = vg%inv_idx(ii, iz)
            if (jdx <= 0) cycle
            dx_node = abs(dble(ii - ix)) * vg%dx
            if (dx_node <= Lh) then
              n_h_neighbors = n_h_neighbors + 1
              beta_sum = beta_sum + exp(-dx_node*dx_node/Lh2)
            end if
          end do

          if (n_h_neighbors > 0) then
            irow = irow + 1
            row_ptr_aug(irow) = innz + 1
            innz = innz + 1
            gindx_aug(innz) = idx
            gval_aug(innz) = -beta_sum * lam_smooth
            do ii = 1, nx_vg
              if (ii == ix) cycle
              jdx = vg%inv_idx(ii, iz)
              if (jdx <= 0) cycle
              dx_node = abs(dble(ii - ix)) * vg%dx
              if (dx_node <= Lh) then
                beta_wt = exp(-dx_node*dx_node/Lh2)
                innz = innz + 1
                gindx_aug(innz) = jdx
                gval_aug(innz) = beta_wt * lam_smooth
              end if
            end do
            b_aug(irow) = 0.0_dp
          end if

          n_v_neighbors = 0
          beta_sum = 0.0_dp
          do kk = 1, nz_vg
            if (kk == iz) cycle
            jdx = vg%inv_idx(ix, kk)
            if (jdx <= 0) cycle
            dz_node = abs(dble(kk - iz)) * vg%dz
            if (dz_node <= Lv) then
              n_v_neighbors = n_v_neighbors + 1
              beta_sum = beta_sum + exp(-dz_node*dz_node/Lv2)
            end if
          end do

          if (n_v_neighbors > 0) then
            irow = irow + 1
            row_ptr_aug(irow) = innz + 1
            innz = innz + 1
            gindx_aug(innz) = idx
            gval_aug(innz) = -beta_sum * lam_smooth
            do kk = 1, nz_vg
              if (kk == iz) cycle
              jdx = vg%inv_idx(ix, kk)
              if (jdx <= 0) cycle
              dz_node = abs(dble(kk - iz)) * vg%dz
              if (dz_node <= Lv) then
                beta_wt = exp(-dz_node*dz_node/Lv2)
                innz = innz + 1
                gindx_aug(innz) = jdx
                gval_aug(innz) = beta_wt * lam_smooth
              end if
            end do
            b_aug(irow) = 0.0_dp
          end if
        end do
      end do
    end do
  end if

  row_ptr_aug(irow+1) = innz + 1

  write(*,'(A,I6,A,I10,A,I6,A,I8)') '    augmented system: ndata_aug=', irow, &
    ' nnz_aug=', innz, ' n_damp=', n_damp, ' n_smooth=', n_smooth_rows

  end if

  bnorm = sqrt(dot_product(b_aug(1:irow), b_aug(1:irow)))
  if (bnorm < 1.0e-30_dp) then
    x = 0.0_dp
    deallocate(gindx_aug, gval_aug, b_aug, row_ptr_aug)
    return
  end if

  allocate(u(irow), v(nparm), w(nparm), Av(irow), Atu(nparm), z(nparm))

  x = 0.0_dp

  u(1:irow) = b_aug(1:irow)
  beta = sqrt(dot_product(u, u))
  if (beta > 0.0_dp) then
    u = u / beta
  else
    x = 0.0_dp
    deallocate(u, v, w, Av, Atu, z, gindx_aug, gval_aug, b_aug, row_ptr_aug)
    return
  end if

  call sparse_ATx(irow, nparm, innz, gindx_aug, gval_aug, row_ptr_aug, u, v)
  alpha = sqrt(dot_product(v, v))
  if (alpha > 0.0_dp) then
    v = v / alpha
  else
    x = 0.0_dp
    deallocate(u, v, w, Av, Atu, z, gindx_aug, gval_aug, b_aug, row_ptr_aug)
    return
  end if

  w = v
  c_bar = 1.0_dp
  s_bar = 0.0_dp
  phi_bar = beta
  rho_bar = alpha
  Anorm2 = alpha**2 + beta**2
  lsqr_atol = 1.0e-8_dp

  do iter = 1, max_iter
    call sparse_Ax(irow, nparm, innz, gindx_aug, gval_aug, row_ptr_aug, v, Av)
    u = Av - alpha * u
    beta = sqrt(dot_product(u, u))
    if (beta > 0.0_dp) u = u / beta

    call sparse_ATx(irow, nparm, innz, gindx_aug, gval_aug, row_ptr_aug, u, Atu)
    v = Atu - beta * v
    alpha = sqrt(dot_product(v, v))
    if (alpha > 0.0_dp) v = v / alpha

    Anorm2 = Anorm2 + alpha**2 + beta**2

    rho = sqrt(rho_bar**2 + beta**2)
    c = rho_bar / rho
    s = beta / rho
    theta = s * alpha
    rho_bar = -c * alpha
    phi = c * phi_bar
    phi_bar = s * phi_bar

    z = w
    w = v - (theta / rho) * w
    x = x + (phi / rho) * z

    rnorm = abs(phi_bar)
    Arnorm = abs(phi_bar * alpha * c)
    Anorm = sqrt(Anorm2)

    if (Anorm > 0.0_dp .and. rnorm > 0.0_dp) then
      test2 = Arnorm / (Anorm * rnorm)
    else
      test2 = 0.0_dp
    end if

    if (mod(iter, 20) == 0) then
      write(*,'(A,I4,A,F12.6,A,F12.6,A,F12.6)') '    LSQR iter=', iter, &
        ' |r|/|b|=', rnorm/bnorm, ' test2=', test2, ' |x|=', sqrt(dot_product(x,x))
    end if

    if (test2 < lsqr_atol) exit
  end do

  write(*,'(A,I4,A,F12.6,A,F12.6,A,F12.6)') '    LSQR: ', iter, &
    ' iterations, |r|/|b|=', rnorm/bnorm, ' test2=', test2, ' |x|=', sqrt(dot_product(x,x))

  deallocate(u, v, w, Av, Atu, z)
  deallocate(gindx_aug, gval_aug, b_aug, row_ptr_aug)
  deallocate(gtg_diag)

end subroutine lsqr_damped_smooth_solver


subroutine apply_gaussian_smooth(nparm, delta_m, lam_smooth, mvscale)
  implicit none
  integer, intent(in) :: nparm
  real(kind=dp), intent(inout) :: delta_m(nparm)
  real(kind=dp), intent(in) :: lam_smooth, mvscale(nparm)

  integer :: ig, ix, iz, vtype, nx_vg, nz_vg, idx, ii, kk, jdx
  real(kind=dp) :: Lh, Lv, Lh2, Lv2, dx_node, dz_node, beta_wt, beta_sum
  real(kind=dp), allocatable :: smooth_x(:), smooth_z(:)
  real(kind=dp) :: smooth_factor
  type(Tvelocity_grid2d), pointer :: vg

  Lh = 114.0_dp
  Lv = 3.0_dp
  Lh2 = Lh * Lh
  Lv2 = Lv * Lv
  smooth_factor = lam_smooth / (1.0_dp + lam_smooth)

  allocate(smooth_x(nparm), smooth_z(nparm))
  smooth_x = 0.0_dp
  smooth_z = 0.0_dp

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx <= 0) cycle

        beta_sum = 0.0_dp
        do ii = 1, nx_vg
          if (ii == ix) cycle
          jdx = vg%inv_idx(ii, iz)
          if (jdx <= 0) cycle
          dx_node = abs(dble(ii - ix)) * vg%dx
          if (dx_node <= Lh) then
            beta_wt = exp(-dx_node*dx_node/Lh2)
            smooth_x(idx) = smooth_x(idx) + beta_wt * delta_m(jdx) * mvscale(jdx)
            beta_sum = beta_sum + beta_wt
          end if
        end do
        if (beta_sum > 0.0_dp) then
          smooth_x(idx) = (smooth_x(idx) / (beta_sum * mvscale(idx))) - delta_m(idx)
        end if

        beta_sum = 0.0_dp
        do kk = 1, nz_vg
          if (kk == iz) cycle
          jdx = vg%inv_idx(ix, kk)
          if (jdx <= 0) cycle
          dz_node = abs(dble(kk - iz)) * vg%dz
          if (dz_node <= Lv) then
            beta_wt = exp(-dz_node*dz_node/Lv2)
            smooth_z(idx) = smooth_z(idx) + beta_wt * delta_m(jdx) * mvscale(jdx)
            beta_sum = beta_sum + beta_wt
          end if
        end do
        if (beta_sum > 0.0_dp) then
          smooth_z(idx) = (smooth_z(idx) / (beta_sum * mvscale(idx))) - delta_m(idx)
        end if
      end do
    end do
  end do

  do jdx = 1, nparm
    delta_m(jdx) = delta_m(jdx) + smooth_factor * (smooth_x(jdx) + smooth_z(jdx))
  end do

  deallocate(smooth_x, smooth_z)

  write(*,'(A,F8.4)') '    Gaussian smooth applied, factor=', smooth_factor

end subroutine apply_gaussian_smooth


subroutine update_model_direct(delta_m, step, mvscale)
  implicit none
  real(kind=dp), intent(in) :: delta_m(n_inv_parms)
  real(kind=dp), intent(in) :: step
  real(kind=dp), intent(in) :: mvscale(n_inv_parms)

  integer :: ig, ix, iz, idx, vtype, nx_vg, nz_vg
  integer :: n_updated, n_skipped
  real(kind=dp) :: old_v, new_v, max_change
  type(Tvelocity_grid2d), pointer :: vg

  n_updated = 0
  n_skipped = 0
  max_change = 0.0_dp

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz

    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0) then
          old_v = vg%velocity(vtype, ix, iz)
          new_v = old_v + step * delta_m(idx)
          if (new_v < 2.0_dp) new_v = 2.0_dp
          if (new_v > 10.0_dp) new_v = 10.0_dp
          vg%velocity(vtype, ix, iz) = new_v
          new_v = vg%velocity(vtype, ix, iz)
          if (abs(new_v - old_v) > max_change) max_change = abs(new_v - old_v)
          if (abs(new_v - old_v) > 1.0e-6_dp) n_updated = n_updated + 1
        end if
      end do
    end do
  end do

  write(*,'(A,I6,A,I6,A,F12.6)') '  update_model_direct: updated=', n_updated, &
    ' skipped_water=', n_skipped, ' max_dv=', max_change

end subroutine update_model_direct


subroutine save_model2d_old()
  implicit none
  integer :: ig, ix, iz, idx, vtype, nx_vg, nz_vg
  type(Tvelocity_grid2d), pointer :: vg

  if (.not. allocated(saved_vel)) then
    allocate(saved_vel(n_inv_parms))
  end if

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0) then
          saved_vel(idx) = vg%velocity(vtype, ix, iz)
        end if
      end do
    end do
  end do

end subroutine save_model2d_old


subroutine restore_model2d_old()
  implicit none
  integer :: ig, ix, iz, idx, vtype, nx_vg, nz_vg
  type(Tvelocity_grid2d), pointer :: vg

  if (.not. allocated(saved_vel)) return

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0) then
          vg%velocity(vtype, ix, iz) = saved_vel(idx)
        end if
      end do
    end do
  end do

end subroutine restore_model2d_old


function compute_actual_rms() result(rms)
  implicit none
  real(kind=dp) :: rms

  integer :: ircv, iray, ndata, n_pos, n_neg
  real(kind=dp) :: res, sumsq, sum_res, avg_res

  call run_forward2d()

  ndata = 0
  sumsq = 0.0_dp
  sum_res = 0.0_dp
  n_pos = 0
  n_neg = 0
  do ircv = 1, n_receivers
    do iray = 1, receiver(ircv)%n_rays
      if (receiver(ircv)%ray(iray)%valid .and. &
          receiver(ircv)%ray(iray)%obs_time >= 0.0_dp) then
        res = receiver(ircv)%ray(iray)%obs_time - receiver(ircv)%ray(iray)%receiver_time
        if (abs(res) > inv_crit_chi) cycle
        sumsq = sumsq + res*res
        sum_res = sum_res + res
        ndata = ndata + 1
        if (res > 0) n_pos = n_pos + 1
        if (res < 0) n_neg = n_neg + 1
      end if
    end do
  end do

  if (ndata > 0) then
    rms = sqrt(sumsq / ndata)
    avg_res = sum_res / ndata
    write(*,'(A,I6,A,F12.6,A,F12.6,A,I5,A,I5)') '    compute_rms: ndata=', ndata, &
      ' rms=', rms, ' avg_res=', avg_res, ' n_pos=', n_pos, ' n_neg=', n_neg
  else
    rms = huge_time
  end if

end function compute_actual_rms


subroutine fd_verify_frechet()
  implicit none
  integer :: ircv, iray, n, k, j, ipt, isec, ivx, ivz
  integer :: nx_vg, nz_vg, vtype, reg_ivgrid, npts
  integer :: idx_test, ivx_test, iz_test, n_test, n_match
  integer :: ix, iz
  integer :: ig_test
  real(kind=dp) :: px, pz, seg_len, vel, vel_p, s_old, s_new, ds
  real(kind=dp) :: wx, wz, t_old, t_new, dt_fd, dt_pred
  real(kind=dp) :: sum_ratio
  type(Tray2d), pointer :: ray
  type(Tvelocity_grid2d), pointer :: vg
  type(Tregion2d), pointer :: reg

  ds = 1.0e-4_dp
  n_test = 0
  n_match = 0
  sum_ratio = 0.0_dp

  write(*,'(A)') '  FD verification of Frechet derivatives (ds=1e-4 s/km):'

  do ircv = 1, min(n_receivers, 20)
    do iray = 1, receiver(ircv)%n_rays
      if (.not. receiver(ircv)%ray(iray)%valid) cycle
      if (receiver(ircv)%ray(iray)%obs_time < 0.0_dp) cycle
      if (receiver(ircv)%ray(iray)%n_pdev <= 0) cycle

      ray => receiver(ircv)%ray(iray)

      idx_test = ray%pdev_indx(1)
      if (idx_test < 1 .or. idx_test > n_inv_parms) cycle

      ig_test = 1
      vg => vgrid(vgrids_to_be_inv(ig_test)%igrid, vgrids_to_be_inv(ig_test)%vtype)
      vtype = vgrids_to_be_inv(ig_test)%vtype
      nx_vg = vg%nx
      nz_vg = vg%nz

      ivx_test = 0
      iz_test = 0
      do iz = 1, nz_vg
        do ix = 1, nx_vg
          if (vg%inv_idx(ix, iz) == idx_test) then
            ivx_test = ix
            iz_test = iz
            exit
          end if
        end do
        if (ivx_test > 0) exit
      end do

      if (ivx_test < 1 .or. ivx_test > nx_vg) cycle
      if (iz_test < 1 .or. iz_test > nz_vg) cycle

      t_old = 0.0_dp
      t_new = 0.0_dp

      do isec = 1, ray%nsections
        if (ray%section(isec)%npoints <= 1) cycle
        reg => ray%section(isec)%reg
        reg_ivgrid = reg%ivgrid

        do k = 1, n_inv_vgrid
          if (vgrids_to_be_inv(k)%igrid /= reg_ivgrid) cycle
          vg => vgrid(vgrids_to_be_inv(k)%igrid, vgrids_to_be_inv(k)%vtype)
          vtype = vgrids_to_be_inv(k)%vtype
          nx_vg = vg%nx
          nz_vg = vg%nz

          npts = ray%section(isec)%npoints
          do ipt = 1, npts - 1
            px = 0.5_dp * (ray%section(isec)%point(1,ipt) + ray%section(isec)%point(1,ipt+1))
            pz = 0.5_dp * (ray%section(isec)%point(2,ipt) + ray%section(isec)%point(2,ipt+1))
            seg_len = sqrt((ray%section(isec)%point(1,ipt+1) - ray%section(isec)%point(1,ipt))**2 + &
                           (ray%section(isec)%point(2,ipt+1) - ray%section(isec)%point(2,ipt))**2)
            if (seg_len < 1.0e-10_dp) cycle

            ivx = floor((px - vg%x(1)) / vg%dx) + 1
            if (ivx < 1) ivx = 1
            if (ivx >= nx_vg) ivx = nx_vg - 1
            if (pz < vg%z(1)) cycle
            if (pz >= vg%z(nz_vg)) then
              ivz = nz_vg - 1
              wz = 1.0_dp
            else
              ivz = floor((pz - vg%z(1)) / vg%dz) + 1
              if (ivz < 1) ivz = 1
              if (ivz >= nz_vg) ivz = nz_vg - 1
              wz = (pz - vg%z(ivz)) / vg%dz
            endif
            wx = (px - vg%x(ivx)) / vg%dx
            if (wx < 0.0_dp) wx = 0.0_dp
            if (wx > 1.0_dp) wx = 1.0_dp
            if (wz < 0.0_dp) wz = 0.0_dp
            if (wz > 1.0_dp) wz = 1.0_dp

            s_old = (1.0_dp-wx)*(1.0_dp-wz)/vg%velocity(vtype,ivx,ivz) + &
                    wx*(1.0_dp-wz)/vg%velocity(vtype,ivx+1,ivz) + &
                    (1.0_dp-wx)*wz/vg%velocity(vtype,ivx,ivz+1) + &
                    wx*wz/vg%velocity(vtype,ivx+1,ivz+1)

            if (ivx == ivx_test .and. ivz == iz_test) then
              s_new = (1.0_dp-wx)*(1.0_dp-wz)*(1.0_dp/vg%velocity(vtype,ivx,ivz)+ds) + &
                      wx*(1.0_dp-wz)/vg%velocity(vtype,ivx+1,ivz) + &
                      (1.0_dp-wx)*wz/vg%velocity(vtype,ivx,ivz+1) + &
                      wx*wz/vg%velocity(vtype,ivx+1,ivz+1)
            else if (ivx+1 == ivx_test .and. ivz == iz_test) then
              s_new = (1.0_dp-wx)*(1.0_dp-wz)/vg%velocity(vtype,ivx,ivz) + &
                      wx*(1.0_dp-wz)*(1.0_dp/vg%velocity(vtype,ivx+1,ivz)+ds) + &
                      (1.0_dp-wx)*wz/vg%velocity(vtype,ivx,ivz+1) + &
                      wx*wz/vg%velocity(vtype,ivx+1,ivz+1)
            else if (ivx == ivx_test .and. ivz+1 == iz_test) then
              s_new = (1.0_dp-wx)*(1.0_dp-wz)/vg%velocity(vtype,ivx,ivz) + &
                      wx*(1.0_dp-wz)/vg%velocity(vtype,ivx+1,ivz) + &
                      (1.0_dp-wx)*wz*(1.0_dp/vg%velocity(vtype,ivx,ivz+1)+ds) + &
                      wx*wz/vg%velocity(vtype,ivx+1,ivz+1)
            else if (ivx+1 == ivx_test .and. ivz+1 == iz_test) then
              s_new = (1.0_dp-wx)*(1.0_dp-wz)/vg%velocity(vtype,ivx,ivz) + &
                      wx*(1.0_dp-wz)/vg%velocity(vtype,ivx+1,ivz) + &
                      (1.0_dp-wx)*wz/vg%velocity(vtype,ivx,ivz+1) + &
                      wx*wz*(1.0_dp/vg%velocity(vtype,ivx+1,ivz+1)+ds)
            else
              s_new = s_old
            end if

            t_old = t_old + seg_len * s_old
            t_new = t_new + seg_len * s_new
          end do
        end do
      end do

      dt_fd = t_new - t_old
      dt_pred = ray%pdev(1) * ds

      if (n_test <= 3) then
        write(*,'(A,I4,A,I3,A,ES12.4,A,ES12.4,A,ES12.4,A,F10.4,A,I6,A,I4,A,I4)') &
          '  ircv=', ircv, ' iray=', iray, &
          ' pdev1=', ray%pdev(1), ' dt_fd=', dt_fd, ' dt_pred=', dt_pred, &
          ' ratio=', merge(dt_pred/dt_fd, 0.0_dp, abs(dt_fd)>1.0e-12_dp), &
          ' idx=', idx_test, ' ix=', ivx_test, ' iz=', iz_test
      end if

      if (abs(dt_fd) > 1.0e-12_dp) then
        n_match = n_match + 1
        sum_ratio = sum_ratio + dt_pred / dt_fd
      end if

      n_test = n_test + 1
      if (n_test <= 10) then
        write(*,'(A,I4,A,I3,A,ES12.4,A,ES12.4,A,F10.4,A,I6,A,I4,A,I4)') &
          '  ircv=', ircv, ' iray=', iray, &
          ' dt_fd=', dt_fd, ' dt_pred=', dt_pred, &
          ' ratio=', merge(dt_pred/dt_fd, 0.0_dp, abs(dt_fd)>1.0e-12_dp), &
          ' idx=', idx_test, ' ix=', ivx_test, ' iz=', iz_test
      end if

      if (n_test >= 20) exit
    end do
    if (n_test >= 20) exit
  end do

  if (n_match > 0) then
    write(*,'(A,F12.6,A,I4)') '  FD verify: avg(dt_pred/dt_fd)=', sum_ratio/n_match, ' n_match=', n_match
  end if

end subroutine fd_verify_frechet


subroutine verify_frechet_prediction(ndata, nparm, nnz, gindx, gval, &
  row_ptr, delta_m, mvscale)
  implicit none
  integer, intent(in) :: ndata, nparm, nnz
  integer, intent(in) :: gindx(nnz), row_ptr(ndata+1)
  real(kind=dp), intent(in) :: gval(nnz), delta_m(nparm), mvscale(nparm)

  integer :: i, k, j, ircv, iray, idata, n_check, n_valid_after
  real(kind=dp) :: step_test, dt_pred, sum_pred2, sum_act2, sum_ratio
  real(kind=dp), allocatable :: t_old(:), dt_pred_vec(:), dt_act_vec(:)
  real(kind=dp), allocatable :: data_wt_local(:)
  integer, allocatable :: old_rcv(:), old_iray(:)

  step_test = 0.01_dp

  allocate(t_old(ndata), old_rcv(ndata), old_iray(ndata))
  idata = 0
  do ircv = 1, n_receivers
    do iray = 1, receiver(ircv)%n_rays
      if (receiver(ircv)%ray(iray)%valid .and. &
          receiver(ircv)%ray(iray)%obs_time >= 0.0_dp) then
        idata = idata + 1
        if (idata > ndata) exit
        t_old(idata) = receiver(ircv)%ray(iray)%receiver_time
        old_rcv(idata) = ircv
        old_iray(idata) = iray
      end if
    end do
    if (idata > ndata) exit
  end do

  call save_model2d_old()
  call update_model2d(delta_m, step_test, mvscale)
  call run_forward2d()

  n_valid_after = 0
  do ircv = 1, n_receivers
    do iray = 1, receiver(ircv)%n_rays
      if (receiver(ircv)%ray(iray)%valid .and. &
          receiver(ircv)%ray(iray)%obs_time >= 0.0_dp) then
        n_valid_after = n_valid_after + 1
      end if
    end do
  end do

  if (n_valid_after /= ndata) then
    write(*,'(A,I6,A,I6,A)') '  WARNING: valid ray count changed from ', ndata, &
      ' to ', n_valid_after, ' during verify_frechet_prediction'
    write(*,'(A)') '  Skipping Frechet prediction verification due to ray count mismatch'
    deallocate(t_old, old_rcv, old_iray)
    call restore_model2d_old()
    return
  end if

  allocate(data_wt_local(ndata))
  do i = 1, ndata
    ircv = old_rcv(i)
    iray = old_iray(i)
    if (receiver(ircv)%ray(iray)%obs_time > 0.0_dp) then
      data_wt_local(i) = 1.0_dp
    else
      data_wt_local(i) = 1.0_dp
    end if
  end do

  allocate(dt_pred_vec(ndata), dt_act_vec(ndata))
  sum_pred2 = 0.0_dp
  sum_act2 = 0.0_dp
  sum_ratio = 0.0_dp
  n_check = 0

  do i = 1, ndata
    ircv = old_rcv(i)
    iray = old_iray(i)

    if (receiver(ircv)%ray(iray)%valid .and. &
        receiver(ircv)%ray(iray)%obs_time >= 0.0_dp) then
      dt_act_vec(i) = receiver(ircv)%ray(iray)%receiver_time - t_old(i)
    else
      dt_act_vec(i) = 0.0_dp
    end if

    dt_pred = 0.0_dp
    do k = row_ptr(i), row_ptr(i+1) - 1
      j = gindx(k)
      if (j >= 1 .and. j <= nparm) then
        dt_pred = dt_pred + gval(k) * delta_m(j)
      end if
    end do
    dt_pred = dt_pred * step_test / data_wt_local(i)
    dt_pred_vec(i) = dt_pred
    sum_pred2 = sum_pred2 + dt_pred**2
    sum_act2 = sum_act2 + dt_act_vec(i)**2
    if (abs(dt_act_vec(i)) > 1.0e-10_dp) then
      sum_ratio = sum_ratio + dt_pred / dt_act_vec(i)
      n_check = n_check + 1
    end if
  end do

  write(*,'(A)') '  Frechet prediction verification (step=0.01):'
  write(*,'(A,F12.8,A,F12.8)') '    RMS(dt_pred)=', sqrt(sum_pred2/ndata), &
    ' RMS(dt_act)=', sqrt(sum_act2/ndata)
  if (sum_act2 > 0.0_dp) then
    write(*,'(A,F12.6)') '    ratio RMS(dt_pred)/RMS(dt_act)=', &
      sqrt(sum_pred2/sum_act2)
  end if
  if (n_check > 0) then
    write(*,'(A,F12.6,A,I6)') '    avg(dt_pred/dt_act)=', sum_ratio/n_check, &
      ' n_check=', n_check
  end if

  do i = 1, min(ndata, 10)
    write(*,'(A,I4,A,ES12.4,A,ES12.4,A,F8.4)') '    ray=', i, &
      ' dt_pred=', dt_pred_vec(i), ' dt_act=', dt_act_vec(i), &
      ' ratio=', merge(dt_pred_vec(i)/dt_act_vec(i), 0.0_dp, abs(dt_act_vec(i))>1.0e-10_dp)
  end do

  deallocate(t_old, old_rcv, old_iray, data_wt_local, dt_pred_vec, dt_act_vec)
  call restore_model2d_old()

end subroutine verify_frechet_prediction


subroutine verify_frechet_fd(ndata, nparm, nnz, gindx, gval, row_ptr, residual)
  implicit none
  integer, intent(in) :: ndata, nparm, nnz
  integer, intent(in) :: gindx(nnz), row_ptr(ndata+1)
  real(kind=dp), intent(in) :: gval(nnz), residual(ndata)

  integer :: ig, ix, iz, idx, vtype, nx_vg, nz_vg, test_idx
  integer :: ircv, iray, k, n_test, n_match
  integer :: n_valid_before, n_valid_after, n_lost, n_new
  real(kind=dp) :: dv, t_old, t_new, dt_actual, dt_predicted, rel_err
  real(kind=dp) :: v_test
  real(kind=dp), allocatable :: old_times(:,:), old_pdev(:,:)
  logical, allocatable :: was_valid(:,:)
  type(Tvelocity_grid2d), pointer :: vg

  dv = 0.001_dp
  v_test = 0.0_dp

  test_idx = -1
  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    do iz = 20, nz_vg - 5
      do ix = nx_vg/3, 2*nx_vg/3
        test_idx = vg%inv_idx(ix, iz)
        if (test_idx > 0) then
          v_test = vg%velocity(vtype, ix, iz)
          exit
        end if
      end do
      if (test_idx > 0) exit
    end do
    if (test_idx > 0) exit
  end do

  if (test_idx <= 0 .or. test_idx > nparm) then
    write(*,'(A)') '  FD verify: could not find test node'
    return
  end if

  write(*,'(A,I6,A,F10.4,A,F10.4)') '  FD verify: test_idx=', test_idx, ' v_test=', v_test, ' dv=', dv

  allocate(old_times(n_receivers, 100), was_valid(n_receivers, 100))
  allocate(old_pdev(n_receivers, 100))
  old_times = 0.0_dp
  old_pdev = 0.0_dp
  was_valid = .false.
  n_valid_before = 0

  do ircv = 1, n_receivers
    do iray = 1, receiver(ircv)%n_rays
      if (receiver(ircv)%ray(iray)%valid .and. &
          receiver(ircv)%ray(iray)%obs_time >= 0.0_dp) then
        old_times(ircv, iray) = receiver(ircv)%ray(iray)%receiver_time
        was_valid(ircv, iray) = .true.
        n_valid_before = n_valid_before + 1
        do k = 1, receiver(ircv)%ray(iray)%n_pdev
          if (receiver(ircv)%ray(iray)%pdev_indx(k) == test_idx) then
            old_pdev(ircv, iray) = receiver(ircv)%ray(iray)%pdev(k)
            exit
          end if
        end do
      end if
    end do
  end do

  write(*,'(A,I6)') '  FD verify: n_valid_before=', n_valid_before

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0 .and. idx == test_idx) then
          vg%velocity(vtype, ix, iz) = vg%velocity(vtype, ix, iz) + dv
        end if
      end do
    end do
  end do

  call velocities_to_grid2d(pgrid)
  call run_forward2d()

  n_valid_after = 0
  n_lost = 0
  n_new = 0
  n_test = 0
  n_match = 0
  do ircv = 1, n_receivers
    do iray = 1, receiver(ircv)%n_rays
      if (receiver(ircv)%ray(iray)%valid .and. &
          receiver(ircv)%ray(iray)%obs_time >= 0.0_dp) then
        n_valid_after = n_valid_after + 1

        if (.not. was_valid(ircv, iray)) then
          n_new = n_new + 1
          cycle
        end if

        t_new = receiver(ircv)%ray(iray)%receiver_time
        t_old = old_times(ircv, iray)
        dt_actual = t_new - t_old

        dt_predicted = old_pdev(ircv, iray) * dv

        if (abs(old_pdev(ircv, iray)) > 1.0e-10_dp) then
          n_test = n_test + 1
          if (abs(dt_predicted) > 1.0e-10_dp) then
            rel_err = (dt_actual - dt_predicted) / dt_predicted
          else
            rel_err = 9999.0_dp
          end if
          if (n_test <= 10) then
            write(*,'(A,I5,A,I3,A,ES12.4,A,ES12.4,A,F8.2)') '      rcv=', ircv, ' ray=', iray, &
              ' dt_act=', dt_actual, ' dt_pred=', dt_predicted, ' rel=', rel_err * 100.0_dp
          end if
          if (abs(dt_predicted) > 1.0e-10_dp .and. abs(rel_err) < 0.1_dp) n_match = n_match + 1
        end if
      else
        if (was_valid(ircv, iray)) n_lost = n_lost + 1
      end if
    end do
  end do

  write(*,'(A,I6,A,I6,A,I3,A,I3)') '  FD verify: n_valid_before=', n_valid_before, &
    ' n_valid_after=', n_valid_after, ' n_lost=', n_lost, ' n_new=', n_new
  write(*,'(A,I5,A,I5,A)') '  FD verify: ', n_match, ' / ', n_test, &
    ' rays match within 10%'

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0 .and. idx == test_idx) then
          vg%velocity(vtype, ix, iz) = vg%velocity(vtype, ix, iz) - dv
        end if
      end do
    end do
  end do

  call run_forward2d()

  deallocate(old_times, was_valid, old_pdev)

end subroutine verify_frechet_fd


subroutine compute_path_integral_time2d(ray)
  implicit none
  type(Tray2d) :: ray

  integer :: isec, ipt, npts, ivx, ivz
  integer :: nx_vg, nz_vg, vtype
  real(kind=dp) :: px, pz, seg_len, vel
  real(kind=dp) :: wx, wz
  real(kind=dp) :: path_time
  type(Tvelocity_grid2d), pointer :: vg
  type(Tregion2d), pointer :: reg

  path_time = 0.0_dp
  vtype = VTYPE_P

  do isec = 1, ray%nsections
    if (ray%section(isec)%npoints <= 1) cycle
    reg => ray%section(isec)%reg

    vg => vgrid(reg%ivgrid, vtype)
    nx_vg = vg%nx
    nz_vg = vg%nz

    npts = ray%section(isec)%npoints

    do ipt = 1, npts - 1
      px = 0.5_dp * (ray%section(isec)%point(1,ipt) + ray%section(isec)%point(1,ipt+1))
      pz = 0.5_dp * (ray%section(isec)%point(2,ipt) + ray%section(isec)%point(2,ipt+1))

      seg_len = sqrt((ray%section(isec)%point(1,ipt+1) - ray%section(isec)%point(1,ipt))**2 + &
                     (ray%section(isec)%point(2,ipt+1) - ray%section(isec)%point(2,ipt))**2)

      if (seg_len < 1.0e-10_dp) cycle

      ivx = floor((px - vg%x(1)) / vg%dx) + 1
      if (ivx < 1) ivx = 1
      if (ivx >= nx_vg) ivx = nx_vg - 1

      if (pz < vg%z(1)) then
         cycle
      else if (pz >= vg%z(nz_vg)) then
         ivz = nz_vg - 1
         wz = 1.0_dp
      else
         ivz = floor((pz - vg%z(1)) / vg%dz) + 1
         if (ivz < 1) ivz = 1
         if (ivz >= nz_vg) ivz = nz_vg - 1
         wz = (pz - vg%z(ivz)) / vg%dz
      endif

      if (ivx < 1) ivx = 1
      if (ivx >= nx_vg) ivx = nx_vg - 1
      if (ivz < 1) ivz = 1
      if (ivz >= nz_vg) ivz = nz_vg - 1

      wx = (px - vg%x(ivx)) / vg%dx
      if (wx < 0.0_dp) wx = 0.0_dp
      if (wx > 1.0_dp) wx = 1.0_dp
      if (wz < 0.0_dp) wz = 0.0_dp
      if (wz > 1.0_dp) wz = 1.0_dp

      vel = (1.0_dp-wx)*(1.0_dp-wz)*vg%velocity(vtype,ivx,ivz) + &
            wx*(1.0_dp-wz)*vg%velocity(vtype,ivx+1,ivz) + &
            (1.0_dp-wx)*wz*vg%velocity(vtype,ivx,ivz+1) + &
            wx*wz*vg%velocity(vtype,ivx+1,ivz+1)

      if (vel < 1.0e-6_dp) vel = 1.0e-6_dp

      path_time = path_time + seg_len / vel
    end do
  end do

end subroutine compute_path_integral_time2d


subroutine rebuild_inv_idx_from_coverage()
  implicit none
  integer :: ig, ix, iz, idx, vtype, nx_vg, nz_vg, n_new, n_old
  integer :: ircv, iray, k
  logical, allocatable :: covered(:)
  type(Tvelocity_grid2d), pointer :: vg

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx
    nz_vg = vg%nz
    n_old = n_inv_parms

    allocate(covered(n_old))
    covered = .false.

    do ircv = 1, n_receivers
      do iray = 1, receiver(ircv)%n_rays
        if (.not. receiver(ircv)%ray(iray)%valid) cycle
        if (receiver(ircv)%ray(iray)%obs_time < 0.0_dp) cycle
        do k = 1, receiver(ircv)%ray(iray)%n_pdev
          idx = receiver(ircv)%ray(iray)%pdev_indx(k)
          if (idx >= 1 .and. idx <= n_old) then
            covered(idx) = .true.
          end if
        end do
      end do
    end do

    n_new = 0
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0 .and. idx <= n_old) then
          if (covered(idx)) then
            n_new = n_new + 1
            vg%inv_idx(ix, iz) = n_new
          else
            vg%inv_idx(ix, iz) = 0
          end if
        end if
      end do
    end do

    write(*,'(A,I8,A,I8,A,I8)') '    rebuild_inv_idx: n_covered=', count(covered), &
      ' n_old=', n_old, ' n_new=', n_new

    n_inv_parms = n_new
    n_inv_active = n_new

    deallocate(covered)
  end do

end subroutine rebuild_inv_idx_from_coverage


subroutine diagnose_normal_equations(ndata, nparm, nnz, gindx, gval, row_ptr, b)
  implicit none
  integer, intent(in) :: ndata, nparm, nnz
  integer, intent(in) :: gindx(nnz), row_ptr(ndata+1)
  real(kind=dp), intent(in) :: gval(nnz), b(ndata)

  real(kind=dp), allocatable :: gtg_diag(:), gtb(:), g_col_norm2(:)
  real(kind=dp) :: rsum, b_norm, gtb_norm, gtg_diag_norm
  integer :: i, j, k

  allocate(gtg_diag(nparm), gtb(nparm), g_col_norm2(nparm))
  gtg_diag = 0.0_dp
  gtb = 0.0_dp
  g_col_norm2 = 0.0_dp

  do i = 1, ndata
    do k = row_ptr(i), row_ptr(i+1) - 1
      j = gindx(k)
      if (j >= 1 .and. j <= nparm) then
        gtg_diag(j) = gtg_diag(j) + gval(k) * gval(k)
        gtb(j) = gtb(j) + gval(k) * b(i)
        g_col_norm2(j) = g_col_norm2(j) + gval(k) * gval(k)
      end if
    end do
  end do

  b_norm = sqrt(dot_product(b, b))
  gtb_norm = sqrt(dot_product(gtb, gtb))
  gtg_diag_norm = sqrt(dot_product(gtg_diag, gtg_diag))

  write(*,'(A,ES12.4)') '    |b|=', b_norm
  write(*,'(A,ES12.4)') '    |G^T*b|=', gtb_norm
  write(*,'(A,ES12.4)') '    |diag(G^T*G)|=', gtg_diag_norm
  write(*,'(A,ES12.4)') '    min(diag(G^T*G))=', minval(gtg_diag)
  write(*,'(A,ES12.4)') '    max(diag(G^T*G))=', maxval(gtg_diag)
  write(*,'(A,ES12.4)') '    median(diag(G^T*G))=', gtg_diag(nparm/2)

  call sort_array(gtg_diag, nparm)
  write(*,'(A,ES12.4)') '    10th pctile(diag(G^T*G))=', gtg_diag(max(1, nparm/10))
  write(*,'(A,ES12.4)') '    90th pctile(diag(G^T*G))=', gtg_diag(min(nparm, nparm*9/10))
  write(*,'(A,ES12.4)') '    condition(diag) =', &
    gtg_diag(min(nparm, nparm*9/10)) / max(gtg_diag(max(1, nparm/10)), 1.0e-30_dp)

  write(*,'(A,I8,A,I8)') '    n_covered=', count(gtg_diag > 0.0_dp), ' n_total=', nparm
  write(*,'(A,ES12.4)') '    min_nonzero(diag(G^T*G))=', minval(gtg_diag, mask=gtg_diag > 0.0_dp)

  deallocate(gtg_diag, gtb, g_col_norm2)

end subroutine diagnose_normal_equations

subroutine sort_array(arr, n)
  implicit none
  integer, intent(in) :: n
  real(kind=dp), intent(inout) :: arr(n)
  integer :: i, j
  real(kind=dp) :: tmp
  do i = 1, n-1
    do j = i+1, n
      if (arr(j) < arr(i)) then
        tmp = arr(i)
        arr(i) = arr(j)
        arr(j) = tmp
      end if
    end do
  end do
end subroutine sort_array


subroutine interpolate_iface_z_inv(iface, px, z_val)
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
end subroutine interpolate_iface_z_inv


subroutine find_optimal_entry_inv(src, ircv, ipath, itf)
  type(Tsource2d), target :: src
  integer :: ircv, ipath, itf
  integer :: i, nseafloor, nrefine, iref
  real(kind=dp) :: x_test, z_sf, t_fmm, t_water, t_total, t_best
  real(kind=dp) :: dist_water, x_best, z_best, wt_best
  real(kind=dp) :: x_left, x_right, dx_refine
  real(kind=dp) :: rx, rz, z_sf_init
  type(Ttime_field2d), pointer :: tf
  type(Tregion2d), pointer :: tf_reg

  rx = rec_x_orig(ircv)
  rz = rec_z_orig(ircv)

  tf => src%time_field(itf)
  if (.not. associated(tf%reg)) return

  tf_reg => tf%reg

  if (associated(tf%arrivaltime)) then
     if (.not. associated(tf_reg%loaded_tf, tf)) then
        tf_reg%arrivaltime = tf%arrivaltime
        tf_reg%time_gradient = tf%time_gradient
        tf_reg%loaded_tf => tf
     endif
  endif

  nseafloor = intrface(2)%nx

  t_best = huge_time
  call interpolate_iface_z_inv(intrface(2), rx, z_sf_init)
  x_best = rx
  z_best = z_sf_init
  wt_best = (z_sf_init - rz) / vwater_global

  do i = 1, nseafloor
     x_test = intrface(2)%x(i)
     z_sf = intrface(2)%z_nodes(i)

     call interpolate_arrivaltime2d(tf_reg, x_test, z_sf, t_fmm)
     if (t_fmm > huge_time * 0.5) cycle

     dist_water = sqrt((x_test - rx)**2 + (z_sf - rz)**2)
     t_water = dist_water / vwater_global
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
     call interpolate_iface_z_inv(intrface(2), x_test, z_sf)

     call interpolate_arrivaltime2d(tf_reg, x_test, z_sf, t_fmm)
     if (t_fmm > huge_time * 0.5) cycle

     dist_water = sqrt((x_test - rx)**2 + (z_sf - rz)**2)
     t_water = dist_water / vwater_global
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

end subroutine find_optimal_entry_inv


! ============================================================
! Read pre-computed entry points from forward model
! Format: ircv src_id path_id entry_x entry_z water_time
! ============================================================
subroutine read_entry_points_file(filename)
  character(len=*), intent(in) :: filename
  integer :: iunit, ios, npts, i
  integer :: ircv, src_id, path_id
  real(kind=dp) :: ex, ez, wt

  iunit = 99
  open(unit=iunit, file=trim(filename), status='old', action='read', iostat=ios)
  if (ios /= 0) then
     write(*,'(A,A)') 'WARNING: Cannot open entry points file: ', trim(filename)
     use_entry_points_file = .false.
     return
  endif

  read(iunit, *, iostat=ios) npts
  if (ios /= 0 .or. npts <= 0) then
     write(*,'(A)') 'WARNING: Invalid entry points file format'
     close(iunit)
     use_entry_points_file = .false.
     return
  endif

  n_entry_points = npts
  allocate(entry_point_data(6, npts))

  do i = 1, npts
     read(iunit, *, iostat=ios) ircv, src_id, path_id, ex, ez, wt
     if (ios /= 0) exit
     entry_point_data(1, i) = dble(ircv)
     entry_point_data(2, i) = dble(src_id)
     entry_point_data(3, i) = dble(path_id)
     entry_point_data(4, i) = ex
     entry_point_data(5, i) = ez
     entry_point_data(6, i) = wt
  end do

  close(iunit)
  use_entry_points_file = .true.
  write(*,'(A,I0,A)') '  Read ', npts, ' entry points from ', trim(filename)

end subroutine read_entry_points_file


! ============================================================
! Look up pre-computed entry point for a given receiver/source/raypath
! ============================================================
subroutine get_entry_point(ircv, src_id, path_id, entry_x, entry_z, water_t, found)
  integer, intent(in) :: ircv, src_id, path_id
  real(kind=dp), intent(out) :: entry_x, entry_z, water_t
  logical, intent(out) :: found
  integer :: i

  found = .false.
  entry_x = 0.0_dp
  entry_z = 0.0_dp
  water_t = 0.0_dp

  if (.not. use_entry_points_file .or. .not. allocated(entry_point_data)) return

  do i = 1, n_entry_points
     if (nint(entry_point_data(1,i)) == ircv .and. &
         nint(entry_point_data(2,i)) == src_id .and. &
         nint(entry_point_data(3,i)) == path_id) then
        entry_x = entry_point_data(4,i)
        entry_z = entry_point_data(5,i)
        water_t = entry_point_data(6,i)
        found = .true.
        return
     endif
  end do

end subroutine get_entry_point


! ============================================================
! Helper: compute avg_dmv for a given dm vector
! avg_dmv = sqrt(sum(dm^2) / n_active) where n_active = nodes with gtg > eps
! ============================================================
function calc_avg_dmv(dm, gtg_diag, eps_min, nparm, n_vel_parms) result(avg)
  implicit none
  integer, intent(in) :: nparm, n_vel_parms
  real(kind=dp), intent(in) :: dm(nparm), gtg_diag(nparm), eps_min
  real(kind=dp) :: avg
  integer :: j, n_act

  avg = 0.0_dp
  n_act = 0
  do j = 1, nparm
    if (j <= n_vel_parms) then
      if (gtg_diag(j) >= eps_min * 0.01_dp) then
        n_act = n_act + 1
        avg = avg + dm(j)**2
      end if
    end if
  end do
  if (n_act > 0) then
    avg = sqrt(avg / n_act)
  end if
end function calc_avg_dmv

! ============================================================
! Helper: compute avg_dmd for interface parameters
! avg_dmd = sqrt(sum(dz^2) / n_iface_parms) where dz = mvscale * dm
! ============================================================
function calc_avg_dmd(dm, mvscale_in, nparm, n_vel_parms) result(avg)
  implicit none
  integer, intent(in) :: nparm, n_vel_parms
  real(kind=dp), intent(in) :: dm(nparm), mvscale_in(nparm)
  real(kind=dp) :: avg
  integer :: j, n_act
  real(kind=dp) :: dz

  avg = 0.0_dp
  n_act = 0
  do j = n_vel_parms + 1, nparm
    n_act = n_act + 1
    dz = mvscale_in(j) * dm(j)
    avg = avg + dz**2
  end do
  if (n_act > 0) then
    avg = sqrt(avg / n_act)
  end if
end function calc_avg_dmd


! ============================================================
! Helper: run LSQR with given damping weight wdv
! Solves the augmented system (matching tomo2d):
!   [G_data        ] dm = [residual                ]
!   [w_smooth * R_h]      [-w_smooth*R_h*dm_sum    ]  (jumping)
!   [w_smooth * R_v]      [-w_smooth*R_v*dm_sum     ]  (jumping)
!   [wdv * I       ]      [-wdv*dm_sum              ]  (jumping)
! ============================================================
subroutine run_lsqr(wdv, wdv_iface, w_smooth, w_smooth_iface, ndata, nparm, nnz, gindx, gval, row_ptr, &
  residual, gtg_diag, eps_min, n_vel_parms, ndata_aug, max_lsqr_iter, &
  mvscale, u, v_vec, w_vec, z_vec, Av_aug, Atu_aug, Av_data, &
  smooth_work, smooth_out, dm_out, dm_sum, corr_x_in, corr_z_in, &
  smooth_type_in, lsqr_iter_out, rnorm_out, bnorm_out, test2_out)

  implicit none
  real(kind=dp), intent(in) :: wdv, wdv_iface, w_smooth, w_smooth_iface, eps_min
  integer, intent(in) :: ndata, nparm, nnz, n_vel_parms, ndata_aug, max_lsqr_iter
  integer, intent(in) :: smooth_type_in
  integer, intent(in) :: gindx(nnz), row_ptr(ndata+1)
  real(kind=dp), intent(in) :: gval(nnz), residual(ndata), gtg_diag(nparm)
  real(kind=dp), intent(in) :: mvscale(nparm)
  real(kind=dp), intent(in) :: dm_sum(nparm)
  real(kind=dp), intent(in) :: corr_x_in, corr_z_in
  real(kind=dp), intent(inout) :: u(ndata_aug), v_vec(nparm), w_vec(nparm), z_vec(nparm)
  real(kind=dp), intent(inout) :: Av_aug(ndata_aug), Atu_aug(nparm), Av_data(ndata)
  real(kind=dp), intent(inout) :: smooth_work(nparm), smooth_out(nparm)
  real(kind=dp), intent(out) :: dm_out(nparm)
  integer, intent(out) :: lsqr_iter_out
  real(kind=dp), intent(out) :: rnorm_out, bnorm_out, test2_out

  integer :: i, j, k, lsqr_iter
  real(kind=dp) :: beta_ls, alpha_ls, c_bar, s_bar, phi_bar, rho_bar
  real(kind=dp) :: rho_ls, theta_ls, c_ls, s_ls, phi_ls
  real(kind=dp) :: Anorm2_ls, Arnorm_ls, atol_ls
  real(kind=dp) :: smooth_h(nparm), smooth_v(nparm)
  integer :: t_ls0, t_ls1, t_ls_rate
  real(kind=dp) :: t_lsqr

  ! Layout depends on smooth_type:
  ! smooth_type=1 (Laplacian): [data(ndata) | L(nparm) | damp(nparm)]
  ! smooth_type=2 (Gaussian):  [data(ndata) | R_h(nparm) | R_v(nparm) | damp(nparm)]

  ! Initialize RHS with jumping strategy:
  u(1:ndata) = residual(1:ndata)

  if (smooth_type_in == 1) then
    ! ---- Laplacian smoothing ----
    if (w_smooth > 0.0_dp .or. w_smooth_iface > 0.0_dp) then
      call apply_laplacian_smooth_2d(nparm, dm_sum, smooth_out)
      do j = 1, nparm
        if (j <= n_vel_parms) then
          if (gtg_diag(j) < eps_min * 0.01_dp) then
            smooth_out(j) = 0.0_dp
          else
            smooth_out(j) = w_smooth * smooth_out(j)
          end if
        else
          smooth_out(j) = w_smooth_iface * smooth_out(j)
        end if
      end do
      u(ndata+1:ndata+nparm) = -smooth_out
    else
      u(ndata+1:ndata+nparm) = 0.0_dp
    end if

    ! Jumping: damping RHS = -wdv * dm_sum_vel - wdv_iface * dm_sum_iface
    if (wdv > 0.0_dp .or. wdv_iface > 0.0_dp) then
      do j = 1, nparm
        if (j <= n_vel_parms) then
          if (gtg_diag(j) < eps_min * 0.01_dp) then
            u(ndata+nparm+j) = 0.0_dp
          else
            u(ndata+nparm+j) = -wdv * dm_sum(j)
          end if
        else
          u(ndata+nparm+j) = -wdv_iface * dm_sum(j)
        end if
      end do
    else
      u(ndata+nparm+1:ndata+2*nparm) = 0.0_dp
    end if
  else if (smooth_type_in == 3) then
    ! ---- Hybrid smoothing: Laplacian for velocity, Gaussian R_h for interface ----
    if (w_smooth > 0.0_dp .or. w_smooth_iface > 0.0_dp) then
      ! Velocity: Laplacian
      call apply_laplacian_smooth_2d(nparm, dm_sum, smooth_out)
      do j = 1, n_vel_parms
        if (gtg_diag(j) < eps_min * 0.01_dp) then
          smooth_out(j) = 0.0_dp
        else
          smooth_out(j) = w_smooth * smooth_out(j)
        end if
      end do
      ! Interface: Gaussian R_h (overwrites interface part of smooth_out)
      if (w_smooth_iface > 0.0_dp) then
        call apply_iface_gauss_smooth_h(nparm, dm_sum, smooth_h)
        do j = n_vel_parms + 1, nparm
          smooth_out(j) = w_smooth_iface * smooth_h(j)
        end do
      else
        do j = n_vel_parms + 1, nparm
          smooth_out(j) = 0.0_dp
        end do
      end if
      u(ndata+1:ndata+nparm) = -smooth_out
    else
      u(ndata+1:ndata+nparm) = 0.0_dp
    end if

    ! Jumping: damping RHS = -wdv * dm_sum_vel - wdv_iface * dm_sum_iface
    if (wdv > 0.0_dp .or. wdv_iface > 0.0_dp) then
      do j = 1, nparm
        if (j <= n_vel_parms) then
          if (gtg_diag(j) < eps_min * 0.01_dp) then
            u(ndata+nparm+j) = 0.0_dp
          else
            u(ndata+nparm+j) = -wdv * dm_sum(j)
          end if
        else
          u(ndata+nparm+j) = -wdv_iface * dm_sum(j)
        end if
      end do
    else
      u(ndata+nparm+1:ndata+2*nparm) = 0.0_dp
    end if
  else
    ! ---- Gaussian smoothing (separable R_h, R_v) ----
    if (w_smooth > 0.0_dp) then
      call apply_gauss_smooth_h(nparm, dm_sum, smooth_h)
      call apply_gauss_smooth_v(nparm, dm_sum, smooth_v)
      do j = 1, nparm
        if (j <= n_vel_parms .and. gtg_diag(j) < eps_min * 0.01_dp) then
          smooth_h(j) = 0.0_dp
          smooth_v(j) = 0.0_dp
        end if
      end do
      u(ndata+1:ndata+nparm) = -w_smooth * smooth_h
      u(ndata+nparm+1:ndata+2*nparm) = -w_smooth * smooth_v
    else
      u(ndata+1:ndata+nparm) = 0.0_dp
      u(ndata+nparm+1:ndata+2*nparm) = 0.0_dp
    end if
    ! Interface Gaussian smoothing: placed in R_h block's interface rows
    if (w_smooth_iface > 0.0_dp .and. gsi_built) then
      call apply_iface_gauss_smooth_h(nparm, dm_sum, smooth_out)
      do j = n_vel_parms + 1, nparm
        u(ndata+j) = -w_smooth_iface * smooth_out(j)
      end do
    end if

    ! Jumping: damping RHS = -wdv * dm_sum_vel - wdv_iface * dm_sum_iface
    if (wdv > 0.0_dp .or. wdv_iface > 0.0_dp) then
      do j = 1, nparm
        if (j <= n_vel_parms) then
          if (gtg_diag(j) < eps_min * 0.01_dp) then
            u(ndata+2*nparm+j) = 0.0_dp
          else
            u(ndata+2*nparm+j) = -wdv * dm_sum(j)
          end if
        else
          u(ndata+2*nparm+j) = -wdv_iface * dm_sum(j)
        end if
      end do
    else
      u(ndata+2*nparm+1:ndata+3*nparm) = 0.0_dp
    end if
  end if
  bnorm_out = sqrt(dot_product(u, u))

  if (bnorm_out < 1.0e-30_dp) then
    dm_out = 0.0_dp
    lsqr_iter_out = 0
    rnorm_out = 0.0_dp
    test2_out = 0.0_dp
    return
  end if

  dm_out = 0.0_dp
  beta_ls = bnorm_out
  u = u / beta_ls

  ! v = A^T * u
  call sparse_ATx(ndata, nparm, nnz, gindx, gval, row_ptr, u(1:ndata), v_vec)
  if (smooth_type_in == 1) then
    ! Laplacian: L^T * u_L (L is symmetric)
    if (w_smooth > 0.0_dp .or. w_smooth_iface > 0.0_dp) then
      smooth_work = u(ndata+1:ndata+nparm)
      call apply_laplacian_smooth_2d(nparm, smooth_work, smooth_out)
      do j = 1, nparm
        if (j <= n_vel_parms) then
          if (gtg_diag(j) < eps_min * 0.01_dp) then
            smooth_out(j) = 0.0_dp
          else
            smooth_out(j) = w_smooth * smooth_out(j)
          end if
        else
          smooth_out(j) = w_smooth_iface * smooth_out(j)
        end if
      end do
      v_vec = v_vec + smooth_out
    end if
    if (wdv > 0.0_dp .or. wdv_iface > 0.0_dp) then
      do j = 1, nparm
        if (j <= n_vel_parms) then
          if (gtg_diag(j) >= eps_min * 0.01_dp) then
            v_vec(j) = v_vec(j) + wdv * u(ndata+nparm+j)
          end if
        else
          v_vec(j) = v_vec(j) + wdv_iface * u(ndata+nparm+j)
        end if
      end do
    end if
  else if (smooth_type_in == 3) then
    ! Hybrid: Laplacian^T for velocity, Gaussian R_h^T for interface
    if (w_smooth > 0.0_dp .or. w_smooth_iface > 0.0_dp) then
      smooth_work = u(ndata+1:ndata+nparm)
      ! Velocity: Laplacian^T (symmetric)
      call apply_laplacian_smooth_2d(nparm, smooth_work, smooth_out)
      do j = 1, n_vel_parms
        if (gtg_diag(j) < eps_min * 0.01_dp) then
          smooth_out(j) = 0.0_dp
        else
          smooth_out(j) = w_smooth * smooth_out(j)
        end if
      end do
      ! Interface: Gaussian R_h^T (symmetric)
      if (w_smooth_iface > 0.0_dp) then
        call apply_iface_gauss_smooth_h(nparm, smooth_work, smooth_h)
        do j = n_vel_parms + 1, nparm
          smooth_out(j) = w_smooth_iface * smooth_h(j)
        end do
      else
        do j = n_vel_parms + 1, nparm
          smooth_out(j) = 0.0_dp
        end do
      end if
      v_vec = v_vec + smooth_out
    end if
    if (wdv > 0.0_dp .or. wdv_iface > 0.0_dp) then
      do j = 1, nparm
        if (j <= n_vel_parms) then
          if (gtg_diag(j) >= eps_min * 0.01_dp) then
            v_vec(j) = v_vec(j) + wdv * u(ndata+nparm+j)
          end if
        else
          v_vec(j) = v_vec(j) + wdv_iface * u(ndata+nparm+j)
        end if
      end do
    end if
  else
    ! Gaussian: R_h^T * u_h + R_v^T * u_v
    if (w_smooth > 0.0_dp) then
      smooth_work = u(ndata+1:ndata+nparm)
      call apply_gauss_smooth_h_T(nparm, smooth_work, smooth_h)
      smooth_work = u(ndata+nparm+1:ndata+2*nparm)
      call apply_gauss_smooth_v_T(nparm, smooth_work, smooth_v)
      do j = 1, nparm
        if (j <= n_vel_parms .and. gtg_diag(j) < eps_min * 0.01_dp) then
          smooth_h(j) = 0.0_dp
          smooth_v(j) = 0.0_dp
        end if
      end do
      v_vec = v_vec + w_smooth * smooth_h + w_smooth * smooth_v
    end if
    ! Interface: R_if^T * u_h (symmetric: constant mvscale on interface nodes)
    if (w_smooth_iface > 0.0_dp .and. gsi_built) then
      smooth_work = u(ndata+1:ndata+nparm)
      call apply_iface_gauss_smooth_h(nparm, smooth_work, smooth_out)
      do j = n_vel_parms + 1, nparm
        v_vec(j) = v_vec(j) + w_smooth_iface * smooth_out(j)
      end do
    end if
    if (wdv > 0.0_dp .or. wdv_iface > 0.0_dp) then
      do j = 1, nparm
        if (j <= n_vel_parms) then
          if (gtg_diag(j) >= eps_min * 0.01_dp) then
            v_vec(j) = v_vec(j) + wdv * u(ndata+2*nparm+j)
          end if
        else
          v_vec(j) = v_vec(j) + wdv_iface * u(ndata+2*nparm+j)
        end if
      end do
    end if
  end if

  do j = 1, nparm
    if (j <= n_vel_parms .and. gtg_diag(j) < eps_min * 0.01_dp) then
      v_vec(j) = 0.0_dp
    end if
  end do

  alpha_ls = sqrt(dot_product(v_vec, v_vec))
  if (alpha_ls < 1.0e-30_dp) then
    dm_out = 0.0_dp
    lsqr_iter_out = 0
    rnorm_out = bnorm_out
    test2_out = 0.0_dp
    return
  end if

  v_vec = v_vec / alpha_ls
  w_vec = v_vec
  c_bar = 1.0_dp
  s_bar = 0.0_dp
  phi_bar = beta_ls
  rho_bar = alpha_ls
  Anorm2_ls = alpha_ls**2 + beta_ls**2
  ! Tolerance loosened from 1e-5: LSQR captures the dominant singular
  ! components within the first ~100 iterations, so a tighter tolerance only
  ! buys micro-accuracy at a large cost on high-ray-count datasets.
  atol_ls = 1.0e-4_dp

  call system_clock(t_ls0, t_ls_rate)
  do lsqr_iter = 1, max_lsqr_iter
    ! Av = A * v
    do j = 1, nparm
      if (j <= n_vel_parms .and. gtg_diag(j) < eps_min * 0.01_dp) then
        v_vec(j) = 0.0_dp
      end if
    end do
    call sparse_Ax(ndata, nparm, nnz, gindx, gval, row_ptr, v_vec, Av_data)
    Av_aug(1:ndata) = Av_data
    if (smooth_type_in == 1) then
      ! Laplacian: Av_smooth = w_vel * L_vel * v + w_iface * L_iface * v
      if (w_smooth > 0.0_dp .or. w_smooth_iface > 0.0_dp) then
        call apply_laplacian_smooth_2d(nparm, v_vec, smooth_out)
        do j = 1, nparm
          if (j <= n_vel_parms) then
            if (gtg_diag(j) < eps_min * 0.01_dp) then
              smooth_out(j) = 0.0_dp
            else
              smooth_out(j) = w_smooth * smooth_out(j)
            end if
          else
            smooth_out(j) = w_smooth_iface * smooth_out(j)
          end if
        end do
        Av_aug(ndata+1:ndata+nparm) = smooth_out
      else
        Av_aug(ndata+1:ndata+nparm) = 0.0_dp
      end if
      if (wdv > 0.0_dp .or. wdv_iface > 0.0_dp) then
        do j = 1, nparm
          if (j <= n_vel_parms) then
            if (gtg_diag(j) < eps_min * 0.01_dp) then
              Av_aug(ndata+nparm+j) = 0.0_dp
            else
              Av_aug(ndata+nparm+j) = wdv * v_vec(j)
            end if
          else
            Av_aug(ndata+nparm+j) = wdv_iface * v_vec(j)
          end if
        end do
      else
        Av_aug(ndata+nparm+1:ndata+2*nparm) = 0.0_dp
      end if
    else if (smooth_type_in == 3) then
      ! Hybrid: Av_smooth = w_vel * L_vel * v + w_iface * R_h_iface * v
      if (w_smooth > 0.0_dp .or. w_smooth_iface > 0.0_dp) then
        ! Velocity: Laplacian
        call apply_laplacian_smooth_2d(nparm, v_vec, smooth_out)
        do j = 1, n_vel_parms
          if (gtg_diag(j) < eps_min * 0.01_dp) then
            smooth_out(j) = 0.0_dp
          else
            smooth_out(j) = w_smooth * smooth_out(j)
          end if
        end do
        ! Interface: Gaussian R_h
        if (w_smooth_iface > 0.0_dp) then
          call apply_iface_gauss_smooth_h(nparm, v_vec, smooth_h)
          do j = n_vel_parms + 1, nparm
            smooth_out(j) = w_smooth_iface * smooth_h(j)
          end do
        else
          do j = n_vel_parms + 1, nparm
            smooth_out(j) = 0.0_dp
          end do
        end if
        Av_aug(ndata+1:ndata+nparm) = smooth_out
      else
        Av_aug(ndata+1:ndata+nparm) = 0.0_dp
      end if
      if (wdv > 0.0_dp .or. wdv_iface > 0.0_dp) then
        do j = 1, nparm
          if (j <= n_vel_parms) then
            if (gtg_diag(j) < eps_min * 0.01_dp) then
              Av_aug(ndata+nparm+j) = 0.0_dp
            else
              Av_aug(ndata+nparm+j) = wdv * v_vec(j)
            end if
          else
            Av_aug(ndata+nparm+j) = wdv_iface * v_vec(j)
          end if
        end do
      else
        Av_aug(ndata+nparm+1:ndata+2*nparm) = 0.0_dp
      end if
    else
      ! Gaussian: Av_smooth = w * R_h * v, w * R_v * v
      if (w_smooth > 0.0_dp) then
        call apply_gauss_smooth_h(nparm, v_vec, smooth_h)
        call apply_gauss_smooth_v(nparm, v_vec, smooth_v)
        do j = 1, nparm
          if (j <= n_vel_parms .and. gtg_diag(j) < eps_min * 0.01_dp) then
            smooth_h(j) = 0.0_dp
            smooth_v(j) = 0.0_dp
          end if
        end do
        Av_aug(ndata+1:ndata+nparm) = w_smooth * smooth_h
        Av_aug(ndata+nparm+1:ndata+2*nparm) = w_smooth * smooth_v
      else
        Av_aug(ndata+1:ndata+nparm) = 0.0_dp
        Av_aug(ndata+nparm+1:ndata+2*nparm) = 0.0_dp
      end if
      ! Interface Gaussian smoothing: placed in R_h block's interface rows
      if (w_smooth_iface > 0.0_dp .and. gsi_built) then
        call apply_iface_gauss_smooth_h(nparm, v_vec, smooth_h)
        do j = n_vel_parms + 1, nparm
          Av_aug(ndata+j) = w_smooth_iface * smooth_h(j)
        end do
      else
        do j = n_vel_parms + 1, nparm
          Av_aug(ndata+j) = 0.0_dp
        end do
      end if
      if (wdv > 0.0_dp .or. wdv_iface > 0.0_dp) then
        do j = 1, nparm
          if (j <= n_vel_parms) then
            if (gtg_diag(j) < eps_min * 0.01_dp) then
              Av_aug(ndata+2*nparm+j) = 0.0_dp
            else
              Av_aug(ndata+2*nparm+j) = wdv * v_vec(j)
            end if
          else
            Av_aug(ndata+2*nparm+j) = wdv_iface * v_vec(j)
          end if
        end do
      else
        Av_aug(ndata+2*nparm+1:ndata+3*nparm) = 0.0_dp
      end if
    end if

    u = Av_aug - alpha_ls * u
    beta_ls = sqrt(dot_product(u, u))
    if (beta_ls > 0.0_dp) u = u / beta_ls

    ! v = A^T * u
    call sparse_ATx(ndata, nparm, nnz, gindx, gval, row_ptr, u(1:ndata), Atu_aug)
    if (smooth_type_in == 1) then
      ! Laplacian: L^T * u_L
      if (w_smooth > 0.0_dp .or. w_smooth_iface > 0.0_dp) then
        smooth_work = u(ndata+1:ndata+nparm)
        call apply_laplacian_smooth_2d(nparm, smooth_work, smooth_out)
        do j = 1, nparm
          if (j <= n_vel_parms) then
            if (gtg_diag(j) < eps_min * 0.01_dp) then
              smooth_out(j) = 0.0_dp
            else
              smooth_out(j) = w_smooth * smooth_out(j)
            end if
          else
            smooth_out(j) = w_smooth_iface * smooth_out(j)
          end if
        end do
        Atu_aug = Atu_aug + smooth_out
      end if
      if (wdv > 0.0_dp .or. wdv_iface > 0.0_dp) then
        do j = 1, nparm
          if (j <= n_vel_parms) then
            if (gtg_diag(j) >= eps_min * 0.01_dp) then
              Atu_aug(j) = Atu_aug(j) + wdv * u(ndata+nparm+j)
            end if
          else
            Atu_aug(j) = Atu_aug(j) + wdv_iface * u(ndata+nparm+j)
          end if
        end do
      end if
    else if (smooth_type_in == 3) then
      ! Hybrid: Laplacian^T for velocity, Gaussian R_h^T for interface
      if (w_smooth > 0.0_dp .or. w_smooth_iface > 0.0_dp) then
        smooth_work = u(ndata+1:ndata+nparm)
        ! Velocity: Laplacian^T (symmetric)
        call apply_laplacian_smooth_2d(nparm, smooth_work, smooth_out)
        do j = 1, n_vel_parms
          if (gtg_diag(j) < eps_min * 0.01_dp) then
            smooth_out(j) = 0.0_dp
          else
            smooth_out(j) = w_smooth * smooth_out(j)
          end if
        end do
        ! Interface: Gaussian R_h^T (symmetric)
        if (w_smooth_iface > 0.0_dp) then
          call apply_iface_gauss_smooth_h(nparm, smooth_work, smooth_h)
          do j = n_vel_parms + 1, nparm
            smooth_out(j) = w_smooth_iface * smooth_h(j)
          end do
        else
          do j = n_vel_parms + 1, nparm
            smooth_out(j) = 0.0_dp
          end do
        end if
        Atu_aug = Atu_aug + smooth_out
      end if
      if (wdv > 0.0_dp .or. wdv_iface > 0.0_dp) then
        do j = 1, nparm
          if (j <= n_vel_parms) then
            if (gtg_diag(j) >= eps_min * 0.01_dp) then
              Atu_aug(j) = Atu_aug(j) + wdv * u(ndata+nparm+j)
            end if
          else
            Atu_aug(j) = Atu_aug(j) + wdv_iface * u(ndata+nparm+j)
          end if
        end do
      end if
    else
      ! Gaussian: R_h^T * u_h + R_v^T * u_v
      if (w_smooth > 0.0_dp) then
        smooth_work = u(ndata+1:ndata+nparm)
        call apply_gauss_smooth_h_T(nparm, smooth_work, smooth_h)
        smooth_work = u(ndata+nparm+1:ndata+2*nparm)
        call apply_gauss_smooth_v_T(nparm, smooth_work, smooth_v)
        do j = 1, nparm
          if (j <= n_vel_parms .and. gtg_diag(j) < eps_min * 0.01_dp) then
            smooth_h(j) = 0.0_dp
            smooth_v(j) = 0.0_dp
          end if
        end do
        Atu_aug = Atu_aug + w_smooth * smooth_h + w_smooth * smooth_v
      end if
      ! Interface: R_if^T * u_h (symmetric: constant mvscale on interface nodes)
      if (w_smooth_iface > 0.0_dp .and. gsi_built) then
        smooth_work = u(ndata+1:ndata+nparm)
        call apply_iface_gauss_smooth_h(nparm, smooth_work, smooth_h)
        do j = n_vel_parms + 1, nparm
          Atu_aug(j) = Atu_aug(j) + w_smooth_iface * smooth_h(j)
        end do
      end if
      if (wdv > 0.0_dp .or. wdv_iface > 0.0_dp) then
        do j = 1, nparm
          if (j <= n_vel_parms) then
            if (gtg_diag(j) >= eps_min * 0.01_dp) then
              Atu_aug(j) = Atu_aug(j) + wdv * u(ndata+2*nparm+j)
            end if
          else
            Atu_aug(j) = Atu_aug(j) + wdv_iface * u(ndata+2*nparm+j)
          end if
        end do
      end if
    end if

    v_vec = Atu_aug - beta_ls * v_vec
    do j = 1, nparm
      if (j <= n_vel_parms .and. gtg_diag(j) < eps_min * 0.01_dp) then
        v_vec(j) = 0.0_dp
      end if
    end do
    alpha_ls = sqrt(dot_product(v_vec, v_vec))
    if (alpha_ls > 0.0_dp) v_vec = v_vec / alpha_ls

    Anorm2_ls = Anorm2_ls + alpha_ls**2 + beta_ls**2
    rho_ls = sqrt(rho_bar**2 + beta_ls**2)
    c_ls = rho_bar / rho_ls
    s_ls = beta_ls / rho_ls
    theta_ls = s_ls * alpha_ls
    rho_bar = -c_ls * alpha_ls
    phi_ls = c_ls * phi_bar
    phi_bar = s_ls * phi_bar

    z_vec = w_vec
    do j = 1, nparm
      if (j <= n_vel_parms .and. gtg_diag(j) < eps_min * 0.01_dp) then
        z_vec(j) = 0.0_dp
      end if
    end do
    w_vec = v_vec - (theta_ls / rho_ls) * w_vec
    do j = 1, nparm
      if (j <= n_vel_parms .and. gtg_diag(j) < eps_min * 0.01_dp) then
        w_vec(j) = 0.0_dp
      end if
    end do
    dm_out = dm_out + (phi_ls / rho_ls) * z_vec

    rnorm_out = abs(phi_bar)
    Arnorm_ls = phi_bar * alpha_ls * abs(c_ls)
    test2_out = Arnorm_ls / (sqrt(Anorm2_ls) * rnorm_out + 1.0e-30_dp)

    if (test2_out < atol_ls) exit
    if (1.0_dp + test2_out <= 1.0_dp) exit
  end do

  lsqr_iter_out = lsqr_iter
  call system_clock(t_ls1)
  t_lsqr = real(t_ls1 - t_ls0, kind=dp) / real(t_ls_rate, kind=dp)
  write(*,'(A,I5,A,ES10.3,A,ES10.3,A,F9.2,A)') '    LSQR converged: iter=', lsqr_iter, &
    ' rnorm=', rnorm_out, ' test2=', test2_out, ' time=', t_lsqr, ' s'

end subroutine run_lsqr


!--------------------------------------------------------------
! Apply model update (velocity + interface) with given step sizes
!--------------------------------------------------------------
subroutine apply_model_step(step_len, step_iface, dmodel_total, mvscale, nparm, n_vel_parms, &
    saved_iface_z, n_iface_parms)
  implicit none
  real(kind=dp), intent(in) :: step_len, step_iface
  real(kind=dp), intent(in) :: dmodel_total(:), mvscale(:)
  integer, intent(in) :: nparm, n_vel_parms
  real(kind=dp), intent(in) :: saved_iface_z(:)
  integer, intent(in) :: n_iface_parms

  integer :: ig, ix, iz, idx, vtype, nx_vg, nz_vg
  real(kind=dp) :: cur_s_v, trial_s, trial_v, init_v
  type(Tvelocity_grid2d), pointer :: vg
  integer :: ip, ii, iiface
  real(kind=dp) :: dz
  type(Tinterface2d), pointer :: ifc

  ! Update velocity
  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx; nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0) then
          if (step_len <= 0.0_dp) then
            cycle  ! no velocity update
          end if
          cur_s_v = 1.0_dp / vg%velocity(vtype, ix, iz)
          trial_s = cur_s_v + step_len * mvscale(idx) * dmodel_total(idx)
          trial_v = 2.0_dp
          if (trial_s > 0.0_dp) trial_v = 1.0_dp / trial_s
          if (trial_v < 2.0_dp) trial_v = 2.0_dp
          if (trial_v > 10.0_dp) trial_v = 10.0_dp
          init_v = vg%velocity_init(vtype, ix, iz)
          if (init_v >= 2.0_dp) then
            if (trial_v > init_v * 1.30_dp) trial_v = init_v * 1.30_dp
            if (trial_v < init_v * 0.70_dp) trial_v = init_v * 0.70_dp
          end if
          vg%velocity(vtype, ix, iz) = trial_v
        end if
      end do
    end do
  end do

  ! Update interface depths
  if (n_inv_iface > 0 .and. step_iface > 0.0_dp) then
    ip = 0
    do ii = 1, n_inv_iface
      iiface = ifaces_to_be_inv(ii)
      ifc => intrface(iiface)
      do ix = 1, ifc%nx
        ip = ip + 1
        idx = ifc%start_index + ix
        if (idx >= 1 .and. idx <= nparm) then
          dz = step_iface * mvscale(idx) * dmodel_total(idx)
          if (dz > 0.5_dp) dz = 0.5_dp
          if (dz < -0.5_dp) dz = -0.5_dp
          ifc%z_nodes(ix) = ifc%z_nodes(ix) + dz
        end if
      end do
    end do
    call rebuild_geometry2d()
  end if
end subroutine apply_model_step

!--------------------------------------------------------------
! Restore model state to saved values
!--------------------------------------------------------------
subroutine restore_model_state(saved_vel, saved_iface_z, n_iface_parms)
  implicit none
  real(kind=dp), intent(in) :: saved_vel(:,:)
  real(kind=dp), intent(in) :: saved_iface_z(:)
  integer, intent(in) :: n_iface_parms

  integer :: ig, ix, iz, vtype, nx_vg, nz_vg
  type(Tvelocity_grid2d), pointer :: vg
  integer :: ip, ii, iiface
  type(Tinterface2d), pointer :: ifc

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx; nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        vg%velocity(vtype, ix, iz) = saved_vel(ix, iz)
      end do
    end do
  end do

  if (n_inv_iface > 0) then
    ip = 0
    do ii = 1, n_inv_iface
      iiface = ifaces_to_be_inv(ii)
      ifc => intrface(iiface)
      do ix = 1, ifc%nx
        ip = ip + 1
        ifc%z_nodes(ix) = saved_iface_z(ip)
      end do
    end do
    call rebuild_geometry2d()
  end if
end subroutine restore_model_state

! Helper subroutine to apply model update with a given step length (velocity only)
subroutine apply_model_update(saved_vel, step_len, mvscale, dmodel_total)
  implicit none
  real(kind=dp), intent(in) :: saved_vel(:,:)
  real(kind=dp), intent(in) :: step_len
  real(kind=dp), intent(in) :: mvscale(:)
  real(kind=dp), intent(in) :: dmodel_total(:)

  integer :: ig, ix, iz, idx, vtype, nx_vg, nz_vg
  real(kind=dp) :: cur_s_v, trial_s, trial_v, init_v
  type(Tvelocity_grid2d), pointer :: vg

  do ig = 1, n_inv_vgrid
    vg => vgrid(vgrids_to_be_inv(ig)%igrid, vgrids_to_be_inv(ig)%vtype)
    vtype = vgrids_to_be_inv(ig)%vtype
    nx_vg = vg%nx; nz_vg = vg%nz
    do iz = 1, nz_vg
      do ix = 1, nx_vg
        idx = vg%inv_idx(ix, iz)
        if (idx > 0) then
          if (step_len <= 0.0_dp) then
            ! Revert to original model
            vg%velocity(vtype, ix, iz) = saved_vel(ix, iz)
          else
            cur_s_v = 1.0_dp / saved_vel(ix, iz)
            trial_s = cur_s_v + step_len * mvscale(idx) * dmodel_total(idx)
            trial_v = 2.0_dp
            if (trial_s > 0.0_dp) trial_v = 1.0_dp / trial_s
            if (trial_v < 2.0_dp) trial_v = 2.0_dp
            if (trial_v > 10.0_dp) trial_v = 10.0_dp
            init_v = vg%velocity_init(vtype, ix, iz)
            if (init_v >= 2.0_dp) then
              if (trial_v > init_v * 1.10_dp) trial_v = init_v * 1.10_dp
              if (trial_v < init_v * 0.90_dp) trial_v = init_v * 0.90_dp
            end if
            vg%velocity(vtype, ix, iz) = trial_v
          end if
        end if
      end do
    end do
  end do
end subroutine apply_model_update


end module mod_inversion2d
