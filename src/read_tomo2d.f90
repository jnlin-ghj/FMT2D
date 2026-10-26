module read_tomo2d
  implicit none
  private
  public :: read_smesh, read_refl, read_geom

contains

  subroutine read_smesh(filename, nx, nz, vwater, xpos, topo, zpos_rel, vel, vair)
    character(len=*), intent(in)               :: filename
    integer, intent(out)                        :: nx, nz
    real(8), intent(out)                        :: vwater
    real(8), allocatable, intent(out)           :: xpos(:)
    real(8), allocatable, intent(out)           :: topo(:)
    real(8), allocatable, intent(out)           :: zpos_rel(:)
    real(8), allocatable, intent(out)           :: vel(:,:)
    real(8), optional, intent(out)              :: vair

    integer :: iunit, ios, ix, iz
    real(8) :: vair_in
    character(len=512) :: line_buf

    iunit = 10
    open(unit=iunit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') 'Error: cannot open smesh file: ', trim(filename)
      stop 1
    end if

    read(iunit, *, iostat=ios) nx, nz, vwater, vair_in
    if (ios /= 0) then
      write(*,'(A)') 'Error: reading smesh header (nx, nz, vwater, vair)'
      stop 1
    end if
    if (present(vair)) vair = vair_in

    ! Validate grid dimensions
    if (nx < 2 .or. nz < 2) then
      write(*,'(A,I0,A,I0)') 'Error: smesh grid too small: nx=', nx, ' nz=', nz
      write(*,'(A)') '       require at least nx>=2 and nz>=2'
      stop 1
    end if

    allocate(xpos(nx), topo(nx), zpos_rel(nz), vel(nx, nz))

    read(iunit, *, iostat=ios) xpos(1:nx)
    if (ios /= 0) then
      write(*,'(A)') 'Error: reading xpos from smesh'
      stop 1
    end if

    read(iunit, *, iostat=ios) topo(1:nx)
    if (ios /= 0) then
      write(*,'(A)') 'Error: reading topo from smesh'
      stop 1
    end if

    read(iunit, *, iostat=ios) zpos_rel(1:nz)
    if (ios /= 0) then
      write(*,'(A)') 'Error: reading zpos_rel from smesh'
      stop 1
    end if

    do ix = 1, nx
      read(iunit, *, iostat=ios) (vel(ix, iz), iz = 1, nz)
      if (ios /= 0) then
        write(*,'(A,I0)') 'Error: reading vel row for ix = ', ix
        stop 1
      end if
    end do

    close(iunit)

    write(*,'(A,A)')        '  Read smesh: ', trim(filename)
    write(*,'(A,I0,A,I0)')  '    nx = ', nx, '  nz = ', nz
    write(*,'(A,F8.3)')     '    vwater = ', vwater
    write(*,'(A,F8.3,A,F8.3)') '    x range: ', xpos(1), ' to ', xpos(nx)
    write(*,'(A,F8.3,A,F8.3)') '    topo range: ', minval(topo), ' to ', maxval(topo)
    write(*,'(A,F8.3,A,F8.3)') '    zpos_rel range: ', zpos_rel(1), ' to ', zpos_rel(nz)
    write(*,'(A,F8.3,A,F8.3)') '    vel range: ', minval(vel), ' to ', maxval(vel)

  end subroutine read_smesh

  subroutine read_refl(filename, nx_refl, refl_x, refl_z)
    character(len=*), intent(in)               :: filename
    integer, intent(out)                        :: nx_refl
    real(8), allocatable, intent(out)           :: refl_x(:)
    real(8), allocatable, intent(out)           :: refl_z(:)

    integer :: iunit, ios, nlines, max_lines
    real(8), allocatable :: tmp_x(:), tmp_z(:)
    character(len=512) :: line_buf

    max_lines = 100000

    iunit = 11
    open(unit=iunit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') 'Error: cannot open refl file: ', trim(filename)
      stop 1
    end if

    allocate(tmp_x(max_lines), tmp_z(max_lines))

    nlines = 0
    do
      read(iunit, '(A)', iostat=ios) line_buf
      if (ios /= 0) exit
      if (len_trim(line_buf) == 0) cycle
      nlines = nlines + 1
      if (nlines > max_lines) then
        write(*,'(A)') 'Error: refl file exceeds max_lines'
        stop 1
      end if
      read(line_buf, *, iostat=ios) tmp_x(nlines), tmp_z(nlines)
      if (ios /= 0) then
        write(*,'(A,I0)') 'Warning: skipping malformed line at refl line ', nlines
        nlines = nlines - 1
        cycle
      end if
    end do

    close(iunit)

    nx_refl = nlines
    if (nx_refl < 1) then
      write(*,'(A,A)') 'Error: refl file has no valid data: ', trim(filename)
      stop 1
    end if
    allocate(refl_x(nx_refl), refl_z(nx_refl))
    refl_x(1:nx_refl) = tmp_x(1:nx_refl)
    refl_z(1:nx_refl) = tmp_z(1:nx_refl)
    deallocate(tmp_x, tmp_z)

    write(*,'(A,A)')       '  Read refl: ', trim(filename)
    write(*,'(A,I0)')      '    nx_refl = ', nx_refl
    write(*,'(A,F8.3,A,F8.3)') '    x range: ', refl_x(1), ' to ', refl_x(nx_refl)
    write(*,'(A,F8.3,A,F8.3)') '    z range: ', minval(refl_z), ' to ', maxval(refl_z)

  end subroutine read_refl

  subroutine read_geom(filename, nobs, obs_x, obs_z, nrec_total, rec_x, rec_z, obs_id, rec_phase, obs_time, uncertainty)
    character(len=*), intent(in)               :: filename
    integer, intent(out)                        :: nobs
    real(8), allocatable, intent(out)           :: obs_x(:)
    real(8), allocatable, intent(out)           :: obs_z(:)
    integer, intent(out)                        :: nrec_total
    real(8), allocatable, intent(out)           :: rec_x(:)
    real(8), allocatable, intent(out)           :: rec_z(:)
    integer, allocatable, intent(out)           :: obs_id(:)
    integer, allocatable, intent(out)           :: rec_phase(:)
    real(8), allocatable, intent(out)           :: obs_time(:)
    real(8), allocatable, intent(out)           :: uncertainty(:)

    integer :: iunit, ios, iobs, irec, nrec_i, idx
    character(len=1) :: tag
    real(8) :: rx, rz, t_obs, unc
    integer :: ray_type
    integer, allocatable :: nrec_per_obs(:)
    integer :: max_rec

    iunit = 12
    open(unit=iunit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') 'Error: cannot open geom file: ', trim(filename)
      stop 1
    end if

    read(iunit, *, iostat=ios) nobs
    if (ios /= 0) then
      write(*,'(A)') 'Error: reading nobs from geom file'
      stop 1
    end if

    if (nobs < 1) then
      write(*,'(A,I0)') 'Error: invalid nobs from geom file: ', nobs
      stop 1
    end if

    allocate(obs_x(nobs), obs_z(nobs), nrec_per_obs(nobs))

    max_rec = 0
    do iobs = 1, nobs
      read(iunit, *, iostat=ios) tag, obs_x(iobs), obs_z(iobs), nrec_per_obs(iobs)
      if (ios /= 0) then
        write(*,'(A,I0)') 'Error: reading OBS header for iobs = ', iobs
        stop 1
      end if
      max_rec = max_rec + nrec_per_obs(iobs)
      do irec = 1, nrec_per_obs(iobs)
        read(iunit, *, iostat=ios) tag, rx, rz, ray_type, t_obs, unc
        if (ios /= 0) then
          write(*,'(A,I0,A,I0)') 'Error: reading rec for iobs = ', iobs, ' irec = ', irec
          stop 1
        end if
      end do
    end do

    close(iunit)

    nrec_total = max_rec
    allocate(rec_x(nrec_total), rec_z(nrec_total), obs_id(nrec_total), rec_phase(nrec_total))
    allocate(obs_time(nrec_total), uncertainty(nrec_total))

    open(unit=iunit, file=trim(filename), status='old', action='read', iostat=ios)
    if (ios /= 0) then
      write(*,'(A,A)') 'Error: cannot reopen geom file: ', trim(filename)
      stop 1
    end if

    read(iunit, *, iostat=ios) nobs
  if (ios /= 0) then
    write(*,'(A)') 'Error: re-reading nobs from geom file'
    stop 1
  end if

  idx = 0
  do iobs = 1, nobs
    read(iunit, *, iostat=ios) tag, obs_x(iobs), obs_z(iobs), nrec_per_obs(iobs)
    if (ios /= 0) then
      write(*,'(A,I0)') 'Error: re-reading OBS header for iobs = ', iobs
      stop 1
    end if
    do irec = 1, nrec_per_obs(iobs)
      idx = idx + 1
      if (idx > nrec_total) then
        write(*,'(A)') 'Error: geom file inconsistent between passes (idx > nrec_total)'
        stop 1
      end if
      read(iunit, *, iostat=ios) tag, rec_x(idx), rec_z(idx), ray_type, t_obs, unc
      if (ios /= 0) then
        write(*,'(A,I0,A,I0)') 'Error: re-reading rec for iobs = ', iobs, ' irec = ', irec
        stop 1
      end if
      obs_id(idx) = iobs
      rec_phase(idx) = ray_type
      obs_time(idx) = t_obs
      uncertainty(idx) = unc
    end do
  end do

    close(iunit)
    deallocate(nrec_per_obs)

    write(*,'(A,A)')       '  Read geom: ', trim(filename)
    write(*,'(A,I0)')      '    nobs = ', nobs
    write(*,'(A,I0)')      '    nrec_total = ', nrec_total

  end subroutine read_geom

end module read_tomo2d