module type_definitions2d

  INTEGER, PARAMETER     :: sp = SELECTED_REAL_KIND(6,37)
  INTEGER, PARAMETER     :: dp = SELECTED_REAL_KIND(15,307)

  INTEGER, PARAMETER     :: VTYPE_P = 1
  INTEGER, PARAMETER     :: VTYPE_S = 2
  INTEGER, PARAMETER     :: NVTYPE  = 2

  INTEGER, PARAMETER     :: NODE_FAR   = -1
  INTEGER, PARAMETER     :: NODE_TRIAL =  0
  INTEGER, PARAMETER     :: NODE_ALIVE =  1

! declarations of derived types
!------------------------------------------------------

type Tpointer_to_integer_array

   integer,dimension(:),pointer   :: p

end type Tpointer_to_integer_array


type Tpropagation_grid2d

    integer       :: nx,nz
    REAL(KIND=dp) :: dx,dz
    REAL(KIND=dp) :: x0,z0
    REAL(KIND=dp) :: xmax,zmax
    REAL(KIND=dp) :: tolerance
    logical       :: is_main_grid
    logical       :: is_source_grid

    integer       :: nnode

    integer       :: index_x0,index_z0

    REAL(KIND=dp), DIMENSION (:), pointer :: x
    REAL(KIND=dp), DIMENSION (:), pointer :: z

    REAL(KIND=dp), DIMENSION (:,:,:), pointer     :: velocity
    REAL(KIND=dp), DIMENSION (:,:),   pointer     :: arrivaltime
    REAL(KIND=dp), DIMENSION (:,:,:), pointer     :: time_gradient

    integer, DIMENSION (:,:), pointer :: node_region
    integer, DIMENSION (:,:), pointer :: rnode_id

    logical, dimension(:,:),  pointer :: fully_regular

    type(Tpointer_to_integer_array), DIMENSION (:,:), pointer :: ccind_from_2dc

end type Tpropagation_grid2d

!-----------------------------------------------------------------------------------

type Tvelocity_grid2d

    integer       :: nx,nz
    REAL(KIND=dp) :: dx,dz
    REAL(KIND=dp) :: x0,z0

    integer       :: nnode

    integer       :: start_index
    logical       :: to_be_inverted

    REAL(KIND=dp), DIMENSION (:), pointer :: x
    REAL(KIND=dp), DIMENSION (:), pointer :: z

    REAL(KIND=dp), DIMENSION (:,:,:), pointer     :: velocity
    REAL(KIND=dp), DIMENSION (:,:,:), pointer     :: velocity_init
    logical,dimension(:,:,:),pointer              :: active
    logical,dimension(:,:,:),pointer              :: active_s  ! S-wave active mask (inactive in water layer)
    integer, dimension(:,:), pointer :: inv_idx
    logical                                       :: vs_available  ! .true. if Vs data is loaded

end type Tvelocity_grid2d


!-------------------------------

type Tpath2d

    integer                            :: id
    integer                            :: n_tf
    integer,dimension(:) ,pointer      :: sequence
    integer,dimension(:) ,pointer      :: tf_sequence
    integer,dimension(:) ,pointer      :: vtype_sequence

    logical                            :: valid
    logical                            :: used
    logical                            :: gridsave
    integer                            :: first_tf_to_save

    integer                            :: refstep
    integer                            :: fitting_interface
    integer                            :: refraction  ! 0=reflection (default), 1=refraction/diving wave

end type Tpath2d

!-------------------------------------------------------------------------------------------------------


type Tbackpointer2d

    integer :: i1
    integer :: i2
    integer :: i3

end type Tbackpointer2d

!-------------------------------

type Tsource2d

    integer                  :: id
    REAL(KIND=dp)            :: x,z
    integer                  :: ix,iz
    type(Tbackpointer2d)     :: cnode(100)
    integer                  :: n_cnode

    logical                  :: on_grid,on_interface,on_pinched_interface
    integer                  :: region_id,interface_id
    integer                  :: topreg_id,botreg_id,topint_id,botint_id

    integer                  :: n_tf_init

    integer,dimension(:),pointer   :: first_tf_up
    integer,dimension(:),pointer   :: first_tf_down

    integer                              :: n_paths
    type(Tpath2d),dimension(:),pointer   :: path

    integer                                 :: n_time_fields
    type(Ttime_field2d),dimension(:),pointer  :: time_field

    integer                  :: start_index
    logical                  :: to_be_inverted

end type Tsource2d


!-------------------------------!

type Treceiver2d

    integer                                    :: id
    REAL(KIND=dp)                              :: x,z
    integer                                    :: n_rays
    type(Tray2d),dimension(:),pointer          :: ray

end type Treceiver2d

!-------------------------------!

type Tinterface2d

    integer :: nx,id
    REAL(KIND=dp) :: dx
    REAL(KIND=dp) :: x0
    logical       :: pinched

    integer       :: nnode
    integer       :: start_index
    logical       :: to_be_inverted

    REAL(KIND=dp), DIMENSION (:), pointer :: x
    REAL(KIND=dp), DIMENSION (:), pointer :: z_nodes
    REAL(KIND=dp), DIMENSION (:), pointer :: z_nodes_init

end type Tinterface2d

!-------------------------------!

type Tinteger_coordinates2d

    integer :: ix,iz

end type Tinteger_coordinates2d

!-------------------------------

type Tintersection2d

   integer :: nnode,n_ccells
   integer :: id
   integer :: iface_id
   type(Tpropagation_grid2d),pointer :: grid

   logical :: pinched

   REAL(KIND=dp), DIMENSION (:), pointer              :: x
   REAL(KIND=dp), DIMENSION (:), pointer              :: z
   REAL(KIND=dp), DIMENSION (:,:), pointer            :: normal

   REAL(KIND=dp), DIMENSION (:,:), pointer              :: vel_top
   REAL(KIND=dp), DIMENSION (:,:), pointer              :: vel_bot

   REAL(KIND=dp), DIMENSION (:), pointer              :: arrivaltime
   REAL(KIND=dp), DIMENSION (:), pointer              :: starttime
   REAL(KIND=dp), DIMENSION (:,:), pointer            :: time_gradient

   integer,dimension(:),pointer           :: intype

   type(Tinteger_coordinates2d), DIMENSION (:), pointer :: ccells
   integer,dimension(:), pointer                      :: n_inodes
   integer,dimension(:,:),pointer                     :: inodes
   integer, DIMENSION (:,:), pointer                  :: ccell_from_inode

   type(Tregion2d),pointer                             :: regleft
   type(Tregion2d),pointer                             :: regright

   integer, DIMENSION (:), pointer          :: rleft_node_id
   integer, DIMENSION (:), pointer          :: rright_node_id

   integer, DIMENSION (:), pointer          :: irg_left
   integer, DIMENSION (:), pointer          :: irg_right

end type Tintersection2d

!-------------------------------

type Tregion2d

   integer      :: id
   integer      :: ivgrid
   type(Tintersection2d),pointer :: itop,ibot
   type(Tpropagation_grid2d),pointer :: grid
   integer      :: ngnode

   integer                                     :: nnode
   type(Tbackpointer2d),dimension(:),pointer   :: node
   integer, dimension(:), pointer              :: node_status
   real(kind=dp),dimension(:),pointer          :: arrivaltime
   real(kind=dp),dimension(:,:),pointer        :: time_gradient
   ! Time field whose data currently occupies arrivaltime / time_gradient.
   ! A region is shared by several time fields, so its time arrays are used as
   ! a scratch view of the currently active time field. Recording which time
   ! field is loaded lets callers skip a re-copy that would be a no-op.
   ! Null means the content is FMM scratch and must not be trusted.
   type(Ttime_field2d), pointer                :: loaded_tf => null()
   real(kind=dp),dimension(:,:),pointer        :: velocity
   real(kind=dp),dimension(:),pointer          :: x
   real(kind=dp),dimension(:),pointer          :: z

   integer, dimension(:,:), pointer            :: region_rnode_id

   ! Per-region boundary node IDs (safe from overwrite by later define_region2d calls)
   integer, dimension(:), pointer              :: itop_rnode_id => null()
   integer, dimension(:), pointer              :: ibot_rnode_id => null()

   integer                                     :: n_init
   integer,dimension(:),pointer                :: init_id
   integer,dimension(:),pointer                :: init_type
   real(kind=dp),dimension(:),pointer          :: init_arrivaltime
   real(kind=dp),dimension(:,:),pointer        :: init_time_gradient

end type Tregion2d

!-------------------------------

type Ttime_field2d

   integer                                    :: id
   integer                                    :: vtype
   integer                                    :: nfile
   type(Tregion2d), pointer                   :: reg
   type(Tintersection2d),pointer              :: istart
   type(Tintersection2d),pointer              :: inonstart
   real(kind=dp),dimension(:),pointer         :: arrivaltime
   real(kind=dp),dimension(:,:),pointer       :: time_gradient

   logical                                    :: turning_rays_present
   logical,dimension(:),pointer               :: received_turning_ray

   integer                                    :: prev_tf
   integer                                    :: next_tf(4)

end type Ttime_field2d

!-------------------------------

type Tray_section2d

   type(Tray2d),pointer                       :: ray
   type(Tregion2d), pointer                   :: reg
   type(Tintersection2d),pointer              :: istart
   type(Tintersection2d),pointer              :: iend
   type(Ttime_field2d),pointer                :: tf
   type(Tsource2d),pointer                    :: source

   integer                                    :: npoints
   real(kind=dp),dimension(:,:),pointer       :: point

   integer                                    :: place_in_sequence
   logical                                    :: diffracted
   logical                                    :: headwave

end type Tray_section2d

!-------------------------------

type Tray2d

   integer                                      :: nsections
   type(Tray_section2d),dimension(:),pointer    :: section
   type(Tsource2d),pointer                      :: source
   integer                                      :: raypath_id
   integer                                      :: source_id
   real(kind=dp)                                :: receiver_time
   real(kind=dp),dimension(2)                   :: receiver_time_gradient
   logical                                      :: diffracted
   logical                                      :: headwave
   logical                                      :: direct_water
   logical                                      :: is_multiray
   integer                                      :: n_subrays
   type(Tray2d),dimension(:),pointer            :: subray

   integer                                      :: n_pdev
   integer,dimension(:),pointer                 :: pdev_indx
   real(kind=dp),dimension(:),pointer           :: pdev

   logical                                      :: valid

   real(kind=dp)                                :: obs_time
   real(kind=dp)                                :: water_time   ! per-ray water layer correction

end type Tray2d


type Tgrid_identifier2d

   integer   :: igrid,vtype

end type Tgrid_identifier2d


!**************************************************************************************************************************

contains

  subroutine pgrid2d_defaults(grid)

    type(Tpropagation_grid2d)  :: grid

    grid%is_main_grid = .false.
    grid%is_source_grid  = .false.
    grid%nnode = 0

    nullify(grid%x)
    nullify(grid%z)
    nullify(grid%velocity)
    nullify(grid%arrivaltime)
    nullify(grid%time_gradient)
    nullify(grid%node_region)
    nullify(grid%rnode_id)
    nullify(grid%fully_regular)
    nullify(grid%ccind_from_2dc)

  end subroutine pgrid2d_defaults


  subroutine vgrid2d_defaults(grid)

    type(Tvelocity_grid2d)  :: grid

    grid%to_be_inverted  = .false.
    grid%nnode = 0
    grid%vs_available = .false.
    nullify(grid%x)
    nullify(grid%z)
    nullify(grid%velocity)
    nullify(grid%velocity_init)
    nullify(grid%active)
    nullify(grid%active_s)
    nullify(grid%inv_idx)

  end subroutine vgrid2d_defaults


  subroutine source2d_defaults(source)
    type(Tsource2d)  :: source

    source%on_grid = .false.
    source%on_interface = .false.
    source%on_pinched_interface = .false.
    source%region_id = 0
    source%interface_id = 0
    source%topreg_id = 0
    source%botreg_id = 0
    source%topint_id = 0
    source%botint_id = 0
    source%n_tf_init = 1
    nullify(source%first_tf_up)
    nullify(source%first_tf_down)
    source%n_paths = 0
    source%to_be_inverted  = .false.
    nullify(source%path)
    source%n_time_fields = 0
    nullify(source%time_field)

  end subroutine source2d_defaults


  subroutine path2d_defaults(path)
    type(Tpath2d)  :: path

    path%id = 0
    path%n_tf = 0
    nullify(path%sequence)
    nullify(path%tf_sequence)
    nullify(path%vtype_sequence)
    path%valid = .true.
    path%used  = .false.
    path%gridsave  = .false.
    path%refstep = 0
    path%fitting_interface = 0
    path%refraction = 0  ! 0=reflection (default)

  end subroutine path2d_defaults


  subroutine receiver2d_defaults(rec)
    type(Treceiver2d) :: rec

    rec%id  = 0
    rec%n_rays  = 0
    nullify(rec%ray)

  end subroutine receiver2d_defaults


  subroutine intersection2d_defaults(isec)
    type(Tintersection2d)  :: isec
    isec%nnode = 0
    isec%n_ccells = 0
    isec%id = 0
    isec%iface_id = 0
    nullify(isec%grid)
    isec%pinched = .false.
    nullify(isec%x)
    nullify(isec%z)
    nullify(isec%normal)
    nullify(isec%vel_top)
    nullify(isec%vel_bot)
    nullify(isec%arrivaltime)
    nullify(isec%starttime)
    nullify(isec%time_gradient)
    nullify(isec%intype)
    nullify(isec%ccells)
    nullify(isec%n_inodes)
    nullify(isec%inodes)
    nullify(isec%ccell_from_inode)
    nullify(isec%regleft)
    nullify(isec%regright)
    nullify(isec%rleft_node_id)
    nullify(isec%rright_node_id)
    nullify(isec%irg_left)
    nullify(isec%irg_right)

  end subroutine intersection2d_defaults


  subroutine interface2d_defaults(iface)
    type(Tinterface2d)  :: iface

    iface%pinched = .false.
    iface%nnode = 0
    iface%to_be_inverted = .false.
    nullify(iface%x)
    nullify(iface%z_nodes)

  end subroutine interface2d_defaults


  subroutine region2d_defaults(reg)
    type(Tregion2d)  :: reg

    reg%id = 0
    reg%ivgrid = 0
    nullify(reg%grid)
    reg%ngnode   = 0
    reg%nnode  = 0
    reg%n_init = 0
    nullify(reg%node)
    nullify(reg%node_status)
    nullify(reg%arrivaltime)
    nullify(reg%time_gradient)
    nullify(reg%loaded_tf)
    nullify(reg%velocity)
    nullify(reg%x)
    nullify(reg%z)
    nullify(reg%init_id)
   nullify(reg%init_type)
   nullify(reg%init_arrivaltime)
   nullify(reg%init_time_gradient)

  end subroutine region2d_defaults


  subroutine time_field2d_defaults(tf)
    type(Ttime_field2d)  :: tf

    tf%id = 0
    tf%vtype = 0
    nullify(tf%reg)
    nullify(tf%istart)
    nullify(tf%inonstart)
    nullify(tf%arrivaltime)
    nullify(tf%time_gradient)
    tf%turning_rays_present = .false.
    nullify(tf%received_turning_ray)
    tf%prev_tf = 0
    tf%next_tf(1:4) = 0

  end subroutine time_field2d_defaults


  subroutine ray2d_defaults(ray)
    type(Tray2d)  :: ray

    ray%nsections = 0
    nullify(ray%section)
    nullify(ray%source)
    ray%raypath_id = 0
    ray%source_id = 0
    ray%diffracted = .false.
    ray%headwave = .false.
    ray%direct_water = .false.
    ray%is_multiray = .false.
    ray%n_subrays = 0
    nullify(ray%subray)
    ray%n_pdev = 0
    nullify(ray%pdev_indx)
    nullify(ray%pdev)
    ray%valid = .true.
    ray%obs_time = -1.0_dp
    ray%water_time = 0.0_dp

  end subroutine ray2d_defaults


  subroutine ray_section2d_defaults(raysec)
    type(Tray_section2d) :: raysec

    nullify(raysec%ray)
    nullify(raysec%reg)
    nullify(raysec%istart)
    nullify(raysec%iend)
    nullify(raysec%tf)
    nullify(raysec%source)
    raysec%npoints = 0
    nullify(raysec%point)
    raysec%place_in_sequence = 0
    raysec%diffracted = .false.
    raysec%headwave = .false.

  end subroutine ray_section2d_defaults


  !-----------------------------------------------------------------------
  ! Get velocity for a given vtype (VTYPE_P or VTYPE_S) at grid point (ix,iz)
  ! For S-wave in water layer (where active_s is .false.), returns 0.0
  !-----------------------------------------------------------------------
  function get_velocity2d(vgrid_in, ix, iz, vtype) result(vel)
    type(Tvelocity_grid2d), intent(in) :: vgrid_in
    integer, intent(in)                :: ix, iz, vtype
    real(kind=dp)                      :: vel

    vel = 0.0_dp

    if (vtype == VTYPE_S) then
      ! S-wave: check active_s mask (inactive in water layer)
      if (associated(vgrid_in%active_s)) then
        if (.not. vgrid_in%active_s(vtype, ix, iz)) return
      end if
    end if

    if (associated(vgrid_in%velocity)) then
      vel = vgrid_in%velocity(vtype, ix, iz)
    end if

  end function get_velocity2d

  !-----------------------------------------------------------------------
  ! Check if S-wave is active at grid point (ix,iz)
  ! Returns .false. in water layer (above seafloor)
  !-----------------------------------------------------------------------
  function is_swave_active2d(vgrid_in, ix, iz) result(is_active)
    type(Tvelocity_grid2d), intent(in) :: vgrid_in
    integer, intent(in)                :: ix, iz
    logical                            :: is_active

    is_active = .true.

    if (associated(vgrid_in%active_s)) then
      is_active = vgrid_in%active_s(VTYPE_S, ix, iz)
    end if

  end function is_swave_active2d


end module type_definitions2d

!*************************************************************************************************************************


module global_variables2d

  use type_definitions2d

! global parameters

  REAL(KIND=dp), PARAMETER   :: interface_tolerance = 0.005_dp

  REAL(KIND=dp), PARAMETER   :: huge_time = 1.0e20_dp

  integer   :: refinement_factor
  integer   :: ncell_to_be_refined
  integer   :: global_source_counter
  integer   :: raypoint_counter

  logical   :: file_mode
  logical   :: no_pp_mode
  logical   :: parallel_mode
  logical   :: display_mode
  logical   :: save_rays_mode
  logical   :: save_timefields_mode

!.........................................................................
  integer                                                  :: n_interfaces
  type(Tinterface2d),dimension(:),pointer                  :: intrface


!.........................................................................
  integer                                                  :: n_vgrids
  integer                                                  :: n_vtypes = 1  ! set to 2 when Vs is available
  type(Tvelocity_grid2d),dimension(:,:),pointer            :: vgrid



!.........................................................................
  type(Tpropagation_grid2d),pointer                        :: pgrid

  integer                                                  :: n_intersections
  type(Tintersection2d),dimension(:),pointer               :: intersection

  integer                                                  :: n_regions
  type(Tregion2d),dimension(:),pointer                     :: region


!.........................................................................
  type(Tpropagation_grid2d),pointer                        :: sgrid

  integer                                                  :: n_sintersections
  type(Tintersection2d),dimension(:),pointer               :: sintersection

  integer                                                  :: n_sregions
  type(Tregion2d),dimension(:),pointer                     :: sregion


!.........................................................................
  integer                                                  :: n_receivers
  type(Treceiver2d),dimension(:),pointer                   :: receiver

  integer                                                  :: n_sources
  integer                                                  :: n_sources_ppinc
  type(Tsource2d),dimension(:),pointer                     :: source

!.........................................................................
  integer                                                  :: n_inv_parms
  integer                                                  :: n_inv_active
  integer                                                  :: n_inv_vgrid
  type(Tgrid_identifier2d) ,dimension(:),pointer           :: vgrids_to_be_inv
  integer                                                  :: n_inv_iface
  integer,dimension(:),pointer                             :: ifaces_to_be_inv

  logical                                                  :: locate_source
  integer                                                  :: n_inv_source
  integer,dimension(:),pointer                             :: sources_to_be_inv

  real(kind=dp)                                            :: vwater_global
  real(kind=dp)                                            :: vair_global

  ! Original grid info (before resampling to pgrid)
  integer                                                  :: nz_orig
  real(kind=dp), allocatable                               :: z_orig(:)

  ! Smoothing type: 1=Laplacian, 2=Gaussian (separable R_h/R_v)
  integer                                                  :: smooth_type = 2

  logical, allocatable                                     :: in_water(:)
  real(kind=dp), allocatable                               :: rec_x_orig(:)
  real(kind=dp), allocatable                               :: rec_z_orig(:)
  real(kind=dp), allocatable                               :: water_time(:)
  real(kind=dp), allocatable                               :: entry_x_best(:)
  real(kind=dp), allocatable                               :: entry_z_best(:)

end module global_variables2d

!*************************************************************************************************************


module interface_definitions2d

interface
   subroutine propagate2d(regin,vtype)
     use type_definitions2d
     type(Tregion2d),target    :: regin
     integer                   :: vtype
   end subroutine
end interface

interface
   subroutine find_intersection2d(isec,iface,grid)
     use type_definitions2d
     type(Tintersection2d)               :: isec
     type(Tinterface2d)                  :: iface
     type(Tpropagation_grid2d),target    :: grid
   end subroutine
end interface

interface
   subroutine define_region2d(reg,itop,ibot,grid)
     use type_definitions2d
     type(Tregion2d),target     :: reg
     type(Tintersection2d)      :: itop,ibot
     type(Tpropagation_grid2d),target :: grid
   end subroutine
end interface

interface
   subroutine sweep_region_from_interface2d(reg,istart_in,vtype,s)
     use type_definitions2d
     type(Tregion2d)                          :: reg
     type(Tintersection2d),target             :: istart_in
     integer                                  :: vtype
     type(Tsource2d)                          :: s
   end subroutine
end interface

interface
   subroutine transfer_diving_return_to_intersection(isec,tf,vtype,min_depth)
     use type_definitions2d
     type(Tintersection2d) :: isec
     type(Ttime_field2d) :: tf
     integer :: vtype
     real(kind=dp) :: min_depth
   end subroutine
end interface

interface
   subroutine initialize_region_from_source2d(reg,sx,sz,t0,vtype_in)
     use type_definitions2d
     type(Tregion2d),target :: reg
     real(kind=dp) :: sx,sz,t0
     integer, optional :: vtype_in
   end subroutine
end interface

interface
   subroutine sweep_sregion_from_interface2d(reg,istart_in,vtype)
     use type_definitions2d
     type(Tregion2d)                          :: reg
     type(Tintersection2d),target             :: istart_in
     integer                                  :: vtype
   end subroutine
end interface

interface
   subroutine initialize_refined_source2d(s,sc,grid,reg,itop,ibot)
     use type_definitions2d
     type(Tsource2d) :: s
     type(Tsource2d) :: sc
     type(Tpropagation_grid2d) :: grid
     type(Tregion2d)         :: reg
     type(Tintersection2d),target   :: itop,ibot
   end subroutine
end interface

end module interface_definitions2d

!**********************************************************************************************************

module mod_fmt2d

  use type_definitions2d
  use global_variables2d
  use interface_definitions2d

end module mod_fmt2d

module mod_fmt2d_nointerfaces

  use type_definitions2d
  use global_variables2d

end module mod_fmt2d_nointerfaces
