module mod_prop2d
use mod_fmt2d_nointerfaces

type(Tregion2d),pointer :: reg
type(Tpropagation_grid2d),pointer :: grid
real(kind=dp),dimension(:),pointer :: velocity
type(Tvelocity_grid2d),pointer :: velgrid
integer :: current_vtype

integer :: fallback_count, fallback_nodes

integer, dimension(:),pointer :: node_from_tree_ind
integer, dimension(:),pointer :: node_is_counted
logical,dimension(:),pointer :: suspect_time_gradient

type(Tinteger_coordinates2d),dimension(4) :: concell
integer,dimension(20) :: connode
integer,dimension(4) :: concell_nnode
integer,dimension(4) :: concell_nbase
logical,dimension(4) :: concell_is_regular
integer,dimension(20) :: trialnode

integer :: ntree, newalive, n_connode, n_concell, n_trial
integer :: n_det_negative
real(kind=dp) :: oldtime
integer,dimension(:),pointer :: in
logical :: diag = .false.
logical :: diat = .false.

end module mod_prop2d


subroutine propagate2d(regin,vtype)

use mod_prop2d
implicit none

type(Tregion2d),target :: regin
integer :: vtype
integer :: n,nt
logical :: first_step

reg => regin
grid => reg%grid
velgrid => vgrid(reg%ivgrid,vtype)
current_vtype = vtype

allocate(velocity(reg%nnode))
velocity = reg%velocity(vtype,:)

n_det_negative = 0
fallback_count = 0
fallback_nodes = 0

allocate(suspect_time_gradient(reg%nnode))
suspect_time_gradient = .false.

allocate(node_from_tree_ind(reg%nnode),node_is_counted(reg%nnode))
node_is_counted = 0
ntree = 0

do n=1,reg%nnode
   if (reg%node_status(n) > 0) then
      if (current_vtype == VTYPE_S .and. velocity(n) < 0.01_dp) then
         reg%node_status(n) = -1
         reg%arrivaltime(n) = huge_time
      else
         call add_node_to_tree(n)
      endif
   endif
end do

first_step = .true.
do while (ntree > 0)

   newalive = node_from_tree_ind(1)
   reg%node_status(newalive) = 0
   call remove_root_from_tree

   if (ntree == 0 .and. .not.first_step) exit

   first_step = .false.

   call find_connected_nodes2d(newalive)

   n_trial = n_connode
   trialnode(1:n_trial) = connode(1:n_trial)

   do nt=1,n_trial

      if (reg%node_status(trialnode(nt)) > 0) then

         oldtime = reg%arrivaltime(trialnode(nt))
         call update_time2d(trialnode(nt))

         if (reg%arrivaltime(trialnode(nt)) < oldtime) call update_tree(trialnode(nt))

      endif

      if (reg%node_status(trialnode(nt)) == -1) then

         if (current_vtype == VTYPE_S .and. velocity(trialnode(nt)) < 0.01_dp) cycle

         reg%arrivaltime(trialnode(nt)) = huge_time
         oldtime = reg%arrivaltime(trialnode(nt))
         call update_time2d(trialnode(nt))

         if (reg%arrivaltime(trialnode(nt)) < oldtime) call add_node_to_tree(trialnode(nt))

      endif
   end do

end do

deallocate(node_from_tree_ind,node_is_counted)

if (n_det_negative > 0) print *,'warning: determinant negative ',n_det_negative,' times'

do n=1,reg%nnode
   if (suspect_time_gradient(n)) then
      call fit_gradient_at_node2d(reg,n,reg%time_gradient(1,n),reg%time_gradient(2,n))
   endif
end do

deallocate(suspect_time_gradient,velocity)

end subroutine propagate2d


subroutine find_connected_nodes2d(centernode)
use mod_prop2d
implicit none

integer :: centernode
integer :: i1,i2,i3
integer :: n,m,i,j,ii,jj,icell,ix,ixx,iz_s,izz,n_before
type(Tintersection2d),pointer :: isec

i1 = reg%node(centernode)%i1
i2 = reg%node(centernode)%i2
i3 = reg%node(centernode)%i3

if (i1 /= 0) then
   if (grid%fully_regular(i1,i2)) then

      n_connode = 0
   if (reg%region_rnode_id(i1-1,i2) > 0) then
      n_connode = n_connode + 1
      connode(n_connode) = reg%region_rnode_id(i1-1,i2)
   endif
   if (reg%region_rnode_id(i1+1,i2) > 0) then
      n_connode = n_connode + 1
      connode(n_connode) = reg%region_rnode_id(i1+1,i2)
   endif
   if (reg%region_rnode_id(i1,i2-1) > 0) then
      n_connode = n_connode + 1
      connode(n_connode) = reg%region_rnode_id(i1,i2-1)
   endif
   if (reg%region_rnode_id(i1,i2+1) > 0) then
      n_connode = n_connode + 1
      connode(n_connode) = reg%region_rnode_id(i1,i2+1)
   endif

      return

   endif
endif

if (i1 /= 0) then

   n_concell = 0
   do i=0,1
      if ((i1-i > 0) .and. (i1-i < grid%nx)) then
         do j=0,1
            if ((i2-j > 0) .and. (i2-j < grid%nz)) then
               n_concell = n_concell + 1
               concell(n_concell)%ix = i1 - i
               concell(n_concell)%iz = i2 - j
            endif
         end do
      end if
   end do

else

   if (i2 == reg%itop%id) then
      isec => reg%itop
   else
      isec => reg%ibot
   endif
   n_concell = 0
   do n=1,4
      icell = isec%ccell_from_inode(n,i3)
      if (icell > 0) then
         n_concell = n_concell + 1
         concell(n_concell)%ix = isec%ccells(icell)%ix
         concell(n_concell)%iz = isec%ccells(icell)%iz
      endif
   enddo

endif

node_is_counted(centernode) = centernode

n_connode = 0

do n=1,n_concell

   concell_is_regular(n) = (i1 /= 0)

   if (associated(grid%ccind_from_2dc(concell(n)%ix,concell(n)%iz)%p)) then

      icell = grid%ccind_from_2dc(concell(n)%ix,concell(n)%iz)%p(reg%itop%iface_id)
      if (icell /= 0) then

         concell_is_regular(n) = .false.

         isec => reg%itop
         do jj=1,isec%n_inodes(icell)

            m = reg%itop_rnode_id(isec%inodes(jj,icell))

            if (node_is_counted(m) /= centernode) then
               n_connode = n_connode + 1
               connode(n_connode) = m
               node_is_counted(m) = centernode
            endif

         end do

      end if

      icell = grid%ccind_from_2dc(concell(n)%ix,concell(n)%iz)%p(reg%ibot%iface_id)
      if (icell /= 0) then

         concell_is_regular(n) = .false.

         isec => reg%ibot
         do jj=1,isec%n_inodes(icell)

            m = reg%ibot_rnode_id(isec%inodes(jj,icell))

            if (node_is_counted(m) /= centernode) then
               n_connode = n_connode + 1
               connode(n_connode) = m
               node_is_counted(m) = centernode
            endif

         end do

      endif

   endif

   do i=0,1
      ii = concell(n)%ix + i
      do j=0,1
         jj = concell(n)%iz + j

         if (ii >= 1 .and. ii <= grid%nx .and. jj >= 1 .and. jj <= grid%nz) then
            if ((i1 /= 0 .and. abs(i1-ii)+abs(i2-jj) == 1) .or. (.not.concell_is_regular(n))) then

               m = reg%region_rnode_id(ii,jj)
               if (m > 0) then
                  if (node_is_counted(m) /= centernode) then
                     n_connode = n_connode + 1
                     connode(n_connode) = m
                     node_is_counted(m) = centernode
                  endif
               endif

            endif
         endif

      end do
   end do

end do

! Fallback for interface nodes: directly search for grid nodes
! adjacent to the interface that belong to the current region.
! This is needed because ccell_from_inode may point to cells whose
! corner nodes belong to a different region (overlapping z-boundaries).
! Always do this for interface nodes - node_is_counted prevents duplicates.
if (i1 == 0) then

   fallback_count = fallback_count + 1
   n_before = n_connode

   ! Get the grid position of this interface node
   ix = nint((reg%x(centernode) - grid%x0) / grid%dx) + 1
   if (ix < 1) ix = 1
   if (ix > grid%nx) ix = grid%nx

   ! Get the grid z-index closest to the interface
   iz_s = nint((reg%z(centernode) - grid%z0) / grid%dz) + 1
   if (iz_s < 1) iz_s = 1
   if (iz_s > grid%nz) iz_s = grid%nz

   ! Determine if this is a top or bottom interface node
   if (i2 == reg%ibot%id) then
      ! Bottom interface: region is above, search upward (decreasing iz)
      do ixx = max(1, ix-1), min(grid%nx, ix+1)
         do izz = iz_s, 1, -1
            m = reg%region_rnode_id(ixx, izz)
            if (m > 0) then
               if (node_is_counted(m) /= centernode) then
                  n_connode = n_connode + 1
                  connode(n_connode) = m
                  node_is_counted(m) = centernode
               endif
               exit  ! found the nearest node above the interface
            endif
         end do
      end do
   else
      ! Top interface: region is below, search downward (increasing iz)
      do ixx = max(1, ix-1), min(grid%nx, ix+1)
         do izz = iz_s, grid%nz
            m = reg%region_rnode_id(ixx, izz)
            if (m > 0) then
               if (node_is_counted(m) /= centernode) then
                  n_connode = n_connode + 1
                  connode(n_connode) = m
                  node_is_counted(m) = centernode
               endif
               exit  ! found the nearest node below the interface
            endif
         end do
      end do
   endif

   fallback_nodes = fallback_nodes + (n_connode - n_before)

endif

return

end subroutine find_connected_nodes2d


subroutine find_connected_cells2d(centernode)
use mod_prop2d
implicit none

integer :: centernode
integer :: i1,i2,i3
integer :: n,m,i,j,ii,jj,icell
type(Tintersection2d),pointer :: isec

i1 = reg%node(centernode)%i1
i2 = reg%node(centernode)%i2
i3 = reg%node(centernode)%i3

if (i1 /= 0) then

   n_concell = 0
   do i=0,1
      if ((i1-i > 0) .and. (i1-i < grid%nx)) then
         do j=0,1
            if ((i2-j > 0) .and. (i2-j < grid%nz)) then
               n_concell = n_concell + 1
               concell(n_concell)%ix = i1 - i
               concell(n_concell)%iz = i2 - j
            endif
         end do
      end if
   end do

else

   if (i2 == reg%itop%id) then
      isec => reg%itop
   else
      isec => reg%ibot
   endif
   n_concell = 0
   do n=1,4
      icell = isec%ccell_from_inode(n,i3)
      if (icell > 0) then
         n_concell = n_concell + 1
         concell(n_concell)%ix = isec%ccells(icell)%ix
         concell(n_concell)%iz = isec%ccells(icell)%iz
      endif
   enddo

endif

n_connode = 0

do n=1,n_concell

   concell_is_regular(n) = (i1 /= 0)

   if (associated(grid%ccind_from_2dc(concell(n)%ix,concell(n)%iz)%p)) then

      icell = grid%ccind_from_2dc(concell(n)%ix,concell(n)%iz)%p(reg%itop%iface_id)
      if (icell /= 0) then

         concell_is_regular(n) = .false.

         isec => reg%itop
         do jj=1,isec%n_inodes(icell)

            m = reg%itop_rnode_id(isec%inodes(jj,icell))

            if (reg%node_status(m) == 0) then
               n_connode = n_connode + 1
               connode(n_connode) = m
            endif

         end do

      end if

      icell = grid%ccind_from_2dc(concell(n)%ix,concell(n)%iz)%p(reg%ibot%iface_id)
      if (icell /= 0) then

         concell_is_regular(n) = .false.

         isec => reg%ibot
         do jj=1,isec%n_inodes(icell)

            m = reg%ibot_rnode_id(isec%inodes(jj,icell))
            if (reg%node_status(m) == 0) then
               n_connode = n_connode + 1
               connode(n_connode) = m
            endif

         end do

      endif

   endif

   do i=0,1
      ii = concell(n)%ix + i
      do j=0,1
         jj = concell(n)%iz + j

         ! Boundary check (defensive: find_connected_nodes2d has explicit checks)
         if (ii < 1 .or. ii > grid%nx .or. jj < 1 .or. jj > grid%nz) cycle

         if ((i1 /= 0 .and. abs(i1-ii)+abs(i2-jj) == 1) .or. (.not.concell_is_regular(n))) then

            if (associated(reg%region_rnode_id)) then
            if (reg%region_rnode_id(ii,jj) > 0) then

               m = reg%region_rnode_id(ii,jj)
               if (m > 0 .and. m <= reg%nnode) then
                 if (reg%node_status(m) == 0) then
                    n_connode = n_connode + 1
                    connode(n_connode) = m
                 endif
               endif

            endif
            endif

         endif

      end do
   end do

   concell_nnode(n) = n_connode

end do

concell_nbase(1) = 1
do n=1,n_concell

   concell_nnode(n) = concell_nnode(n) - concell_nbase(n) + 1
   if (n < n_concell) concell_nbase(n+1) = concell_nbase(n) + concell_nnode(n)

end do

return

end subroutine find_connected_cells2d


subroutine update_time2d(update_node)
use mod_prop2d
implicit none

integer :: update_node
integer :: n,i,j,ii,jj
integer :: i1,i2,i3
integer :: anode(30),nanode
integer :: id

i1 = reg%node(update_node)%i1
i2 = reg%node(update_node)%i2
i3 = reg%node(update_node)%i3

if (i1 /= 0) then

   if (grid%fully_regular(i1,i2)) then

      do j=-1,1,2
         do i=-1,1,2

            nanode = 0

            id = reg%region_rnode_id(i1+i,i2)
            if (id > 0 .and. id <= reg%nnode) then
              if (reg%node_status(id) == 0) then
                 nanode = nanode + 1
                 anode(nanode) = id
              endif
            endif
            id = reg%region_rnode_id(i1,i2+j)
            if (id > 0 .and. id <= reg%nnode) then
              if (reg%node_status(id) == 0) then
                 nanode = nanode + 1
                 anode(nanode) = id
              endif
            endif

            if (nanode > 0) then

               do n=1,nanode
                  if (anode(n) == newalive) then
                     anode(n) = anode(1)
                     anode(1) = newalive
                  endif
               end do

               if (anode(1) == newalive) call regular_update2d(update_node,anode(1:nanode),nanode)

            endif

         end do
      end do

      return

   endif

endif

call find_connected_cells2d(update_node)

do n=1,n_concell

   if (concell_nnode(n) > 0) then

      nanode = concell_nnode(n)

      if (nanode > 30) then
         print *,'nanode =',nanode
         stop 'subroutine update_time2d: nanode > 30'
      endif

      anode(1:nanode) = connode(concell_nbase(n):concell_nbase(n)+nanode-1)

      do i=1,nanode
         if (anode(i) == newalive) then
            j = anode(1)
            anode(1) = anode(i)
            anode(i) = j
         endif
      end do

      if (anode(1) == newalive) then

         if (concell_is_regular(n)) then

            call regular_update2d(update_node,anode(1:nanode),nanode)

         else

            select case (nanode)

            case(1)

               call time_from_1_node2d(update_node,anode(1))

            case(2)

               call time_from_1_node2d(update_node,anode(1))
               call time_from_2_nodes2d(update_node,anode(1),anode(2))

            case(3:)

               call time_from_1_node2d(update_node,anode(1))

               do ii=2,nanode
                  call time_from_2_nodes2d(update_node,anode(1),anode(ii))
               end do

            end select

         endif

      endif

   endif

end do

end subroutine update_time2d


subroutine regular_update2d(unode,anode,n_anode)
use mod_prop2d
implicit none

integer :: unode,n_anode,anode(n_anode)
integer :: iu,ju,ia,ja,ias,jas,di,dj,secnode,ioff,joff
real(kind=dp) :: a,b,c,t1,d2,det,tx,tz,t2x,t2z,velx,velz
integer :: n
logical :: x_info,z_info
real(kind=dp) :: tmax

if (n_anode > 2) stop 'too many nodes for regular update2d'

iu = reg%node(unode)%i1
ju = reg%node(unode)%i2

a = 0.0_dp
b = 0.0_dp
c = 0.0_dp

t2x = 0.0_dp
t2z = 0.0_dp

x_info = .false.
z_info = .false.

tmax = 0.0_dp

do n=1,n_anode

   tmax = max(tmax,reg%arrivaltime(anode(n)))

   ia = reg%node(anode(n))%i1
   ja = reg%node(anode(n))%i2

   if (ia == 0) cycle  ! skip interface nodes in regular update

   if (ia /= iu) then

      x_info = .true.
      ioff = ia
      di = ia - iu
      ias = iu + 2*di

      velx = 0.5_dp*(velocity(unode)+velocity(anode(n)))

      if (ias >= 1 .and. ias <= grid%nx) then

         if (reg%region_rnode_id(ias,ja) > 0) then

            secnode = reg%region_rnode_id(ias,ja)

            if (secnode > 0 .and. secnode <= reg%nnode) then
            if (reg%node_status(secnode) == 0 .and. reg%arrivaltime(anode(n)) >= reg%arrivaltime(secnode)) then

               if (velocity(anode(n))*velocity(secnode) > tiny(1.0_dp) .and. &
                   dot_product(reg%time_gradient(1:2,anode(n)),reg%time_gradient(1:2,secnode)) >= &
                    0.71_dp/(velocity(anode(n))*velocity(secnode))) then

                  tx = reg%arrivaltime(anode(n))
                  t2x = reg%arrivaltime(secnode)
                  d2 = (velx/grid%dx)**2
                  a = a + 2.25_dp*d2
                  b = b + (-6.0_dp*tx + 1.5_dp*t2x)*d2
                  c = c + (4.0_dp*tx*tx - 2.0_dp*tx*t2x + 0.25_dp*t2x*t2x)*d2

                  cycle

               endif

            endif
            endif

         endif

      endif

      tx = reg%arrivaltime(anode(n))
      d2 = (velx/grid%dx)**2
      a = a + d2
      b = b + (-2.0_dp*tx)*d2
      c = c + (tx*tx)*d2

   else if (ja /= ju) then

      z_info = .true.
      joff = ja
      dj = ja - ju
      jas = ju + 2*dj

      velz = 0.5_dp*(velocity(unode)+velocity(anode(n)))

      if (jas >= 1 .and. jas <= grid%nz) then

         if (reg%region_rnode_id(ia,jas) > 0) then

            secnode = reg%region_rnode_id(ia,jas)

            if (secnode > 0 .and. secnode <= reg%nnode) then
            if (reg%node_status(secnode) == 0 .and. reg%arrivaltime(anode(n)) >= reg%arrivaltime(secnode)) then

               if (velocity(anode(n))*velocity(secnode) > tiny(1.0_dp) .and. &
                   dot_product(reg%time_gradient(1:2,anode(n)),reg%time_gradient(1:2,secnode)) >= &
                    0.71_dp/(velocity(anode(n))*velocity(secnode))) then

                  tz = reg%arrivaltime(anode(n))
                  t2z = reg%arrivaltime(secnode)
                  d2 = (velz/grid%dz)**2
                  a = a + 2.25_dp*d2
                  b = b + (-6.0_dp*tz + 1.5_dp*t2z)*d2
                  c = c + (4.0_dp*tz*tz - 2.0_dp*tz*t2z + 0.25_dp*t2z*t2z)*d2

                  cycle

               endif

            endif
            endif

         endif

      endif

      tz = reg%arrivaltime(anode(n))
      d2 = (velz/grid%dz)**2
      a = a + d2
      b = b + (-2.0_dp*tz)*d2
      c = c + (tz*tz)*d2

   endif

end do

c = c - 1.0_dp

! Guard against division by zero when a==0 (degenerate case: no valid neighbors)
if (abs(a) < tiny(1.0_dp)) return

det = b*b - 4.0_dp*a*c

if (det < 0.0_dp) then

   if (0.5_dp*sqrt(abs(det))/abs(a) < 1.0e20_dp) then
      det = 0.0_dp
      n_det_negative = n_det_negative + 1
   else
      print *,'negative determinant too large in regular_update2d'
      print *,a,b,c
      stop 'negative determinant too large'
   endif

endif

t1 = 0.5_dp*(-b + sqrt(det))/a

t1 = max(t1,tmax)

if (t1 < reg%arrivaltime(unode)) then

   reg%arrivaltime(unode) = t1

   reg%time_gradient(1:2,unode) = 0.0_dp

   if (x_info) then
      if (t2x == 0.0_dp) then
         reg%time_gradient(1,unode) = (t1-tx)/(grid%x(iu)-grid%x(ioff))
      else
         reg%time_gradient(1,unode) = 0.5_dp*(3.0_dp*t1-4.0_dp*tx+t2x)/(grid%x(iu)-grid%x(ioff))
      endif
   endif

   if (z_info) then
      if (t2z == 0.0_dp) then
         reg%time_gradient(2,unode) = (t1-tz)/(grid%z(ju)-grid%z(joff))
      else
         reg%time_gradient(2,unode) = 0.5_dp*(3.0_dp*t1-4.0_dp*tz+t2z)/(grid%z(ju)-grid%z(joff))
      endif
   endif

   suspect_time_gradient(unode) = .not.(x_info .and. z_info)

endif

return

end subroutine regular_update2d


subroutine time_from_1_node2d(unode,node1)
use mod_prop2d
implicit none

integer :: unode,node1
real(kind=dp) :: atime,ttime,dist
real(kind=dp) :: xu(2),x1(2)

xu(1) = reg%x(unode)
xu(2) = reg%z(unode)

x1(1) = reg%x(node1)
x1(2) = reg%z(node1)

dist = sqrt(sum((xu-x1)**2))

if (dist < 0.01_dp*grid%tolerance) then
   reg%arrivaltime(unode) = reg%arrivaltime(node1)
   reg%time_gradient(1:2,unode) = reg%time_gradient(1:2,node1)
   return
endif

ttime = dist/(0.5_dp*(velocity(unode)+velocity(node1)))
atime = reg%arrivaltime(node1) + ttime

if (atime < reg%arrivaltime(unode)) then

   reg%arrivaltime(unode) = atime

   reg%time_gradient(1:2,unode) = (xu-x1)/(dist*velocity(unode))

endif

return

end subroutine time_from_1_node2d


subroutine time_from_2_nodes2d(unode,anode1,anode2)
use mod_prop2d
implicit none

integer :: unode,node1,node2,anode1,anode2
real(kind=dp) :: atime

real(kind=dp) :: wn1(2),wn2(2)
real(kind=dp) :: d(2),p(2),q(2),u1,u2,det,ae,be,ce
real(kind=dp) :: n1(2),n2(2),r1(2),r2(2),xu(2),x1(2),x2(2)
real(kind=dp) :: v12

integer :: vscount

node1 = anode1
node2 = anode2

if (reg%arrivaltime(node1) < reg%arrivaltime(node2)) then
   vscount = node1
   node1 = node2
   node2 = vscount
endif

xu(1) = reg%x(unode)
xu(2) = reg%z(unode)

x1(1) = reg%x(node1)
x1(2) = reg%z(node1)

x2(1) = reg%x(node2)
x2(2) = reg%z(node2)

r1 = x1 - xu
r2 = x2 - xu

if (sqrt(sum(r1**2)) < 0.01_dp*grid%tolerance) then
   reg%arrivaltime(unode) = reg%arrivaltime(node1)
   reg%time_gradient(1:2,unode) = reg%time_gradient(1:2,node1)
   return
endif

if (sqrt(sum(r2**2)) < 0.01_dp*grid%tolerance) then
   reg%arrivaltime(unode) = reg%arrivaltime(node2)
   reg%time_gradient(1:2,unode) = reg%time_gradient(1:2,node2)
   return
endif

if (sqrt(sum((x2-x1)**2)) < 0.01_dp*grid%tolerance) then
   call time_from_1_node2d(unode,node1)
   return
endif

v12 = 0.5_dp*(velocity(node1)+velocity(node2))

d = (x1-x2)/v12

p = (reg%arrivaltime(node1)-reg%arrivaltime(node2)) / sum(d*d) * d

q(1) = -d(2)
q(2) = d(1)

ae = sum(q*q)
be = 2.0_dp*sum(p*q)
ce = sum(p*p) - 1.0_dp

det = be*be - 4.0_dp*ae*ce

if (det < 0.0_dp) return

u1 = 0.5_dp*(-be+sqrt(det))/ae
u2 = 0.5_dp*(-be-sqrt(det))/ae

wn1 = p + u1*q
wn2 = p + u2*q

n1(1) = -r1(2)
n1(2) = r1(1)
n1 = n1/sqrt(sum(n1**2))
if (dot_product(n1,r2) < 0.0_dp) n1 = -n1

n2(1) = -r2(2)
n2(2) = r2(1)
n2 = n2/sqrt(sum(n2**2))
if (dot_product(n2,r1) < 0.0_dp) n2 = -n2

vscount = 0

if (dot_product(wn1,n1) <= 0.0_dp .and. dot_product(wn1,n2) <= 0.0_dp) then

   atime = reg%arrivaltime(node1) - dot_product(wn1,r1)/(0.5_dp*(velocity(unode)+velocity(node1)))
   if (atime >= reg%arrivaltime(node1)) then
      if (atime < reg%arrivaltime(unode)) then
         reg%arrivaltime(unode) = atime
         reg%time_gradient(1:2,unode) = wn1/(sqrt(sum(wn1**2))*velocity(node1))
      endif
      vscount = vscount + 1
   endif

end if

if (dot_product(wn2,n1) <= 0.0_dp .and. dot_product(wn2,n2) <= 0.0_dp) then

   atime = reg%arrivaltime(node1) - dot_product(wn2,r1)/(0.5_dp*(velocity(unode)+velocity(node1)))
   if (atime >= reg%arrivaltime(node1)) then
      if (atime < reg%arrivaltime(unode)) then
         reg%arrivaltime(unode) = atime
         reg%time_gradient(1:2,unode) = wn2/(sqrt(sum(wn2**2))*velocity(node1))
      endif
      vscount = vscount + 1
   endif

end if

return

end subroutine time_from_2_nodes2d


subroutine fit_gradient_at_node2d(reg,centernode,dtdx,dtdz)
use mod_fmt2d_nointerfaces
implicit none

type(Tregion2d) :: reg

integer,parameter :: maxfitnodes = 30
integer :: n,m,i,j,i1,i2,i3,n_cnode
integer :: ii,jj,nfit,icell
integer :: centernode,n_concell_loc,mstore(maxfitnodes)
type(Tinteger_coordinates2d) :: concell_loc(4)
logical :: concell_irreg(4)
real(kind=dp) :: rx(maxfitnodes),rz(maxfitnodes),atime(maxfitnodes)
real(kind=dp) :: dtdx,dtdz,dx,dz,dw

real(kind=dp) :: ata(2,2),atb(2),sol(2),det_loc

type(Tintersection2d),pointer :: isec
type(Tpropagation_grid2d),pointer :: lgrid

lgrid => reg%grid

i1 = reg%node(centernode)%i1
i2 = reg%node(centernode)%i2
i3 = reg%node(centernode)%i3

if (i1 /= 0) then

   n_concell_loc = 0
   do i=0,1
      if ((i1-i > 0) .and. (i1-i < lgrid%nx)) then
         do j=0,1
            if ((i2-j > 0) .and. (i2-j < lgrid%nz)) then
               n_concell_loc = n_concell_loc + 1
               concell_loc(n_concell_loc)%ix = i1 - i
               concell_loc(n_concell_loc)%iz = i2 - j
            endif
         end do
      end if
   end do

else

   if (i2 == reg%itop%id) then
      isec => reg%itop
   else
      isec => reg%ibot
   endif
   n_concell_loc = 0
   do n=1,4
      icell = isec%ccell_from_inode(n,i3)
      if (icell > 0) then
         n_concell_loc = n_concell_loc + 1
         concell_loc(n_concell_loc)%ix = isec%ccells(icell)%ix
         concell_loc(n_concell_loc)%iz = isec%ccells(icell)%iz
      endif
   enddo

endif

do n=1,n_concell_loc
   concell_irreg(n) = associated(lgrid%ccind_from_2dc(concell_loc(n)%ix,concell_loc(n)%iz)%p)
end do

if (count(concell_irreg(1:n_concell_loc)) == 0) then

   if (i1 /= 0) then

      m = reg%region_rnode_id(max(1,i1-1),i2)
      n = reg%region_rnode_id(min(lgrid%nx,i1+1),i2)
      if (m > 0 .and. m <= reg%nnode .and. n > 0 .and. n <= reg%nnode) then
        dtdx = (reg%arrivaltime(n)-reg%arrivaltime(m))/(reg%x(n)-reg%x(m))
      else if (m > 0 .and. m <= reg%nnode) then
        ! One-sided: only left neighbor available
        dtdx = (reg%arrivaltime(centernode)-reg%arrivaltime(m))/(reg%x(centernode)-reg%x(m))
      else if (n > 0 .and. n <= reg%nnode) then
        ! One-sided: only right neighbor available
        dtdx = (reg%arrivaltime(n)-reg%arrivaltime(centernode))/(reg%x(n)-reg%x(centernode))
      else
        dtdx = 0.0_dp
      endif

      m = reg%region_rnode_id(i1,max(1,i2-1))
      n = reg%region_rnode_id(i1,min(lgrid%nz,i2+1))
      if (m > 0 .and. m <= reg%nnode .and. n > 0 .and. n <= reg%nnode) then
        dtdz = (reg%arrivaltime(n)-reg%arrivaltime(m))/(reg%z(n)-reg%z(m))
      else if (m > 0 .and. m <= reg%nnode) then
        ! One-sided: only upper neighbor available
        dtdz = (reg%arrivaltime(centernode)-reg%arrivaltime(m))/(reg%z(centernode)-reg%z(m))
      else if (n > 0 .and. n <= reg%nnode) then
        ! One-sided: only lower neighbor available
        dtdz = (reg%arrivaltime(n)-reg%arrivaltime(centernode))/(reg%z(n)-reg%z(centernode))
      else
        dtdz = 0.0_dp
      endif

      return

   endif

   if (n_concell_loc == 0) then
      dtdx = 0.0_dp
      dtdz = 0.0_dp
      return
   endif

endif

n_cnode = 0

do n=1,n_concell_loc

   if (concell_irreg(n)) then

      if (associated(lgrid%ccind_from_2dc(concell_loc(n)%ix,concell_loc(n)%iz)%p)) then

         icell = lgrid%ccind_from_2dc(concell_loc(n)%ix,concell_loc(n)%iz)%p(reg%itop%iface_id)
         if (icell /= 0) then
            isec => reg%itop
            do jj=1,isec%n_inodes(icell)
               m = reg%itop_rnode_id(isec%inodes(jj,icell))
               if (m == centernode) cycle
               if (n_cnode == 0) then
                  n_cnode = 1
                  mstore(n_cnode) = m
                  rx(n_cnode) = reg%x(m)
                  rz(n_cnode) = reg%z(m)
                  atime(n_cnode) = reg%arrivaltime(m)
               else
                  if (count(mstore(1:n_cnode) == m) == 0) then
                     if (n_cnode < maxfitnodes) then
                     n_cnode = n_cnode + 1
                     mstore(n_cnode) = m
                     rx(n_cnode) = reg%x(m)
                     rz(n_cnode) = reg%z(m)
                     atime(n_cnode) = reg%arrivaltime(m)
                     endif
                  endif
               endif
            end do
         endif

         icell = lgrid%ccind_from_2dc(concell_loc(n)%ix,concell_loc(n)%iz)%p(reg%ibot%iface_id)
         if (icell /= 0) then
            isec => reg%ibot
            do jj=1,isec%n_inodes(icell)
               m = reg%ibot_rnode_id(isec%inodes(jj,icell))
               if (m == centernode) cycle
               if (n_cnode == 0) then
                  n_cnode = 1
                  mstore(n_cnode) = m
                  rx(n_cnode) = reg%x(m)
                  rz(n_cnode) = reg%z(m)
                  atime(n_cnode) = reg%arrivaltime(m)
               else
                  if (count(mstore(1:n_cnode) == m) == 0) then
                     if (n_cnode < maxfitnodes) then
                     n_cnode = n_cnode + 1
                     mstore(n_cnode) = m
                     rx(n_cnode) = reg%x(m)
                     rz(n_cnode) = reg%z(m)
                     atime(n_cnode) = reg%arrivaltime(m)
                     endif
                  endif
               endif
            end do
         endif

      endif

   endif

   do i=0,1
      ii = concell_loc(n)%ix + i
      if (ii < 1 .or. ii > lgrid%nx) cycle
      do j=0,1
         jj = concell_loc(n)%iz + j
         if (jj < 1 .or. jj > lgrid%nz) cycle

         if (reg%region_rnode_id(ii,jj) > 0) then

            m = reg%region_rnode_id(ii,jj)
            if (m <= 0 .or. m > reg%nnode) cycle
            if (m == centernode) cycle

            if (n_cnode == 0) then
               n_cnode = 1
               mstore(n_cnode) = m
               rx(n_cnode) = reg%x(m)
               rz(n_cnode) = reg%z(m)
               atime(n_cnode) = reg%arrivaltime(m)
            else
               if (count(mstore(1:n_cnode) == m) == 0) then
                  if (n_cnode < maxfitnodes) then
                  n_cnode = n_cnode + 1
                  mstore(n_cnode) = m
                  rx(n_cnode) = reg%x(m)
                  rz(n_cnode) = reg%z(m)
                  atime(n_cnode) = reg%arrivaltime(m)
                  endif
               endif
            endif

         endif

      end do
   end do

end do

nfit = 2

ata = 0.0_dp
atb = 0.0_dp

do i=1,n_cnode

   dx = rx(i) - reg%x(centernode)
   dz = rz(i) - reg%z(centernode)
   dw = 1.0_dp/(sqrt(dx*dx+dz*dz) + lgrid%tolerance)

   ata(1,1) = ata(1,1) + dx*dx*dw*dw
   ata(1,2) = ata(1,2) + dx*dz*dw*dw
   ata(2,2) = ata(2,2) + dz*dz*dw*dw
   atb(1) = atb(1) + (atime(i)-reg%arrivaltime(centernode))*dx*dw*dw
   atb(2) = atb(2) + (atime(i)-reg%arrivaltime(centernode))*dz*dw*dw

end do

ata(2,1) = ata(1,2)

det_loc = ata(1,1)*ata(2,2) - ata(1,2)*ata(2,1)

if (abs(det_loc) < 1.0e-30_dp) then
   dtdx = 0.0_dp
   dtdz = 0.0_dp
   return
endif

sol(1) = (ata(2,2)*atb(1) - ata(1,2)*atb(2))/det_loc
sol(2) = (ata(1,1)*atb(2) - ata(2,1)*atb(1))/det_loc

dtdx = sol(1)
dtdz = sol(2)

return

end subroutine fit_gradient_at_node2d


subroutine add_node_to_tree(nn)
use mod_prop2d
implicit none

integer :: nn,parent_index,child_index,temp

ntree = ntree + 1

reg%node_status(nn) = ntree
node_from_tree_ind(ntree) = nn

child_index = ntree
parent_index = ntree/2

do while (parent_index > 0)

   if (reg%arrivaltime(nn) < reg%arrivaltime(node_from_tree_ind(parent_index))) then

      reg%node_status(nn) = parent_index
      reg%node_status(node_from_tree_ind(parent_index)) = child_index
      temp = node_from_tree_ind(child_index)
      node_from_tree_ind(child_index) = node_from_tree_ind(parent_index)
      node_from_tree_ind(parent_index) = temp
      child_index = parent_index
      parent_index = child_index/2

   else

      parent_index = 0

   endif

end do

return

end subroutine add_node_to_tree


subroutine remove_root_from_tree
use mod_prop2d
implicit none

integer :: parent_index,child_index,temp

if (ntree == 1) then
   ntree = 0
   return
endif

reg%node_status(node_from_tree_ind(ntree)) = 1
node_from_tree_ind(1) = node_from_tree_ind(ntree)

ntree = ntree - 1

parent_index = 1
child_index = 2*parent_index

do while (child_index < ntree)

   if (reg%arrivaltime(node_from_tree_ind(child_index)) > &
       reg%arrivaltime(node_from_tree_ind(child_index+1))) &
       child_index = child_index + 1

   if (reg%arrivaltime(node_from_tree_ind(child_index)) < &
       reg%arrivaltime(node_from_tree_ind(parent_index))) then

      reg%node_status(node_from_tree_ind(parent_index)) = child_index
      reg%node_status(node_from_tree_ind(child_index)) = parent_index
      temp = node_from_tree_ind(child_index)
      node_from_tree_ind(child_index) = node_from_tree_ind(parent_index)
      node_from_tree_ind(parent_index) = temp
      parent_index = child_index
      child_index = 2*parent_index

   else

      child_index = ntree + 1

   endif

end do

if (child_index == ntree) then

   if (reg%arrivaltime(node_from_tree_ind(child_index)) < &
       reg%arrivaltime(node_from_tree_ind(parent_index))) then

      reg%node_status(node_from_tree_ind(parent_index)) = child_index
      reg%node_status(node_from_tree_ind(child_index)) = parent_index
      temp = node_from_tree_ind(child_index)
      node_from_tree_ind(child_index) = node_from_tree_ind(parent_index)
      node_from_tree_ind(parent_index) = temp

   endif

endif

return

end subroutine remove_root_from_tree


subroutine update_tree(nn)
use mod_prop2d
implicit none

integer :: nn,temp,parent_index,child_index

child_index = reg%node_status(nn)
parent_index = child_index/2

do while (parent_index > 0)
   if (reg%arrivaltime(nn) < reg%arrivaltime(node_from_tree_ind(parent_index))) then

      reg%node_status(nn) = parent_index
      reg%node_status(node_from_tree_ind(parent_index)) = child_index
      temp = node_from_tree_ind(child_index)
      node_from_tree_ind(child_index) = node_from_tree_ind(parent_index)
      node_from_tree_ind(parent_index) = temp
      child_index = parent_index
      parent_index = child_index/2

   else

      parent_index = 0

   endif

end do

return

end subroutine update_tree


subroutine find_intersection2d(isec,iface,grid)
use mod_fmt2d_nointerfaces
implicit none

type(Tintersection2d) :: isec
type(Tinterface2d) :: iface
type(Tpropagation_grid2d),target :: grid

integer :: ix,iz,n,icell,ncc,inode_count
integer :: iz_cell
real(kind=dp) :: z_val,dzdx,norm_len
integer,allocatable :: temp_ccell(:,:)
integer :: max_ccells

isec%grid => grid
isec%id = iface%id
isec%iface_id = iface%id
isec%pinched = iface%pinched

isec%nnode = iface%nx

allocate(isec%x(isec%nnode))
allocate(isec%z(isec%nnode))
allocate(isec%normal(2,isec%nnode))

isec%x = iface%x
isec%z = iface%z_nodes

do ix=1,isec%nnode

   if (ix > 1 .and. ix < isec%nnode) then
      dzdx = (isec%z(ix+1) - isec%z(ix-1))/(isec%x(ix+1) - isec%x(ix-1))
   else if (ix == 1 .and. isec%nnode >= 2) then
      dzdx = (isec%z(ix+1) - isec%z(ix))/(isec%x(ix+1) - isec%x(ix))
   else if (ix == isec%nnode .and. isec%nnode >= 2) then
      dzdx = (isec%z(ix) - isec%z(ix-1))/(isec%x(ix) - isec%x(ix-1))
   else
      dzdx = 0.0_dp
   endif

   norm_len = sqrt(1.0_dp + dzdx**2)
   isec%normal(1,ix) = -dzdx/norm_len
   isec%normal(2,ix) = 1.0_dp/norm_len

end do

max_ccells = isec%nnode
allocate(temp_ccell(2,max_ccells))

ncc = 0
do ix=1,isec%nnode

   z_val = isec%z(ix)

   iz_cell = 0
   do iz=1,grid%nz-1
      if (grid%z(iz) <= z_val .and. z_val < grid%z(iz+1)) then
         iz_cell = iz
         exit
      endif
   end do

   if (iz_cell == 0) then
      if (z_val <= grid%z(1)) iz_cell = 1
      if (z_val >= grid%z(grid%nz)) iz_cell = grid%nz - 1
   endif

   if (iz_cell > 0) then
      ncc = ncc + 1
      if (ncc > max_ccells) then
         print *,'find_intersection2d: too many cut cells'
         stop
      endif
      temp_ccell(1,ncc) = ix
      temp_ccell(2,ncc) = iz_cell
   endif

end do

isec%n_ccells = ncc

if (ncc > 0) then
   allocate(isec%ccells(ncc))
   allocate(isec%n_inodes(ncc))
   allocate(isec%inodes(4,ncc))

   do icell=1,ncc
      isec%ccells(icell)%ix = temp_ccell(1,icell)
      isec%ccells(icell)%iz = temp_ccell(2,icell)
   end do

   do icell=1,ncc
      isec%n_inodes(icell) = 0
   end do

   do ix=1,isec%nnode
      z_val = isec%z(ix)

      do icell=1,ncc
         if (isec%ccells(icell)%ix == ix) then
            iz_cell = isec%ccells(icell)%iz
            if (grid%z(iz_cell) <= z_val .and. z_val <= grid%z(iz_cell+1)) then
               isec%n_inodes(icell) = isec%n_inodes(icell) + 1
               if (isec%n_inodes(icell) > 4) then
                  print *,'find_intersection2d: too many inodes in cell'
                  stop
               endif
               isec%inodes(isec%n_inodes(icell),icell) = ix
            endif
         endif
      end do
   end do

else

   allocate(isec%ccells(1))
   allocate(isec%n_inodes(1))
   allocate(isec%inodes(4,1))
   isec%n_inodes(1) = 0

endif

allocate(isec%ccell_from_inode(4,isec%nnode))
isec%ccell_from_inode = 0

do ix=1,isec%nnode
   inode_count = 0
   do icell=1,ncc
      if (isec%ccells(icell)%ix == ix) then
         inode_count = inode_count + 1
         if (inode_count <= 4) then
            isec%ccell_from_inode(inode_count,ix) = icell
         endif
      endif
   end do
end do

if (.not.associated(grid%ccind_from_2dc)) then
   allocate(grid%ccind_from_2dc(grid%nx,grid%nz))
   do ix=1,grid%nx
      do iz=1,grid%nz
         nullify(grid%ccind_from_2dc(ix,iz)%p)
      end do
   end do
endif

do icell=1,ncc
   ix = isec%ccells(icell)%ix
   iz = isec%ccells(icell)%iz

   if (ix >= 1 .and. ix < grid%nx .and. iz >= 1 .and. iz < grid%nz) then

      if (.not.associated(grid%ccind_from_2dc(ix,iz)%p)) then
         allocate(grid%ccind_from_2dc(ix,iz)%p(n_interfaces))
         grid%ccind_from_2dc(ix,iz)%p = 0
      endif

      grid%ccind_from_2dc(ix,iz)%p(isec%iface_id) = icell

   endif

end do

deallocate(temp_ccell)

return

end subroutine find_intersection2d


subroutine define_region2d(reg,itop,ibot,grid)
use mod_fmt2d_nointerfaces
implicit none

type(Tregion2d),target :: reg
type(Tintersection2d),target :: itop,ibot
type(Tpropagation_grid2d),target :: grid

integer :: ix,iz,n,m
integer :: ngrid_nodes,ntop_nodes,nbot_nodes
real(kind=dp) :: z_top,z_bot

reg%grid => grid
reg%itop => itop
reg%ibot => ibot

ngrid_nodes = 0
do ix=1,grid%nx
   do iz=1,grid%nz
      z_top = interpolate_iface_z2d(itop,grid%x(ix))
      z_bot = interpolate_iface_z2d(ibot,grid%x(ix))
      if (z_top <= grid%z(iz) .and. grid%z(iz) <= z_bot) then
         ngrid_nodes = ngrid_nodes + 1
      endif
   end do
end do

ntop_nodes = itop%nnode
nbot_nodes = ibot%nnode

reg%ngnode = ngrid_nodes
reg%nnode = ngrid_nodes + ntop_nodes + nbot_nodes

allocate(reg%node(reg%nnode))
allocate(reg%node_status(reg%nnode))
allocate(reg%arrivaltime(reg%nnode))
allocate(reg%time_gradient(2,reg%nnode))
allocate(reg%velocity(NVTYPE,reg%nnode))
allocate(reg%x(reg%nnode))
allocate(reg%z(reg%nnode))
allocate(reg%region_rnode_id(grid%nx,grid%nz))

reg%node_status = -1
reg%arrivaltime = huge_time
reg%time_gradient = 0.0_dp
reg%loaded_tf => null()
reg%velocity = 0.0_dp
reg%region_rnode_id = 0

if (.not.associated(grid%rnode_id)) then
   allocate(grid%rnode_id(grid%nx,grid%nz))
   grid%rnode_id = 0
endif

if (.not.associated(grid%node_region)) then
   allocate(grid%node_region(grid%nx,grid%nz))
   grid%node_region = 0
endif

n = 0
do ix=1,grid%nx
   do iz=1,grid%nz
      z_top = interpolate_iface_z2d(itop,grid%x(ix))
      z_bot = interpolate_iface_z2d(ibot,grid%x(ix))
      if (z_top <= grid%z(iz) .and. grid%z(iz) <= z_bot) then
         n = n + 1
         reg%node(n)%i1 = ix
         reg%node(n)%i2 = iz
         reg%node(n)%i3 = 0
         reg%x(n) = grid%x(ix)
         reg%z(n) = grid%z(iz)
         grid%rnode_id(ix,iz) = n
         grid%node_region(ix,iz) = reg%id
         reg%region_rnode_id(ix,iz) = n
      else
         if (grid%node_region(ix,iz) == 0) grid%node_region(ix,iz) = -1
      endif
   end do
end do

allocate(itop%rright_node_id(itop%nnode))
allocate(reg%itop_rnode_id(itop%nnode))
do m=1,itop%nnode
   n = n + 1
   reg%node(n)%i1 = 0
   reg%node(n)%i2 = itop%id
   reg%node(n)%i3 = m
   reg%x(n) = itop%x(m)
   reg%z(n) = itop%z(m)
   itop%rright_node_id(m) = n
   reg%itop_rnode_id(m) = n
end do

allocate(ibot%rleft_node_id(ibot%nnode))
allocate(reg%ibot_rnode_id(ibot%nnode))
do m=1,ibot%nnode
   n = n + 1
   reg%node(n)%i1 = 0
   reg%node(n)%i2 = ibot%id
   reg%node(n)%i3 = m
   reg%x(n) = ibot%x(m)
   reg%z(n) = ibot%z(m)
   ibot%rleft_node_id(m) = n
   reg%ibot_rnode_id(m) = n
end do

return

contains

function interpolate_iface_z2d(isec,px) result(z_val)
type(Tintersection2d) :: isec
real(kind=dp) :: px,z_val
integer :: ixx
real(kind=dp) :: w

if (isec%nnode < 2) then
   z_val = isec%z(1)
   return
endif

ixx = floor((px - isec%x(1))/isec%grid%dx) + 1

if (ixx < 1) then
   z_val = isec%z(1)
   return
endif
if (ixx >= isec%nnode) then
   z_val = isec%z(isec%nnode)
   return
endif

w = (px - isec%x(ixx))/isec%grid%dx
if (w < 0.0_dp) w = 0.0_dp
if (w > 1.0_dp) w = 1.0_dp

z_val = (1.0_dp - w)*isec%z(ixx) + w*isec%z(ixx+1)

return

end function interpolate_iface_z2d

end subroutine define_region2d


subroutine velocities_to_grid2d(pgrid_in)
use mod_fmt2d_nointerfaces
implicit none

type(Tpropagation_grid2d),target :: pgrid_in

integer :: ix,iz,vtype
type(Tvelocity_grid2d),pointer :: vg

do vtype=1,NVTYPE

   vg => vgrid(pgrid_in%index_x0,vtype)

   if (.not.associated(vg%velocity)) cycle

   if (.not.associated(pgrid_in%velocity)) then
      allocate(pgrid_in%velocity(NVTYPE,pgrid_in%nx,pgrid_in%nz))
      pgrid_in%velocity = 0.0_dp
   endif

   if (vg%nx == pgrid_in%nx .and. vg%nz == pgrid_in%nz .and. &
       abs(vg%dx - pgrid_in%dx) < 1.0e-10_dp .and. &
       abs(vg%dz - pgrid_in%dz) < 1.0e-10_dp .and. &
       abs(vg%x(1) - pgrid_in%x(1)) < 1.0e-10_dp .and. &
       abs(vg%z(1) - pgrid_in%z(1)) < 1.0e-10_dp) then

      do ix=1,pgrid_in%nx
         do iz=1,pgrid_in%nz
            pgrid_in%velocity(vtype,ix,iz) = vg%velocity(vtype,ix,iz)
         end do
      end do

   else

      do ix=1,pgrid_in%nx
         do iz=1,pgrid_in%nz
            pgrid_in%velocity(vtype,ix,iz) = vg%velocity(vtype,1,1)
         end do
      end do

   endif

end do

return

end subroutine velocities_to_grid2d


subroutine velocities_to_intersection2d(isec)
use mod_fmt2d_nointerfaces
implicit none

type(Tintersection2d) :: isec

integer :: m,vtype
integer :: ivx,ivz
real(kind=dp) :: wx,wz,px,pz,dz_lower,dz_upper
real(kind=dp) :: pz_off, delta_z
type(Tvelocity_grid2d),pointer :: vg
type(Tpropagation_grid2d),pointer :: lgrid

lgrid => isec%grid

if (.not.associated(isec%vel_top)) then
   allocate(isec%vel_top(NVTYPE,isec%nnode))
   allocate(isec%vel_bot(NVTYPE,isec%nnode))
endif

do vtype=1,NVTYPE

   vg => vgrid(lgrid%index_x0,vtype)

   if (.not.associated(vg%velocity)) cycle

   delta_z = vg%dz

   do m=1,isec%nnode

      px = isec%x(m)
      pz = isec%z(m)

      ivx = floor((px - vg%x(1))/vg%dx) + 1
      if (ivx < 1) ivx = 1
      if (ivx >= vg%nx) ivx = vg%nx - 1

      wx = (px - vg%x(ivx))/vg%dx
      if (wx < 0.0_dp) wx = 0.0_dp
      if (wx > 1.0_dp) wx = 1.0_dp

      pz_off = pz - delta_z
      if (pz_off < vg%z(1)) then
         ivz = 1
         wz = 0.0_dp
      else if (pz_off >= vg%z(vg%nz)) then
         ivz = vg%nz - 1
         wz = 1.0_dp
      else
         do ivz = 1, vg%nz - 1
            if (pz_off >= vg%z(ivz) .and. pz_off < vg%z(ivz+1)) exit
         end do
         if (ivz >= vg%nz) ivz = vg%nz - 1
         dz_lower = pz_off - vg%z(ivz)
         dz_upper = vg%z(ivz+1) - vg%z(ivz)
         if (dz_upper > 0.0_dp) then
            wz = dz_lower / dz_upper
         else
            wz = 0.0_dp
         endif
      endif
      if (wz < 0.0_dp) wz = 0.0_dp
      if (wz > 1.0_dp) wz = 1.0_dp

      isec%vel_top(vtype,m) = &
           (1.0_dp-wx)*(1.0_dp-wz)*vg%velocity(vtype,ivx,ivz) + &
           wx*(1.0_dp-wz)*vg%velocity(vtype,ivx+1,ivz) + &
           (1.0_dp-wx)*wz*vg%velocity(vtype,ivx,ivz+1) + &
           wx*wz*vg%velocity(vtype,ivx+1,ivz+1)

      pz_off = pz + delta_z
      if (pz_off < vg%z(1)) then
         ivz = 1
         wz = 0.0_dp
      else if (pz_off >= vg%z(vg%nz)) then
         ivz = vg%nz - 1
         wz = 1.0_dp
      else
         do ivz = 1, vg%nz - 1
            if (pz_off >= vg%z(ivz) .and. pz_off < vg%z(ivz+1)) exit
         end do
         if (ivz >= vg%nz) ivz = vg%nz - 1
         dz_lower = pz_off - vg%z(ivz)
         dz_upper = vg%z(ivz+1) - vg%z(ivz)
         if (dz_upper > 0.0_dp) then
            wz = dz_lower / dz_upper
         else
            wz = 0.0_dp
         endif
      endif
      if (wz < 0.0_dp) wz = 0.0_dp
      if (wz > 1.0_dp) wz = 1.0_dp

      isec%vel_bot(vtype,m) = &
           (1.0_dp-wx)*(1.0_dp-wz)*vg%velocity(vtype,ivx,ivz) + &
           wx*(1.0_dp-wz)*vg%velocity(vtype,ivx+1,ivz) + &
           (1.0_dp-wx)*wz*vg%velocity(vtype,ivx,ivz+1) + &
           wx*wz*vg%velocity(vtype,ivx+1,ivz+1)

   end do

end do

return

end subroutine velocities_to_intersection2d


subroutine velocities_to_region2d(reg,pgrid_in)
use mod_fmt2d_nointerfaces
implicit none

type(Tregion2d) :: reg
type(Tpropagation_grid2d),target :: pgrid_in

integer :: n,vtype
integer :: i1,i2,i3
integer :: ivx,ivz
real(kind=dp) :: wx,wz,px,pz
type(Tvelocity_grid2d),pointer :: vg

do vtype=1,NVTYPE

   vg => vgrid(reg%ivgrid,vtype)

   if (.not.associated(vg%velocity)) cycle

   do n=1,reg%nnode

      i1 = reg%node(n)%i1
      i2 = reg%node(n)%i2
      i3 = reg%node(n)%i3

      if (i1 /= 0) then

         if (associated(pgrid_in%velocity)) then
            reg%velocity(vtype,n) = pgrid_in%velocity(vtype,i1,i2)
         else

            px = reg%x(n)
            pz = reg%z(n)

            ivx = floor((px - vg%x(1))/vg%dx) + 1
            ivz = floor((pz - vg%z(1))/vg%dz) + 1

            if (ivx < 1) ivx = 1
            if (ivx >= vg%nx) ivx = vg%nx - 1
            if (ivz < 1) ivz = 1
            if (ivz >= vg%nz) ivz = vg%nz - 1

            wx = (px - vg%x(ivx))/vg%dx
            wz = (pz - vg%z(ivz))/vg%dz

            if (wx < 0.0_dp) wx = 0.0_dp
            if (wx > 1.0_dp) wx = 1.0_dp
            if (wz < 0.0_dp) wz = 0.0_dp
            if (wz > 1.0_dp) wz = 1.0_dp

            reg%velocity(vtype,n) = &
                 (1.0_dp-wx)*(1.0_dp-wz)*vg%velocity(vtype,ivx,ivz) + &
                 wx*(1.0_dp-wz)*vg%velocity(vtype,ivx+1,ivz) + &
                 (1.0_dp-wx)*wz*vg%velocity(vtype,ivx,ivz+1) + &
                 wx*wz*vg%velocity(vtype,ivx+1,ivz+1)

         endif

      else

         if (i2 == reg%itop%id) then
            reg%velocity(vtype,n) = reg%itop%vel_bot(vtype,i3)
         else if (i2 == reg%ibot%id) then
            reg%velocity(vtype,n) = reg%ibot%vel_top(vtype,i3)
         endif

      endif

   end do

end do

! Second pass: fix velocity at boundary grid nodes.
! Grid nodes at the bottom boundary (no node below in this region) get
! vel_top of the bottom interface (velocity from the region's side).
! Grid nodes at the top boundary (no node above in this region) get
! vel_bot of the top interface (velocity from the region's side).
! This is critical for the water layer where the shared velocity grid
! assigns sediment velocity to seafloor-level nodes.

do n=1,reg%nnode
   i1 = reg%node(n)%i1
   if (i1 == 0) cycle  ! interface nodes already correct

   i2 = reg%node(n)%i2

   ! Check if this is a bottom boundary grid node
   if (i2 < pgrid_in%nz) then
      if (reg%region_rnode_id(i1, i2+1) <= 0) then
         ! Bottom boundary: use vel_top of bottom interface
         if (associated(reg%ibot%vel_top)) then
            i3 = nint((reg%x(n) - reg%ibot%x(1)) / pgrid_in%dx) + 1
            if (i3 < 1) i3 = 1
            if (i3 > reg%ibot%nnode) i3 = reg%ibot%nnode
            do vtype=1,NVTYPE
               if (associated(vgrid(reg%ivgrid,vtype)%velocity)) then
                  reg%velocity(vtype,n) = reg%ibot%vel_top(vtype,i3)
               endif
            end do
         endif
      endif
   endif

   ! Check if this is a top boundary grid node
   if (i2 > 1) then
      if (reg%region_rnode_id(i1, i2-1) <= 0) then
         ! Top boundary: use vel_bot of top interface
         if (associated(reg%itop%vel_bot)) then
            i3 = nint((reg%x(n) - reg%itop%x(1)) / pgrid_in%dx) + 1
            if (i3 < 1) i3 = 1
            if (i3 > reg%itop%nnode) i3 = reg%itop%nnode
            do vtype=1,NVTYPE
               if (associated(vgrid(reg%ivgrid,vtype)%velocity)) then
                  reg%velocity(vtype,n) = reg%itop%vel_bot(vtype,i3)
               endif
            end do
         endif
      endif
   endif
end do

return

end subroutine velocities_to_region2d


subroutine tag_regular_nodes2d(pgrid_in)
use mod_fmt2d_nointerfaces
implicit none

type(Tpropagation_grid2d),target :: pgrid_in

integer :: ix,iz,i,j
logical :: is_regular

if (.not.associated(pgrid_in%fully_regular)) then
   allocate(pgrid_in%fully_regular(pgrid_in%nx,pgrid_in%nz))
endif

pgrid_in%fully_regular = .true.

do ix=1,pgrid_in%nx
   do iz=1,pgrid_in%nz

      if (pgrid_in%node_region(ix,iz) <= 0) then
         pgrid_in%fully_regular(ix,iz) = .false.
         cycle
      endif

      is_regular = .true.

      do i=0,1
         do j=0,1
            if ((ix-i > 0) .and. (ix-i < pgrid_in%nx) .and. &
                (iz-j > 0) .and. (iz-j < pgrid_in%nz)) then
               if (associated(pgrid_in%ccind_from_2dc(ix-i,iz-j)%p)) then
                  is_regular = .false.
               endif
            endif
         end do
      end do

      if (ix == 1 .or. ix == pgrid_in%nx) is_regular = .false.
      if (iz == 1 .or. iz == pgrid_in%nz) is_regular = .false.

      pgrid_in%fully_regular(ix,iz) = is_regular

   end do
end do

return

end subroutine tag_regular_nodes2d


subroutine sweep_region_from_interface2d(reg,istart_in,vtype,s)
use mod_fmt2d_nointerfaces
implicit none

type(Tregion2d),target :: reg
type(Tintersection2d),target :: istart_in
integer :: vtype
type(Tsource2d) :: s

integer :: n,m,nt
integer :: istart_nnode
integer,allocatable :: istart_rnode(:)
integer :: connode(20)
integer :: n_connode
integer :: n_activated, total_activated
integer :: ix, iz, ixx, iz_if, iz_start, iz_end, iz_step
integer :: side
logical :: is_bottom, mapped_interface
real(kind=dp) :: dist, v_loc, iface_time, dt, z_iface
real(kind=dp) :: dvec(2), grad_iface(2), norm_iface(2)
type(Tpropagation_grid2d),pointer :: lgrid

lgrid => reg%grid

reg%node_status = -1
reg%arrivaltime = huge_time
reg%time_gradient = 0.0_dp

istart_nnode = istart_in%nnode
allocate(istart_rnode(istart_nnode))
istart_rnode = 0
mapped_interface = .false.

! Use per-region boundary node IDs when this interface is an actual
! boundary of the target region.  Do not reuse mappings that belong to
! another region; combined regions can contain internal interfaces.
if (istart_in%id == reg%itop%id .and. associated(reg%itop_rnode_id)) then
   istart_rnode = reg%itop_rnode_id
   mapped_interface = .true.
else if (istart_in%id == reg%ibot%id .and. associated(reg%ibot_rnode_id)) then
   istart_rnode = reg%ibot_rnode_id
   mapped_interface = .true.
else if (associated(istart_in%regright) .and. associated(istart_in%rright_node_id)) then
   if (istart_in%regright%id == reg%id) then
      istart_rnode = istart_in%rright_node_id
      mapped_interface = .true.
   endif
else if (associated(istart_in%regleft) .and. associated(istart_in%rleft_node_id)) then
   if (istart_in%regleft%id == reg%id) then
      istart_rnode = istart_in%rleft_node_id
      mapped_interface = .true.
   endif
endif

if (.not. mapped_interface) then
   total_activated = 0

   do n=1,istart_nnode
      if (associated(istart_in%arrivaltime)) then
         iface_time = istart_in%arrivaltime(n)
      else
         iface_time = 0.0_dp
      endif
      if (iface_time >= huge_time*0.5_dp) cycle

      if (associated(istart_in%time_gradient)) then
         grad_iface = istart_in%time_gradient(1:2,n)
      else
         grad_iface = 0.0_dp
      endif

      if (associated(istart_in%normal)) then
         norm_iface = istart_in%normal(1:2,n)
      else
         norm_iface(1) = 0.0_dp
         norm_iface(2) = 1.0_dp
      endif

      if (dot_product(grad_iface, norm_iface) < 0.0_dp) then
         side = -1
      else
         side = 1
      endif

      ix = nint((istart_in%x(n) - lgrid%x0) / lgrid%dx) + 1
      if (ix < 1) ix = 1
      if (ix > lgrid%nx) ix = lgrid%nx

      iz_if = 0
      z_iface = istart_in%z(n)
      do iz = 1, lgrid%nz - 1
         if (lgrid%z(iz) <= z_iface .and. z_iface < lgrid%z(iz+1)) then
            iz_if = iz
            exit
         endif
      end do
      if (iz_if == 0) then
         if (z_iface <= lgrid%z(1)) iz_if = 1
         if (z_iface >= lgrid%z(lgrid%nz)) iz_if = lgrid%nz - 1
      endif
      if (iz_if == 0) cycle

      if (side > 0) then
         iz_start = iz_if + 1
         iz_end = min(lgrid%nz, iz_if + 3)
         iz_step = 1
      else
         iz_start = iz_if
         iz_end = max(1, iz_if - 2)
         iz_step = -1
      endif

      do iz = iz_start, iz_end, iz_step
         do ixx = max(1, ix-1), min(lgrid%nx, ix+1)
            nt = reg%region_rnode_id(ixx, iz)
            if (nt <= 0 .or. nt > reg%nnode) cycle
            if (vtype == VTYPE_S .and. reg%velocity(vtype,nt) < 0.01_dp) cycle

            dvec(1) = reg%x(nt) - istart_in%x(n)
            dvec(2) = reg%z(nt) - istart_in%z(n)
            dist = sqrt(sum(dvec**2))
            if (dist <= 1.0e-10_dp) cycle

            v_loc = reg%velocity(vtype, nt)
            if (v_loc < 0.01_dp) v_loc = 0.01_dp

            if (sqrt(sum(grad_iface**2)) > 1.0e-10_dp) then
               dt = dot_product(dvec, grad_iface)
               if (dt <= 0.0_dp) dt = dist / v_loc
            else
               grad_iface = dvec / (dist * v_loc)
               dt = dist / v_loc
            endif

            if (reg%node_status(nt) == -1 .or. iface_time + dt < reg%arrivaltime(nt)) then
               reg%node_status(nt) = 1
               reg%arrivaltime(nt) = iface_time + dt
               reg%time_gradient(1:2,nt) = grad_iface
               total_activated = total_activated + 1
            endif
         end do
      end do
   end do

   deallocate(istart_rnode)
   return
endif

do n=1,istart_nnode
   m = istart_rnode(n)
   if (m > 0 .and. m <= reg%nnode) then
      if (vtype == VTYPE_S .and. reg%velocity(vtype,m) < 0.01_dp) cycle
      reg%node_status(m) = 1
      if (associated(istart_in%arrivaltime)) then
         reg%arrivaltime(m) = istart_in%arrivaltime(n)
      else
         reg%arrivaltime(m) = 0.0_dp
      endif
      if (associated(istart_in%time_gradient)) then
         reg%time_gradient(1:2,m) = istart_in%time_gradient(1:2,n)
      else
         reg%time_gradient(1:2,m) = 0.0_dp
      endif
   endif
end do

total_activated = 0

! Activate grid nodes adjacent to interface nodes.
! For each activated grid node, compute an initial arrival time
! from the nearest interface node using distance/velocity.
! This is critical for bottom interfaces where FMM cannot propagate
! from interface nodes to grid nodes through the normal cell structure.

! Determine if starting interface is itop or ibot
is_bottom = (istart_in%id == reg%ibot%id)

do n=1,istart_nnode
   m = istart_rnode(n)
   if (m <= 0 .or. m > reg%nnode) cycle
   if (reg%node_status(m) /= 1) cycle
   if (reg%node(m)%i1 /= 0) cycle

   ! Get the grid x-index for this interface node
   ix = nint((reg%x(m) - lgrid%x0) / lgrid%dx) + 1
   if (ix < 1 .or. ix > lgrid%nx) ix = nint((reg%x(m) - lgrid%x0) / lgrid%dx) + 1
   if (ix < 1) ix = 1
   if (ix > lgrid%nx) ix = lgrid%nx

   ! Get the grid z-index closest to the interface
   iz_if = 0
   do iz = 1, lgrid%nz - 1
      if (lgrid%z(iz) <= reg%z(m) .and. reg%z(m) < lgrid%z(iz+1)) then
         iz_if = iz
         exit
      endif
   end do
   if (iz_if == 0) then
      if (reg%z(m) <= lgrid%z(1)) iz_if = 1
      if (reg%z(m) >= lgrid%z(lgrid%nz)) iz_if = lgrid%nz - 1
   endif
   if (iz_if == 0) cycle

   ! Search for grid nodes in this region adjacent to the interface
   if (is_bottom) then
      ! Bottom interface: search upward (decreasing iz) into the region
      iz_start = iz_if
      iz_end = max(1, iz_if - 2)
      iz_step = -1
   else
      ! Top interface: search downward (increasing iz) into the region
      iz_start = iz_if + 1
      iz_end = min(lgrid%nz, iz_if + 3)
      iz_step = 1
   endif

   do iz = iz_start, iz_end, iz_step
      do ixx = max(1, ix-1), min(lgrid%nx, ix+1)
         if (reg%region_rnode_id(ixx, iz) > 0) then
            nt = reg%region_rnode_id(ixx, iz)
            if (nt > 0 .and. nt <= reg%nnode) then
               if (reg%node_status(nt) == -1) then
                  if (vtype == VTYPE_S .and. reg%velocity(vtype,nt) < 0.01_dp) cycle
                  reg%node_status(nt) = 1
                  ! Compute initial time from interface node
                  dist = sqrt((reg%x(nt) - reg%x(m))**2 + &
                       (reg%z(nt) - reg%z(m))**2)
                  v_loc = reg%velocity(vtype, nt)
                  if (v_loc < 0.01_dp) v_loc = 0.01_dp
                  reg%arrivaltime(nt) = reg%arrivaltime(m) + dist / v_loc
                  ! Set gradient from interface node
                  reg%time_gradient(1,nt) = (reg%x(nt) - reg%x(m)) / (dist * v_loc)
                  reg%time_gradient(2,nt) = (reg%z(nt) - reg%z(m)) / (dist * v_loc)
                  total_activated = total_activated + 1
               endif
            endif
         endif
      end do
   end do

end do

n = 0
do m=1,reg%nnode
   if (reg%node_status(m) == 1 .and. reg%node(m)%i1 == 0) n = n + 1
end do

deallocate(istart_rnode)

return

contains

subroutine find_connected_nodes2d_local(lreg,centernode,lgrid2,cnodes,ncnodes)
type(Tregion2d) :: lreg
integer :: centernode
type(Tpropagation_grid2d),pointer :: lgrid2
integer :: cnodes(20),ncnodes
integer :: i1,i2,i3,ix,iz,ii,jj,ixx,izz,m,n,icell
type(Tintersection2d),pointer :: isec

i1 = lreg%node(centernode)%i1
i2 = lreg%node(centernode)%i2
i3 = lreg%node(centernode)%i3

ncnodes = 0

if (i1 /= 0) then

   do ix=i1-1,i1+1,2
      if (ix >= 1 .and. ix <= lgrid2%nx) then
         m = lreg%region_rnode_id(ix, i2)
         if (m > 0 .and. m <= lreg%nnode) then
            ncnodes = ncnodes + 1
            cnodes(ncnodes) = m
         endif
      endif
   end do

   do iz=i2-1,i2+1,2
      if (iz >= 1 .and. iz <= lgrid2%nz) then
         m = lreg%region_rnode_id(i1, iz)
         if (m > 0 .and. m <= lreg%nnode) then
            ncnodes = ncnodes + 1
            cnodes(ncnodes) = m
         endif
      endif
   end do

else

   if (i2 == lreg%itop%id) then
      isec => lreg%itop
   else
      isec => lreg%ibot
   endif

   do n=1,4
      icell = isec%ccell_from_inode(n,i3)
      if (icell > 0) then
         ix = isec%ccells(icell)%ix
         iz = isec%ccells(icell)%iz
         do ii=0,1
            do jj=0,1
               ixx = ix + ii
               izz = iz + jj
               if (ixx >= 1 .and. ixx <= lgrid2%nx .and. izz >= 1 .and. izz <= lgrid2%nz) then
                  m = lreg%region_rnode_id(ixx,izz)
                  if (m > 0 .and. m <= lreg%nnode) then
                     ncnodes = ncnodes + 1
                     cnodes(ncnodes) = m
                  endif
               endif
            end do
         end do
      endif
   end do

endif

end subroutine find_connected_nodes2d_local

end subroutine sweep_region_from_interface2d


subroutine initialize_region_from_source2d(reg,sx,sz,t0,vtype_in)
use mod_fmt2d_nointerfaces
implicit none

type(Tregion2d),target :: reg
real(kind=dp) :: sx,sz,t0
integer, optional :: vtype_in

integer :: ix,iz,n,m
integer :: ix0,iz0
integer :: init_vtype
real(kind=dp) :: dist,v_avg
type(Tpropagation_grid2d),pointer :: lgrid

lgrid => reg%grid
init_vtype = VTYPE_P
if (present(vtype_in)) init_vtype = vtype_in

reg%node_status = -1
reg%arrivaltime = huge_time
reg%time_gradient = 0.0_dp

ix0 = floor((sx - lgrid%x0)/lgrid%dx) + 1
iz0 = floor((sz - lgrid%z0)/lgrid%dz) + 1

if (ix0 < 1) ix0 = 1
if (ix0 >= lgrid%nx) ix0 = lgrid%nx - 1
if (iz0 < 1) iz0 = 1
if (iz0 >= lgrid%nz) iz0 = lgrid%nz - 1

do ix=ix0,ix0+1
   do iz=iz0,iz0+1

      if (ix >= 1 .and. ix <= lgrid%nx .and. iz >= 1 .and. iz <= lgrid%nz) then
         if (reg%region_rnode_id(ix,iz) > 0) then

            m = reg%region_rnode_id(ix,iz)
            if (m > 0 .and. m <= reg%nnode) then

               dist = sqrt((lgrid%x(ix)-sx)**2 + (lgrid%z(iz)-sz)**2)
               v_avg = reg%velocity(init_vtype,m)
               if (v_avg < 0.01_dp) v_avg = 0.01_dp

               if (dist < 0.01_dp*lgrid%tolerance) then
                  reg%arrivaltime(m) = t0
                  reg%node_status(m) = 0
                  reg%time_gradient(1:2,m) = 0.0_dp
               else
                  reg%arrivaltime(m) = t0 + dist/v_avg
                  reg%node_status(m) = 1
                  reg%time_gradient(1,m) = (lgrid%x(ix)-sx)/(dist*v_avg)
                  reg%time_gradient(2,m) = (lgrid%z(iz)-sz)/(dist*v_avg)
               endif

            endif

         endif
      endif

   end do
end do

do n=1,reg%nnode
   if (reg%node_status(n) == 0) then
      i1: do ix=-1,1,2
         do iz=-1,1,2
            m = get_neighbor2d(reg,n,ix,iz,lgrid)
            if (m > 0 .and. reg%node_status(m) == -1) then
               reg%node_status(m) = 1
            endif
         end do
      end do i1
   endif
end do

return

contains

function get_neighbor2d(lreg,node,dix,diz,lgrid2) result(nbr)
type(Tregion2d) :: lreg
integer :: node,dix,diz,nbr
type(Tpropagation_grid2d),pointer :: lgrid2
integer :: i1,i2

nbr = 0

i1 = lreg%node(node)%i1
i2 = lreg%node(node)%i2

if (i1 == 0) return

if (i1+dix >= 1 .and. i1+dix <= lgrid2%nx .and. &
    i2+diz >= 1 .and. i2+diz <= lgrid2%nz) then
   if (lreg%region_rnode_id(i1+dix,i2+diz) > 0) then
      nbr = lreg%region_rnode_id(i1+dix,i2+diz)
      if (nbr <= 0 .or. nbr > lreg%nnode) nbr = 0
   endif
endif

return

end function get_neighbor2d

end subroutine initialize_region_from_source2d


subroutine transfer_times_to_intersection(isec, tf)
use mod_fmt2d_nointerfaces
implicit none

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

return

contains

subroutine sample_region_time(lreg, px, pz, atime, grad)
type(Tregion2d) :: lreg
real(kind=dp) :: px, pz, atime, grad(2)
type(Tpropagation_grid2d), pointer :: lgrid
integer :: ix, iz, n1, n2, n3, n4, nn
real(kind=dp) :: w1, w2

lgrid => lreg%grid

ix = floor((px - lgrid%x0) / lgrid%dx) + 1
iz = floor((pz - lgrid%z0) / lgrid%dz) + 1

if (ix < 1) ix = 1
if (ix >= lgrid%nx) ix = lgrid%nx - 1
if (iz < 1) iz = 1
if (iz >= lgrid%nz) iz = lgrid%nz - 1

w1 = (px - lgrid%x(ix)) / lgrid%dx
w2 = (pz - lgrid%z(iz)) / lgrid%dz
if (w1 < 0.0_dp) w1 = 0.0_dp
if (w1 > 1.0_dp) w1 = 1.0_dp
if (w2 < 0.0_dp) w2 = 0.0_dp
if (w2 > 1.0_dp) w2 = 1.0_dp

n1 = lreg%region_rnode_id(ix, iz)
n2 = lreg%region_rnode_id(ix+1, iz)
n3 = lreg%region_rnode_id(ix, iz+1)
n4 = lreg%region_rnode_id(ix+1, iz+1)

if (n1 > 0 .and. n2 > 0 .and. n3 > 0 .and. n4 > 0) then
   atime = (1.0_dp-w1)*(1.0_dp-w2)*lreg%arrivaltime(n1) + &
        w1*(1.0_dp-w2)*lreg%arrivaltime(n2) + &
        (1.0_dp-w1)*w2*lreg%arrivaltime(n3) + &
        w1*w2*lreg%arrivaltime(n4)
   grad = (1.0_dp-w1)*(1.0_dp-w2)*lreg%time_gradient(1:2,n1) + &
        w1*(1.0_dp-w2)*lreg%time_gradient(1:2,n2) + &
        (1.0_dp-w1)*w2*lreg%time_gradient(1:2,n3) + &
        w1*w2*lreg%time_gradient(1:2,n4)
else
   nn = nearest_region_node(lreg, px, pz)
   if (nn > 0) then
      atime = lreg%arrivaltime(nn)
      grad = lreg%time_gradient(1:2, nn)
   else
      atime = huge_time
      grad = 0.0_dp
   endif
endif

end subroutine sample_region_time

function nearest_region_node(lreg, px, pz) result(nn)
type(Tregion2d) :: lreg
real(kind=dp) :: px, pz
integer :: nn, i
real(kind=dp) :: dist, distmin

nn = 0
distmin = 1.0e100_dp
do i = 1, lreg%nnode
   dist = (px - lreg%x(i))**2 + (pz - lreg%z(i))**2
   if (dist < distmin) then
      distmin = dist
      nn = i
   endif
end do

end function nearest_region_node

end subroutine transfer_times_to_intersection


subroutine transfer_diving_return_to_intersection(isec, tf, vtype, min_depth)
use mod_fmt2d_nointerfaces
implicit none

type(Tintersection2d) :: isec
type(Ttime_field2d) :: tf
integer :: vtype
real(kind=dp) :: min_depth

type(Tregion2d), pointer :: lreg
integer :: n, m, best_m
real(kind=dp) :: dx, dz, dist, vloc, cand, best_t, z_if
real(kind=dp) :: grad_norm

if (.not. associated(isec%arrivaltime)) then
   allocate(isec%arrivaltime(isec%nnode))
   allocate(isec%time_gradient(2, isec%nnode))
endif

isec%arrivaltime = huge_time
isec%time_gradient = 0.0_dp

if (.not. associated(tf%reg)) return
lreg => tf%reg
if (.not. associated(tf%arrivaltime)) return

do n = 1, isec%nnode
   z_if = isec%z(n)
   best_t = huge_time
   best_m = 0

   do m = 1, lreg%nnode
      if (m > size(tf%arrivaltime)) cycle
      if (tf%arrivaltime(m) >= huge_time*0.5_dp) cycle
      if (lreg%z(m) <= z_if + min_depth) cycle
      if (abs(lreg%x(m) - isec%x(n)) > 6.0_dp * min_depth) cycle
      if (vtype == VTYPE_S .and. lreg%velocity(vtype,m) < 0.01_dp) cycle

      dx = lreg%x(m) - isec%x(n)
      dz = lreg%z(m) - z_if
      dist = sqrt(dx*dx + dz*dz)
      if (dist <= 1.0e-10_dp) cycle
      vloc = lreg%velocity(vtype,m)
      if (vloc < 0.01_dp) vloc = 0.01_dp
      cand = tf%arrivaltime(m) + dist / vloc
      if (cand < best_t) then
         best_t = cand
         best_m = m
      endif
   end do

   if (best_m > 0) then
      isec%arrivaltime(n) = best_t
      isec%time_gradient(1,n) = (isec%x(n) - lreg%x(best_m))
      isec%time_gradient(2,n) = (z_if - lreg%z(best_m))
      grad_norm = sqrt(sum(isec%time_gradient(1:2,n)**2))
      vloc = lreg%velocity(vtype,best_m)
      if (vloc < 0.01_dp) vloc = 0.01_dp
      if (grad_norm > 1.0e-10_dp) then
         isec%time_gradient(1:2,n) = isec%time_gradient(1:2,n) / (grad_norm * vloc)
      else
         isec%time_gradient(1:2,n) = 0.0_dp
      endif
      if (associated(tf%time_gradient) .and. best_m <= size(tf%time_gradient,2)) then
         isec%time_gradient(1:2,n) = 0.5_dp * (isec%time_gradient(1:2,n) + &
              tf%time_gradient(1:2,best_m))
      endif
   endif
end do

return

end subroutine transfer_diving_return_to_intersection


subroutine reflect_gradient_at_interface(isec, target_reg, target_vtype)
use mod_fmt2d_nointerfaces
implicit none

type(Tintersection2d) :: isec
type(Tregion2d) :: target_reg
integer :: target_vtype
integer :: n, target_side
real(kind=dp) :: norm(2), grad_perp, grad_par(2)
real(kind=dp) :: v_target, p_total, p_par2, p_norm2

do n = 1, isec%nnode
   if (associated(isec%normal)) then
      norm = isec%normal(1:2,n)
   else
      norm(1) = 0.0_dp
      norm(2) = 1.0_dp
   endif

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
      isec%time_gradient(1:2,n) = grad_par + target_side * sqrt(p_norm2) * norm
   endif
end do

return

contains

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

end subroutine reflect_gradient_at_interface


subroutine tag_active_vgrid_nodes2d(reg)
  use mod_fmt2d_nointerfaces
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
