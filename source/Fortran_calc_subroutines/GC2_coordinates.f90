!---------------------------------------------------------------------
!     GC2 Coordinate System
!     This set of function computes the GC2 distances
!
!     Main function:	GC2Coordinates
!     Dependencies:	GC2CoordinatesPoint
!			GC2DiscordantSeg
!
!---------------------------------------------------------------------

      function GC2Coordinates(fault_xyz,n_seg,npt_seg,pt_xyz,n_pt)
!     GC2Coordinates returns the U and T distances based on the 
!     GC2 coordinate system
!     Input arguments:
!        fault_xyz: fault coordinates, 3rd dimension for different fault branches
!        n_seg: number of segments
!        npt_seg: number of points per segment
!        pt_xyz: point-station coordinates
!        n_pt: number of points
!     Output arguments:
!        t_pt: fault normal coordinate of points, GC2 coordinate system
!        u_pt: fault parallel coordinate of points, GC2 coordinate system

!     declare variables
      use memory_module
      implicit none
      real, intent(IN), dimension(max_num_pt_seg,3,max_num_seg) :: fault_xyz !fault trace coordinates
      real, dimension(max_num_pt_seg,3,max_num_seg) :: fault_cor_xyz         !fault trace, corrected
      integer, intent(IN) :: n_seg                              !number of segments
      integer, intent(IN), dimension(max_num_seg) :: npt_seg    !number of points per segment
      real, intent(IN), dimension(max_num_sta,3) :: pt_xyz      !staions'/points' coordinates
      integer, intent(IN) :: n_pt                               !number of stations/points
      real, dimension(2) :: pt_tu                               !fualt norm. and par. pt coordinates
      real, dimension(max_num_sta,2) :: GC2Coordinates          !fualt normal and parallel coordinates
      real, dimension(3) :: orig_fault_xyz                      !fault origin
      real, dimension(3) :: b_hat                               !general fault parallel vector
      real, dimension(2) :: pt_xy
      integer :: j                                              !loop index

      interface
        function GC2CoorPoint(fault_xyz,orig_fault_xyz,b_hat,n_seg,npt_seg,pt_xy)
          use memory_module
          real, intent(IN), dimension(max_num_pt_seg,3,max_num_seg) :: fault_xyz
          real, intent(IN), dimension(3) :: orig_fault_xyz
          real, intent(IN), dimension(2) :: b_hat
          integer, intent(IN) :: n_seg
          integer, intent(IN), dimension(max_num_seg) :: npt_seg
          real, intent(IN), dimension(2) :: pt_xy
          real, dimension(2) :: GC2CoorPoint
        end function
      
        subroutine GC2DiscordantSeg(fault_xyz,n_seg,npt_seg,orig_fault_xyz,b_hat)
          use memory_module
          real, intent(INOUT), dimension(max_num_pt_seg,3,max_num_seg) :: fault_xyz
          integer, intent(IN) :: n_seg
          integer, intent(IN), dimension(max_num_seg) :: npt_seg
          real, intent(OUT), dimension(3) :: orig_fault_xyz
          real, intent(OUT), dimension(2) :: b_hat
        end subroutine
      end interface

!     revert discordant segments
      fault_cor_xyz = fault_xyz
      call GC2DiscordantSeg(fault_cor_xyz,n_seg,npt_seg,orig_fault_xyz,b_hat)

!     compute coordinates of all points
      do j = 1,n_pt
         pt_xy = pt_xyz(j,1:2)
         pt_tu = GC2CoorPoint(fault_cor_xyz,orig_fault_xyz,b_hat,n_seg,npt_seg,pt_xy)
         GC2Coordinates(j,:) = pt_tu
      end do 

      end function GC2Coordinates


!--------------------------------------------------------------
!    Computes the U and T coordinates for a point based on the GC2 system

      function GC2CoorPoint(fault_xyz,orig_fault_xyz,b_hat,n_seg,npt_seg,pt_xy)
!     GC2CoordinatesPoint returns the GC2 coordinate system for the point pt_xy
!     Input arguments:
!        fault_xyz: fault coordinates, 3rd dimension for different fault branches
!        orig_fault_xyz: fault origin coordinates
!        b_hat: general fault parallel vector
!        n_seg: number of segments
!        npt_seg: number of points per segment
!        pt_xy: point-station coordinates

!     declare variables
      use memory_module
      implicit none
      real, intent(IN), dimension(max_num_pt_seg,3,max_num_seg) :: fault_xyz !fault trace coordinates
      real, intent(IN), dimension(3) :: orig_fault_xyz       !fault origin
      real, intent(IN), dimension(2) :: b_hat                !fault parallel coordinate
      integer, intent(IN) :: n_seg                           !number of segments
      integer, intent(IN), dimension(max_num_seg) :: npt_seg !number of points per segment
      real, intent(IN), dimension(2) :: pt_xy                !point-station coordinates
      real :: u_shift                                        !segment parallel coordinates
      real, dimension(2) :: pt_rel_xy                        !rel. pt. coor. with respect to seg. j
      real, dimension(2) :: u_vec, t_vec                     !unit parallel and normal to seg. vector
      real, dimension(max_num_pt_seg,max_num_seg) :: l_seg, s_seg !length segment, cumulative dist.
      real, dimension(max_num_pt_seg,max_num_seg) :: u_seg, t_seg !fault normal and par. coordinates
      real, dimension(max_num_pt_seg,max_num_seg) :: w_seg        !segment weights
      real :: w_tot                                          !total weights
      logical :: flag_exit                                   !flag to exit do loop
      integer :: j, k                                        !loop index
      real :: eps                                            !tolerance
      real :: u_pt, t_pt                                     !u & t GC2 coordinates
      real, dimension(2) :: GC2CoorPoint                     !ut GC2 coordinates

      eps = 1e-8;


      do j = 1,n_seg
!       segment parallel coordinate shift
        u_shift = dot_product(fault_xyz(1,1:2,j)-orig_fault_xyz(1:2),b_hat)
        do k = 1,npt_seg(j)-1
!         relative point coordinates
          pt_rel_xy = pt_xy - fault_xyz(k,1:2,j)
!         segment parallel vector
          u_vec = fault_xyz(k+1,1:2,j) - fault_xyz(k,1:2,j)
!         segment length
          l_seg(k,j) = sqrt(sum(u_vec**2))
!         cumulative length + coordinate shift
          s_seg(k,j) = sum(l_seg(1:k-1,j)) + u_shift
!         unit-vector in the segment parallel direction
          u_vec = u_vec/l_seg(k,j)
!         unit-vector segment normal 
          t_vec(1) = u_vec(2)
          t_vec(2) = -1*u_vec(1)
!         fault normal and parallel coordinates
          u_seg(k,j) = dot_product(pt_rel_xy,u_vec)
          t_seg(k,j) = dot_product(pt_rel_xy,t_vec)
        end do
      end do

!     compute segment weights
      flag_exit = .true.
      w_seg(:,:) = 0.0


     do j = 1,n_seg
       do k = 1,npt_seg(j)-1
!        segment weights
         if (abs(t_seg(k,j)) > eps) then
!          point non-colinear with segment
           w_seg(k,j) = 1/t_seg(k,j)*(  atan((l_seg(k,j)-u_seg(k,j))/t_seg(k,j)) &
                                      - atan(-u_seg(k,j)/t_seg(k,j)))
         else if (or(u_seg(k,j) < 0, u_seg(k,j) > l_seg(k,j))) then
!          point on the extension of segment
           w_seg(k,j) = 1/(u_seg(k,j)-l_seg(k,j)) - 1/u_seg(k,j)
         else
!          point is segment j
           w_seg(:,:) = 0
           w_seg(k,j) = 1
           flag_exit = .false.
           exit
         end if
       end do
!      flag exit
       if (.not. flag_exit) then 
         exit
       end if
     end do
!    total weight
     w_tot = sum(w_seg)
!    normalize weights
     w_seg = w_seg / w_tot

!    parallel and normal point coordinate
     t_pt = 0.0
     u_pt = 0.0
     do j = 1,n_seg
       k = npt_seg(j)-1
!      normal fault coordinate
       t_pt = t_pt + dot_product(t_seg(1:k,j),w_seg(1:k,j))
!      parallel fault coordinate
       u_pt = u_pt + dot_product(u_seg(1:k,j)+s_seg(1:k,j),w_seg(1:k,j))
     end do

     GC2CoorPoint(1) = t_pt
     GC2CoorPoint(2) = u_pt
     end function GC2CoorPoint


!--------------------------------------------------------------
!     Identifies and rotates the discordant segments 

      subroutine GC2DiscordantSeg(fault_xyz,n_seg,npt_seg,orig_fault_xyz,b_hat)
!     GC2DiscordantSeg identifies and rotates the discordant segments in the GC2
!     coordinate stystem
!     Input arguments:
!        fault_xyz: fault coordinates, 3rd dimension for different fault branches
!        n_seg: number of segments
!        npt_seg: number of points per segment
!     Output argument:
!        fault_xzy: fault coordinates corrected discordant segments 
!        orig_fault_xyz: fault origin coordinates
!        b_hat: general fault parallel vector

!     declare variables
      use memory_module
      implicit none
      real, intent(INOUT), dimension(max_num_pt_seg,3,max_num_seg) :: fault_xyz !fault trace coor.
      integer, intent(IN) :: n_seg                            !number of segments
      integer, intent(IN), dimension(max_num_seg) :: npt_seg  !number of points per segment
      real, intent(OUT), dimension(3) :: orig_fault_xyz       !fault origin
      real, intent(OUT), dimension(2) :: b_hat                !fault parallel coordinate
      real, dimension(max_num_seg,max_num_seg) :: dist_seg2seg!distance between segments
      integer, dimension(2) :: i_maxdist                      !indices max distance
      real, dimension(2) :: a_hat                             !vect connect. the two most distant seg.
      real, dimension(2) :: v_seg                             !horizontal component of segment
      real, dimension(max_num_pt_seg,3) :: seg_xyz            !segment coordinates
      real, dimension(max_num_seg) :: e_j                     !component of v_seg in a_hat direction
      real :: E                                               !sum of e_j of all segments
      integer :: j, k                                         !loop index
      real :: eps                                             !tolerance

!     epsilon tolerance
      eps = 1e-5

!      find the most distant segments
      do j = 1,n_seg
        do k = j,n_seg
          dist_seg2seg(j,k) = sqrt(  (fault_xyz(1,1,j) - fault_xyz(1,1,k))**2 &
                                   + (fault_xyz(1,2,j) - fault_xyz(1,2,k))**2 )
        end do
      end do
      i_maxdist = maxloc(dist_seg2seg)

!     define alpha vector connecting end segments
      a_hat(1) = fault_xyz(1,1,i_maxdist(2)) - fault_xyz(1,1, i_maxdist(1))
      a_hat(2) = fault_xyz(1,2,i_maxdist(2)) - fault_xyz(1,2, i_maxdist(1))
      a_hat = a_hat/sqrt(sum(a_hat**2)) !unit-vector

!     compute E to identify discordant segments
      do j = 1,n_seg
!       vector connecting the end points of seg j
        v_seg = fault_xyz(npt_seg(j),1:2,j) - fault_xyz(1,1:2,j)
        e_j(j) = dot_product(v_seg,a_hat); 
      end do
      E = sum(e_j)

!     reverse discordant segments
      do j = 1,n_seg
!       coordinates of segment j
        seg_xyz = fault_xyz(:,:,j);
!       check if polarity should be reversed
        if (e_j(j)*E < 0) then 
!         reverse segment
          fault_xyz(1:npt_seg(j),:,j) = seg_xyz(npt_seg(j):-1:1,:);
        end if
      end do

!     find new most distant segments
      do j = 1,n_seg
        do k = j,n_seg
          dist_seg2seg(j,k) = sqrt(  (fault_xyz(1,1,j) - fault_xyz(1,1,k))**2 &
                                   + (fault_xyz(1,2,j) - fault_xyz(1,2,k))**2 )
        end do
      end do
      i_maxdist = maxloc(dist_seg2seg)
!     fault origin point
      orig_fault_xyz = fault_xyz(1,:,i_maxdist(1))
 
!     define beta vector, average strike
      b_hat(:) = 0.
      do j = 1,n_seg
        !vector connecting the end points of seg j
        v_seg = fault_xyz(npt_seg(j),1:2,j) - fault_xyz(1,1:2,j)
        b_hat = b_hat + v_seg
      end do
      b_hat = b_hat/sqrt(sum(b_hat**2)) !unit-vector

      end subroutine GC2DiscordantSeg
