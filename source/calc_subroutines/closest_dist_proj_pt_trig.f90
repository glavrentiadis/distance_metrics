!---------------------------------------------------------------------
!     Projection point on triangle
!     This set of functions and routines returns the coordinates of the
!     projection on a point on a triangle in 3D Euclidean space
! 
!     Author: G. Lavrentiadis
!     Date:   03/13/2019
!     Rev:    0.0
! 
!     Main function: ClosDistTrig3D
!     Dependencies:    ClosDistLineSeg3D
!                      InsideTriag
!                      TrigArea
!                      RotMatBetweenVect
! 
!---------------------------------------------------------------------


      subroutine ClosDistProjPtTrig3D(v1,v2,v3,pt,clos_dist,pt_prj)
!     ClosDistProjPtTrig3D returns the closest distance and projection
!     point on a triangle
!     Input arguments:
!        v1: xyz coordinates of first triangle vertex
!        v2: xyz coordinates of first triangle vertex
!        v3: xyz coordinates of first triangle vertex
!        pt: xyz coordinates of point
!     Output arguments:
!        prj_pt: projection point


!     declare variables
      implicit none
      real, intent(IN), dimension(3) :: v1, v2, v3         !vertices coordinates
      real, intent(IN), dimension(3) :: pt                 !point coordinates
      real, dimension(3) :: v1_pr, v2_pr, v3_pr            !prime vector coordinates
      real, dimension(3) :: pt_pr                          !prime point coordinates
      real, dimension(3) :: orig_off                       !new orign coordinates
      real, dimension(3) :: axis_x, axis_y                 !axes unit vectors
      real, dimension(3) :: v3_pr_yz
      real, dimension(3,3) :: rot_mat1, rot_mat2, rot_mat3 !rotation matrices
      real, dimension(3,3) :: inv_rot_mat                  !inverse rotation matrix
      real :: Atrig, eps
      logical :: flag_inside                               !flag for projection point being in triangle
      real, dimension(3) :: clos_dist2edge                 !closest distance to triangle edges
      integer, dimension(1) :: edge_id                     !edge id with the smallest distance
      real,intent(OUT) :: clos_dist                        !closest distance
      real,intent(OUT), dimension(3) :: pt_prj             !projection point

!     dependencies
      interface
        function RotMatBetweenVects(v_orig, v_rot)
          real, intent(IN), dimension(3) :: v_orig, v_rot
          real, dimension(3,3) :: RotMatBetweenVects
        end function

        function InsideTriag(v1, v2, v3, pt)
          logical :: InsideTriag
          real, intent(IN), dimension(2) :: v1, v2, v3, pt
        end function

        function ClosDistLineSeg3D(v1, v2, pt)
          real, intent(IN), dimension(3) :: v1, v2, pt
          real :: ClosDistLineSeg3D
        end function

        function PrjPtLineSeg3D(v1, v2, pt)
          real, intent(IN), dimension(3) :: v1, v2, pt
          real, dimension(3) :: PrjPtLineSeg3D
        end function

        function TrigArea(v1x, v1y, v2x, v2y, v3x, v3y)
          real, intent(IN) :: v1x, v2x, v3x
          real, intent(IN) :: v1y, v2y, v3y
          real :: TrigArea
        end function

        function matinv3(A)
          real, intent(IN), dimension(3,3) :: A
          real, dimension(3,3) :: matinv3
        end function

        function eye(n)
          integer, intent(IN) :: n
          real, dimension(3,3) :: eye
        end function
      end interface

!     initialize axes vectors
      axis_x = (/ 1., 0., 0. /)
      axis_y = (/ 0., 1., 0. /)
!     initialize projection point
      pt_prj(:) = 0.
!     initialize distance to edges
      clos_dist2edge = (/0., 0., 0./)

!     origin offset
      orig_off = v1

!     compute triangle and point coordinates for new origin
      v1_pr = v1 - orig_off
      v2_pr = v2 - orig_off
      v3_pr = v3 - orig_off
      pt_pr = pt - orig_off 

!     check if points are co-linear
      Atrig = TrigArea(v1_pr(1),v1_pr(2),v2_pr(1),v2_pr(2),v3_pr(1),v3_pr(2))
      eps = 1e-4 !tolerance for zero area
      if (Atrig < eps) then
        flag_inside = .false. !colinear points
        rot_mat1 = eye(3)
        rot_mat3 = eye(3)
      else
!       1st rotation: rotate v1-v2 side to x-x' axis
        rot_mat1 = RotMatBetweenVects(v2_pr,axis_x);            ! rotation matrix
!       compute rotated coordinates
        v1_pr = matmul(rot_mat1, v1_pr)
        v2_pr = matmul(rot_mat1, v2_pr)
        v3_pr = matmul(rot_mat1, v3_pr)
        pt_pr = matmul(rot_mat1, pt_pr)

!       2nd rotation: rotate x-x' axis to v1-v2 side
        rot_mat2 = RotMatBetweenVects(v2_pr,axis_x);
!       vertex 2, yz component
        v3_pr_yz = v3_pr
        v3_pr_yz(1) = 0.
!       rotation matrix
        rot_mat3 = RotMatBetweenVects(v3_pr_yz,axis_y)

!       compute rotated coordinates
        v1_pr = matmul(rot_mat3, v1_pr)
        v2_pr = matmul(rot_mat3, v2_pr)
        v3_pr = matmul(rot_mat3, v3_pr)
        pt_pr = matmul(rot_mat3, pt_pr)

!       check if triangle projection is in triangle
        flag_inside = InsideTriag(v1_pr(1:2),v2_pr(1:2),v3_pr(1:2),pt_pr(1:2))
      end if

!     closest distance and projection point
      if (flag_inside) then
!       closest distance
        clos_dist = abs(pt_pr(3))
!       projection point
        pt_prj(1:2) = pt_pr(1:2)
      else
        clos_dist2edge(1) = ClosDistLineSeg3D(v1_pr,v2_pr,pt_pr)
        clos_dist2edge(2) = ClosDistLineSeg3D(v2_pr,v3_pr,pt_pr)
        clos_dist2edge(3) = ClosDistLineSeg3D(v1_pr,v3_pr,pt_pr)
!       closest distance
        clos_dist = minval(clos_dist2edge)
!       projection point
        edge_id = minloc(clos_dist2edge)
        if (edge_id(1) == 1) then
          pt_prj = PrjPtLineSeg3D(v1_pr,v2_pr,pt_pr)
        elseif (edge_id(1) == 2) then
          pt_prj = PrjPtLineSeg3D(v2_pr,v3_pr,pt_pr)
        else
          pt_prj = PrjPtLineSeg3D(v1_pr,v3_pr,pt_pr)
        end if
      end if

!     remove effects of rotation
!     3rd rotation
      inv_rot_mat = matinv3(rot_mat3)
      pt_prj = matmul(inv_rot_mat,pt_prj)
!     1st rotation
      inv_rot_mat = matinv3(rot_mat1)
      pt_prj = matmul(inv_rot_mat,pt_prj)

!     remove effect of offset
      pt_prj = pt_prj + orig_off 

      end subroutine ClosDistProjPtTrig3D


!--------------------------------------------------------------
!     Closest distance between a point and line segment

      function PrjPtLineSeg3D(v1, v2, pt)
!     PrjPtLineSeg3D retruns the prjection point of pt on a line segment v1-v2
!     Input arguments:
!        v1: xyz coordinates of start of line segment
!        v2: xyz coordinates of end of line segment
!     Output arguments:
!        pt_prj: projection point

!     declare variables
      implicit none
      real, intent(IN), dimension(3) :: v1, v2    !vertices coordinates
      real, intent(IN), dimension(3) :: pt        !point coordinates
      real, dimension(3) :: v1_pr, v2_pr, pt_pr   !prime point and vertex coordinates
      real, dimension(3) :: norm_v1v2             !unit vect v1-v2
      real :: len_seg, len_prj, eps               !v1-v2 segment len, length v1-pt_prj 
      real, dimension(3) :: pt_prj                !proj. point of pt  in v1-v2
      real, dimension(3) :: PrjPtLineSeg3D
      
      eps = 1e-6 !tolerance for zero length

!     set orign to v1
      pt_pr = pt - v1
      v2_pr = v2 - v1
      v1_pr = v1 - v1

!     length of line segment v1-v2
      len_seg = sqrt(sum( (v1_pr-v2_pr)**2 ))

      if (len_seg < eps) then
!       v1 and v2 points are collocated
        pt_prj = v1_pr
      else
!       unit vector v1-v2
        norm_v1v2 = v2_pr/len_seg

!       projection point on line segment
        len_prj = dot_product(norm_v1v2,pt_pr)
!       limit projection length between 0 (projection point is v1) and
!       len_seg (projection point is v2)
        len_prj = max(len_prj,0.0)
        len_prj = min(len_prj,len_seg)

!       projection point coordinates
        pt_prj = len_prj * norm_v1v2

      end if

      !correct for coordinate offset
      pt_prj = pt_prj + v1
      PrjPtLineSeg3D = pt_prj

      end function PrjPtLineSeg3D



