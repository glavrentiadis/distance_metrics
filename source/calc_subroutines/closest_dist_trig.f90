!---------------------------------------------------------------------
!     Closest distance to triangle
!     This set of functions and routines returns the closest distance
!     between a triangle and point in 3D Euclidean space
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


      function ClosDistTrig3D(v1,v2,v3,pt)
!     ClosDistTrig3D returns the closest distance to a triangle
!     Input arguments:
!        v1: xyz coordinates of first triangle vertex
!        v2: xyz coordinates of first triangle vertex
!        v3: xyz coordinates of first triangle vertex
!        pt: xyz coordinates of point
!     Output arguments:
!        clos_dist: closest distance


!     declare variables
      implicit none
      real, intent(IN), dimension(3) :: v1, v2, v3         !vertices coordinates
      real, intent(IN), dimension(3) :: pt                 !point coordinates
      real, dimension(3) :: v1_pr, v2_pr, v3_pr            !prime vector coordinates
      real, dimension(3) :: pt_pr                          !prime point coordinates
      real, dimension(3) :: orig_off                       !new orign coordinates
      real, dimension(3) :: axis_x, axis_y                 !axes unit vectors
      real, dimension(3) :: v3_pr_yz
      real, dimension(3,3) :: rot_mat                      !rotation matrix
      real :: Atrig, eps
      logical :: flag_inside                               !flag for projection point being in triangle
      real :: clos_dist1, clos_dist2, clos_dist3           !closest distance to triangle edges
      real :: clos_dist                                    !closest distance
      real :: ClosDistTrig3D                               !closest distance, output

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

        function TrigArea(v1x, v1y, v2x, v2y, v3x, v3y)
          real, intent(IN) :: v1x, v2x, v3x
          real, intent(IN) :: v1y, v2y, v3y
          real :: TrigArea
        end function

      end interface

!     initialize axes vectors
      axis_x = (/ 1., 0., 0. /)
      axis_y = (/ 0., 1., 0. /)

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
      else
!       1st rotation: rotate v1-v2 side to x-x' axis
        rot_mat = RotMatBetweenVects(v2_pr,axis_x);            ! rotation matrix
!       compute rotated coordinates
        v1_pr = matmul(rot_mat, v1_pr)
        v2_pr = matmul(rot_mat, v2_pr)
        v3_pr = matmul(rot_mat, v3_pr)
        pt_pr = matmul(rot_mat, pt_pr)

!       2nd rotation: rotate x-x' axis to v1-v2 side
        rot_mat = RotMatBetweenVects(v2_pr,axis_x);
!       vertex 2, yz component
        v3_pr_yz = v3_pr
        v3_pr_yz(1) = 0.
!       rotation matrix
        rot_mat = RotMatBetweenVects(v3_pr_yz,axis_y)

!       compute rotated coordinates
        v1_pr = matmul(rot_mat, v1_pr)
        v2_pr = matmul(rot_mat, v2_pr)
        v3_pr = matmul(rot_mat, v3_pr)
        pt_pr = matmul(rot_mat, pt_pr)

!       check if triangle projection is in triangle
        flag_inside = InsideTriag(v1_pr(1:2),v2_pr(1:2),v3_pr(1:2),pt_pr(1:2))
      end if

!       determine closest distance
      if (flag_inside) then 
        clos_dist = abs(pt_pr(3))
      else
        clos_dist1 = ClosDistLineSeg3D(v1_pr,v2_pr,pt_pr)
        clos_dist2 = ClosDistLineSeg3D(v2_pr,v3_pr,pt_pr)
        clos_dist3 = ClosDistLineSeg3D(v1_pr,v3_pr,pt_pr)
!       closest distance equal to the minimum side distance
        clos_dist = min(min(clos_dist1,clos_dist2),clos_dist3);
      end if
      ClosDistTrig3D = clos_dist

      end function ClosDistTrig3D


!--------------------------------------------------------------
!     Closest distance between a point and line segment

      function ClosDistLineSeg3D(v1, v2, pt)
!     ClosDistLineSeg retruns the closest distance between a line segment v1-v2 and pt
!     Input arguments:
!        v1: xyz coordinates of start of line segment
!        v2: xyz coordinates of end of line segment
!     Output arguments:
!        clos_dist: closest distance

!     declare variables
      implicit none
      real, intent(IN), dimension(3) :: v1, v2    !vertices coordinates
      real, intent(IN), dimension(3) :: pt        !point coordinates
      real, dimension(3) :: v1_pr, v2_pr, pt_pr   !prime point and vertex coordinates
      real, dimension(3) :: norm_v1v2             !unit vect v1-v2
      real, dimension(3) :: pt_prj                !proj. point of pt  in v1-v2
      real :: len_seg, len_prj, eps               !v1-v2 segment len, length v1-pt_prj 
      real :: clos_dist, ClosDistLineSeg3D        !closest distance output

      eps = 1e-6 !tolerance for zero length

!     set orign to v1
      pt_pr = pt - v1
      v2_pr = v2 - v1
      v1_pr = v1 - v1

!     length of line segment v1-v2
      len_seg = sqrt(sum( (v1_pr-v2_pr)**2 ))

      if (len_seg < eps) then
!       v1 and v2 points are collocated
        clos_dist = sqrt(sum( (pt_pr - v1_pr)**2 ))
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

!       closest distance
        clos_dist = sqrt(sum( (pt_pr - pt_prj)**2 ))
      end if

      ClosDistLineSeg3D = clos_dist

      end function ClosDistLineSeg3D


!--------------------------------------------------------------
!     Creates rotation matrix between to vectors

      function RotMatBetweenVects(a,b)
!     returns the matrix that applies the rotation of vector a to b (b = R * a)
!     Input arguments:
!        a: original vector
!        b: rotated vector
!     Output arguments:
!        rot_mat: rotation matrix

!     declare variables
      implicit none
      real, intent(IN), dimension(3) :: a        !a is the original vector
      real, intent(IN), dimension(3) :: b        !b is the final-rotated vector
      real :: norm_a, norm_b                     !norm of a and b vectors
      real, dimension(3) :: n_a, n_b             !normalized a and b vectors
      real, dimension (3) :: v                   !rotation axis
      real :: c                                  !cosine of rotation angle
      real, dimension(3,3) :: skp_v, skp_v2      !skewed-sym cross prod. and squared of it
      real, dimension(3,3) :: rot_mat            !rotation matrix
      real, dimension(3,3) :: RotMatBetweenVects !rotation matrix, output
      integer :: k                               !loop index


!     compute a and b unit vectors
      norm_a = sqrt(sum( a**2 ))
      norm_b = sqrt(sum( b**2 ))

      n_a = a/norm_a
      n_b = b/norm_b

!     rotation axis, cross product
      v(1) = n_a(2)*n_b(3) - n_a(3)*n_b(2) 
      v(2) = n_a(3)*n_b(1) - n_a(1)*n_b(3)
      v(3) = n_a(1)*n_b(2) - n_a(2)*n_b(1)

!     cosine of rotation angle
      c = dot_product(n_a,n_b)

!     skewed-symmetric cross product matrix
!     first line
      skp_v(1,1) =  0.0
      skp_v(1,2) = -v(3)
      skp_v(1,3) =  v(2)
!     second line
      skp_v(2,1) =  v(3)
      skp_v(2,2) =  0.0
      skp_v(2,3) = -v(1)
!     thrid line
      skp_v(3,1) = -v(2)
      skp_v(3,2) =  v(1)
      skp_v(3,3) =  0.0

!     squared skewed-symmetric cross product matrix
      skp_v2 = matmul(skp_v,skp_v)

!     rotation matrix
!     equation: rot_mat = eye(3) + skp_v + 1./(1+c) * skp_v^2
      rot_mat = skp_v + 1./(1.+c) * skp_v2
      do k = 1,3
         rot_mat(k,k) = rot_mat(k,k) + 1. 
      end do

!     use identify matrix if vectors a and b are parallel
      if (abs(1+c) < 1e-5) then
        rot_mat(:,:) = 0.
        do k = 1,3
          rot_mat(k,k) = 1.
        end do
      end if 

      RotMatBetweenVects = rot_mat

      end function RotMatBetweenVects


!--------------------------------------------------------------
!     Area of triangle computed from it vertex coordinates

      function TrigArea(x1,y1,x2,y2,x3,y3)
!     TrigArea returns the area of a triangle based on the coordinates of the verticies
!     input arguments:
!        x1: 1st vertex, x-coordinate
!        y1: 1st vertex, y-coordinate
!        x2: 2nd vertex, x-coordinate
!        y2: 2nd vertex, y-coordinate
!        x3: 3rd vertex, x-coordinate
!        y3: 3rd vertex, y-coordinate
!     output arguments
!        A: area of triangle

!     declare variables
      implicit none
      real, intent(IN) :: x1, y1, x2, y2, x3, y3 !vertex coordinates 
      real :: A                                  !triangle area
      real :: TrigArea

      A = abs((x1 * (y2 - y3) + x2 * (y3 - y1)  + x3 * (y1 - y2)) / 2.0)
      TrigArea = A

      end function TrigArea


!--------------------------------------------------------------
!     Checks if a point is inside a triangle in 2D space

      function InsideTriag(v1, v2, v3, pt)
!     InsideTriag determines if pt is inside a triangle defined by the v1, v2 and v3 verticies

!     declare variables
      implicit none
      real, intent(IN), dimension(2) :: v1, v2, v3      !coordinates of v1, v2 and v3 verticies
      real, intent(IN), dimension(2) :: pt              !point coordinates
      real :: A, Asub1, Asub2, Asub3, Atot_sub          !Triangle area, area of sub-triangles 
      real :: eps                                       !tolerance in the area comparison
      logical :: flag_pt, InsideTriag                   !flag is true if point is inside the triangle

      interface
        function TrigArea(v1x, v1y, v2x, v2y, v3x, v3y)
          real, intent(IN) :: v1x, v2x, v3x
          real, intent(IN) :: v1y, v2y, v3y
          real :: TrigArea
        end function
      end interface

!     initialize parameters
      flag_pt = .false.
      eps = 1e-4 !tolerance in area calculation 

!     A: triangle
      A = TrigArea(v1(1),v1(2),v2(1),v2(2),v3(1),v3(2))

!     A1: area between v1, v2 and pt
      Asub1 = TrigArea(v1(1),v1(2),v2(1),v2(2),pt(1),pt(2))
!     A2: area between v2, v3 and pt
      Asub2 = TrigArea(v2(1),v2(2),v3(1),v3(2),pt(1),pt(2))
!     A3: area between v1, v3 and pt
      Asub3 = TrigArea(v1(1),v1(2),v3(1),v3(2),pt(1),pt(2))

!     total area of 3 sub-triangles
      Atot_sub = Asub1 + Asub2 + Asub3;

!     check if pt is inside the v1, v2 and v3 triangle
      if (abs(Atot_sub - A) / A  < eps) then
         flag_pt = .true.
      end if
      InsideTriag = flag_pt

      end function InsideTriag
