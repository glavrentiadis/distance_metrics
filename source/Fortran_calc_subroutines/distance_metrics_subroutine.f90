!------------------------------------------------------------------
!     Distance Metrics Main
!     It is a fortran subroutine that computes common distance metrics used in GMPES
! 
!-------------------------------------------------------------------

! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
!     Module to define the size of the different arrays

      module memory_module
!       memory_module declares the size of the different arrays

!       fault arrays
        integer, parameter :: max_num_seg = 1        !maximum number of segments
        integer, parameter :: max_num_pt_seg = 100000!maximum number of points per segment
!       station arrays
        integer, parameter :: max_num_sta = 200000   !maximum number of stations

      end module memory_module


!-------------------------------------------------------------!
! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- !
! -   Main subroutine for calculating distance metrics      - ! 
! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- !
!-------------------------------------------------------------!

      subroutine dist_metrics(num_seg,num_pt_seg,flt_cor_top,flt_cor_base,num_sta,sta_cor, &
                             & r_rup,r_jb,r_x,r_y,r_y0,u_pt,t_pt)

!     declare variables
      use memory_module
      implicit none
      integer :: i, k                                           !loop indices
      integer, intent(IN) :: num_seg, num_sta                   !num of segmetns, num of stations
      integer, intent(IN), dimension(max_num_seg) :: num_pt_seg !num of points per segment
      !coor. at the top of the fault, 3rd dim for diff. seg.
      real, intent(IN), dimension(max_num_pt_seg,3,max_num_seg) :: flt_cor_top
      !coor. at the base of the fault, 3rd dim for diff. seg.
      real, intent(IN), dimension(max_num_pt_seg,3,max_num_seg) :: flt_cor_base 
      real, dimension(max_num_sta,3) :: flt_seg_k_cor           !segment coordinates
      real, intent(IN), dimension(max_num_sta,3) :: sta_cor     !stations' coordinates
!-----GC2 station coordinates
      real, dimension(max_num_sta,2) :: tu_pt                   !GC2 station coordinates
!-----GC2 fault coordinates
      real, dimension(max_num_pt_seg,2,max_num_seg) :: tu_fault 
      real, dimension(max_num_sta,2) :: tu_flt_seg_k            !segment coordinates
      real :: len_fault                                         !fault length
!-----distance metrics arrays
      real, dimension(max_num_sta) :: r_rup_k, r_jb_k           !segment distances
      real, intent(OUT), dimension(max_num_sta) :: r_rup, r_jb
      real, intent(OUT), dimension(max_num_sta) :: r_x, r_y, r_y0
      real, intent(OUT), dimension(max_num_sta) :: u_pt, t_pt

      interface
         function FaultRrup(fault_top, fualt_base, sta, n_flt, n_sta)
           use memory_module
           integer, intent(IN) :: n_flt, n_sta
           real, intent(IN), dimension(max_num_pt_seg,3) :: fault_top, fualt_base
           real, intent(IN), dimension(max_num_sta,3) :: sta
           real, dimension(max_num_sta) :: FaultRrup
         end function

         function FaultRJB(fault_top, fualt_base, sta, n_flt, n_sta)
           use memory_module
           integer, intent(IN) :: n_flt, n_sta
           real, intent(IN), dimension(max_num_pt_seg,3) :: fault_top, fualt_base
           real, intent(IN), dimension(max_num_sta,3) :: sta
           real, dimension(max_num_sta) :: FaultRJB
         end function

         function GC2Coordinates(fault_xyz,n_seg,npt_seg,pt_xyz,n_pt)
           use memory_module
           real, intent(IN), dimension(max_num_pt_seg,3,max_num_seg) :: fault_xyz
           integer, intent(IN) :: n_seg
           integer, intent(IN), dimension(max_num_seg) :: npt_seg
           real, intent(IN), dimension(max_num_sta,3) :: pt_xyz
           integer, intent(IN) :: n_pt
           real, dimension(max_num_sta,2) :: GC2Coordinates
         end function
      end interface 

!     Compute distance metrics
! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
!     compute closest distance (R_rup)
      r_rup = FaultRrup(flt_cor_top(:,:,1),flt_cor_base(:,:,1),sta_cor,num_pt_seg(1),num_sta)
!     compute Joyner-Boore distance (R_jb)
      r_jb  =  FaultRJB(flt_cor_top(:,:,1),flt_cor_base(:,:,1),sta_cor,num_pt_seg(1),num_sta)
      do k = 2,num_seg
!       compute closest distance (R_rup)
        r_rup_k = FaultRrup(flt_cor_top(:,:,k),flt_cor_base(:,:,k),sta_cor, &
                            & num_pt_seg(k),num_sta)
        r_rup = min(r_rup,r_rup_k)

!       compute Joyner-Boore distance (R_jb)
        r_jb_k  = FaultRJB(flt_cor_top(:,:,k),flt_cor_base(:,:,k),sta_cor, &
                           & num_pt_seg(k),num_sta)
        r_jb = min(r_jb,r_jb_k)
      end do

!     compute GC2 station coordinates
      tu_pt = GC2Coordinates(flt_cor_top,num_seg,num_pt_seg,sta_cor,num_sta)
      u_pt = tu_pt(:,2)
      t_pt = tu_pt(:,1)
!     compute GC2 fault coordinates
      do k = 1,num_seg
!       copy the coordinates of the k^th segment in flt_seg_k_cor       
        do i = 1,num_pt_seg(k)
          flt_seg_k_cor(i,:) = flt_cor_top(i,:,k)
        end do
!       compute the GC2 coordinates of the k^th segment
        tu_flt_seg_k = GC2Coordinates(flt_cor_top,num_seg,num_pt_seg, &
                                         & flt_cor_top(:,:,k),num_pt_seg(k))
!       copy the GC2 coordinates of the k^th array in tu_fualt
        do i = 1,num_pt_seg(k)
          tu_fault(i,:,k) = tu_flt_seg_k(i,:)
        end do
      end do
      len_fault = maxval(tu_fault(:,2,:))
   
!     compute  distance metrics (Rx, Ry, R_y0)
      r_x = t_pt
      !parallel to strike distance from the center of rupture
      r_y = abs(u_pt - len_fault/2.)
      !parallel to strike distance from the end of rupture
      do k = 1,num_sta
        if ((u_pt(k) >= 0.) .and. (u_pt(k) <= len_fault)) then
!         point is located in the rupture projection
          r_y0(k) = 0.
        else
          r_y0(k) = min(abs(u_pt(k)), abs(u_pt(k) - len_fault)) 
        end if
      end do

      end subroutine dist_metrics


!---------------------------------------------------------------------------------!
! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- !
! -   Main subroutine for calculating distance metrics and projection points    - ! 
! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- !
!---------------------------------------------------------------------------------!

      subroutine dist_metrics_and_prjpt(num_seg,num_pt_seg,flt_cor_top,flt_cor_base, & 
                                      & num_sta,sta_cor, &
                                      & r_rup,r_jb,r_x,r_y,r_y0,u_pt,t_pt,pt_prj)

!     declare variables
      use memory_module
      implicit none
      integer :: i, k                                           !loop indices
      integer, intent(IN) :: num_seg, num_sta                   !num of segmetns, num of stations
      integer, intent(IN), dimension(max_num_seg) :: num_pt_seg !num of points per segment
      !coor. at the top of the fault, 3rd dim for diff. seg.
      real, intent(IN), dimension(max_num_pt_seg,3,max_num_seg) :: flt_cor_top
      !coor. at the base of the fault, 3rd dim for diff. seg.
      real, intent(IN), dimension(max_num_pt_seg,3,max_num_seg) :: flt_cor_base 
      real, dimension(max_num_sta,3) :: flt_seg_k_cor           !segment coordinates
      real, intent(IN), dimension(max_num_sta,3) :: sta_cor     !stations' coordinates
!-----GC2 station coordinates
      real, dimension(max_num_sta,2) :: tu_pt                   !GC2 station coordinates
!-----GC2 fault coordinates
      real, dimension(max_num_pt_seg,2,max_num_seg) :: tu_fault 
      real, dimension(max_num_sta,2) :: tu_flt_seg_k            !segment coordinates
      real :: len_fault                                         !fault length
!-----distance metrics arrays
      real, dimension(max_num_sta) :: r_rup_k, r_jb_k           !segment distances
      real, intent(OUT), dimension(max_num_sta) :: r_rup, r_jb
      real, intent(OUT), dimension(max_num_sta) :: r_x, r_y, r_y0
      real, intent(OUT), dimension(max_num_sta) :: u_pt, t_pt
!-----projection points
      real, dimension(max_num_sta,3) :: pt_prj_k
      real, intent(OUT), dimension(max_num_sta,3) :: pt_prj

      interface
         subroutine FaultRrupPtPrj(fault_top, fualt_base, sta, n_flt, n_sta, r_rup, pt_prj)
           use memory_module
           integer, intent(IN) :: n_flt, n_sta
           real, intent(IN), dimension(max_num_pt_seg,3) :: fault_top, fualt_base
           real, intent(IN), dimension(max_num_sta,3) :: sta
           real, intent(OUT), dimension(max_num_sta) :: r_rup
           real, intent(OUT), dimension(max_num_sta,3) :: pt_prj
         end subroutine

         function FaultRJB(fault_top, fualt_base, sta, n_flt, n_sta)
           use memory_module
           integer, intent(IN) :: n_flt, n_sta
           real, intent(IN), dimension(max_num_pt_seg,3) :: fault_top, fualt_base
           real, intent(IN), dimension(max_num_sta,3) :: sta
           real, dimension(max_num_sta) :: FaultRJB
         end function

         function GC2Coordinates(fault_xyz,n_seg,npt_seg,pt_xyz,n_pt)
           use memory_module
           real, intent(IN), dimension(max_num_pt_seg,3,max_num_seg) :: fault_xyz
           integer, intent(IN) :: n_seg
           integer, intent(IN), dimension(max_num_seg) :: npt_seg
           real, intent(IN), dimension(max_num_sta,3) :: pt_xyz
           integer, intent(IN) :: n_pt
           real, dimension(max_num_sta,2) :: GC2Coordinates
         end function
      end interface 

!     Compute distance metrics
! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
!     compute closest distance and projection pt (R_rup and pr_prj)
      k = 1
      call FaultRrupPtPrj(flt_cor_top(:,:,k),flt_cor_base(:,:,k),sta_cor, &
                          num_pt_seg(k),num_sta,r_rup,pt_prj)
!     compute Joyner-Boore distance (R_jb)
      r_jb  =  FaultRJB(flt_cor_top(:,:,1),flt_cor_base(:,:,1),sta_cor,num_pt_seg(1),num_sta)
      do k = 2,num_seg
!       compute closest distance and projection pt (R_rup and pr_prj)
        call FaultRrupPtPrj(flt_cor_top(:,:,k),flt_cor_base(:,:,k),sta_cor, &
                            & num_pt_seg(k),num_sta,r_rup_k,pt_prj_k)
        do i = 1,num_pt_seg(k)
          if (r_rup_k(i) < r_rup(i) ) then
            r_rup(i) = r_rup_k(i)
            pt_prj(i,:) = pt_prj_k(i,:)
          end if
        end do

!       compute Joyner-Boore distance (R_jb)
        r_jb_k  = FaultRJB(flt_cor_top(:,:,k),flt_cor_base(:,:,k),sta_cor, &
                           & num_pt_seg(k),num_sta)
        r_jb = min(r_jb,r_jb_k)
      end do

!     compute GC2 station coordinates
      tu_pt = GC2Coordinates(flt_cor_top,num_seg,num_pt_seg,sta_cor,num_sta)
      u_pt = tu_pt(:,2)
      t_pt = tu_pt(:,1)


!     compute GC2 fault coordinates
      do k = 1,num_seg
!       copy the coordinates of the k^th segment in flt_seg_k_cor       
        do i = 1,num_pt_seg(k)
          flt_seg_k_cor(i,:) = flt_cor_top(i,:,k)
        end do
!       compute the GC2 coordinates of the k^th segment
        tu_flt_seg_k = GC2Coordinates(flt_cor_top,num_seg,num_pt_seg, &
                                         & flt_cor_top(:,:,k),num_pt_seg(k))
!       copy the GC2 coordinates of the k^th array in tu_fualt
        do i = 1,num_pt_seg(k)
          tu_fault(i,:,k) = tu_flt_seg_k(i,:)
        end do
      end do
      len_fault = maxval(tu_fault(:,2,:))
   
!     compute  distance metrics (Rx, Ry, R_y0)
      r_x = t_pt
      !parallel to strike distance from the center of rupture
      r_y = abs(u_pt - len_fault/2.)
      !parallel to strike distance from the end of rupture
      do k = 1,num_sta
        if ((u_pt(k) >= 0.) .and. (u_pt(k) <= len_fault)) then
!         point is located in the rupture projection
          r_y0(k) = 0.
        else
          r_y0(k) = min(abs(u_pt(k)), abs(u_pt(k) - len_fault)) 
        end if
      end do

      end subroutine dist_metrics_and_prjpt

!--------------------------------------------------------------
!     functions for computing the distance metrics
!--------------------------------------------------------------

!--------------------------------------------------------------
!     Closest distance (R_rup)

      function FaultRrup(fault_top, fault_base, sta, n_flt, n_sta)
!     FaultRrup comptues and returns r_rup at for all stations

      use memory_module
      implicit none
      integer :: i, j
      integer, intent(IN) :: n_flt, n_sta
      real, intent(IN), dimension(max_num_pt_seg,3) :: fault_top  !fault coordinates top and base
      real, intent(IN), dimension(max_num_pt_seg,3) :: fault_base
      real, intent(IN), dimension(max_num_sta,3) :: sta           !station coordinates
      real, dimension(max_num_sta) :: r_rup, FaultRrup            !closest distance
      real :: r_rup_i
      real, dimension(3) :: v1, v2, v3, pt

      interface
         function ClosDistTrig3D(v1,v2,v3,pt)
           real, intent(IN), dimension(3) :: v1, v2, v3, pt
           real :: ClosDistTrig3D
         end function
      end interface

!     iterate over all stations
      do i = 1,n_sta
!        station coordinates
         pt = sta(i,:)

!        vertices of first sub-triangle, 1st segment
         v1 = fault_top(1,:)
         v2 = fault_top(2,:)
         v3 = fault_base(1,:)
!        closest distance to first sub-triangle
         r_rup(i) = ClosDistTrig3D(v1,v2,v3,pt)

!        iterate over all segments (2nd to last)
         do j = 1,(n_flt - 1)

!           vertices of first sub-triangle, segment i
            v1 = fault_top(j,:)
            v2 = fault_top(j+1,:)
            v3 = fault_base(j,:)
!           closest distance to first sub-triangle
            r_rup_i = ClosDistTrig3D(v1,v2,v3,pt)
            r_rup(i) = min(r_rup(i),r_rup_i)

!           vertices of second sub-triangle, segment i
            v1 = fault_top(j+1,:)
            v2 = fault_base(j+1,:)
            v3 = fault_base(j,:)
!           closest distance to first sub-triangle
            r_rup_i = ClosDistTrig3D(v1,v2,v3,pt)
            r_rup(i) = min(r_rup(i),r_rup_i)

         end do
      end do

      FaultRrup = r_rup

      end function FaultRrup

!--------------------------------------------------------------
!     Joyner-Boore Distance (R_JB)

      function FaultRJB(fault_top, fault_base, sta, n_flt, n_sta)
!     FaultRJB comptues and returns r_jb at for all stations

      use memory_module
      implicit none
      integer, intent(IN) :: n_flt, n_sta
      real, intent(IN), dimension(max_num_pt_seg,3) :: fault_top  !fault coordinates top and base
      real, intent(IN), dimension(max_num_pt_seg,3) :: fault_base
      real, intent(IN), dimension(max_num_sta,3) :: sta           !station coordinates
      real, dimension(max_num_pt_seg,3) :: fault_top_elev0, fault_base_elev0
      real, dimension(max_num_sta,3) :: sta_elev0
      real, dimension(max_num_sta) :: FaultRJB                    !Joyner-Boore distance

      interface
         function FaultRrup(fault_top, fault_base, sta, n_flt, n_sta)
           use memory_module
           integer, intent(IN) :: n_flt, n_sta
           real, intent(IN), dimension(max_num_pt_seg,3) :: fault_top, fault_base
           real, intent(IN), dimension(max_num_sta,3) :: sta
           real, dimension(max_num_sta) :: FaultRrup
         end function
      end interface

!     set elevation to zero to compute the horizontal projection dist (r_jb)
      fault_top_elev0 = fault_top
      fault_base_elev0 = fault_base
      sta_elev0 = sta

      fault_top_elev0(:,3) = 0.0
      fault_base_elev0(:,3) = 0.0
      sta_elev0(:,3) = 0.0

!     compute R_jb distance
      FaultRJB = FaultRrup(fault_top_elev0, fault_base_elev0, sta_elev0, n_flt, n_sta)

      end function FaultRJB

!--------------------------------------------------------------
!     Closest distance and projection point (R_rup, pt_prj)

      subroutine FaultRrupPtPrj(fault_top, fault_base, sta, n_flt, n_sta, r_rup, pt_prj)
!     FaultRrup comptues and returns r_rup at for all stations

      use memory_module
      implicit none
      integer :: i, j
      integer, intent(IN) :: n_flt, n_sta
      real, intent(IN), dimension(max_num_pt_seg,3) :: fault_top  !fault coordinates top and base
      real, intent(IN), dimension(max_num_pt_seg,3) :: fault_base
      real, intent(IN), dimension(max_num_sta,3) :: sta           !station coordinates
      real, dimension(3) :: v1, v2, v3, pt
      real :: r_rup_i
      real, dimension(3) :: pt_prj_i
      real, intent(OUT), dimension(max_num_sta) :: r_rup          !closest distance
      real, intent(OUT), dimension(max_num_sta,3) :: pt_prj       !closest distance

      interface
         subroutine ClosDistProjPtTrig3D(v1,v2,v3,pt,clos_dist,pt_prj)
           real, intent(IN), dimension(3) :: v1, v2, v3, pt
           real, intent(OUT) :: clos_dist
           real, intent(OUT), dimension(3) :: pt_prj
         end subroutine
      end interface

!     iterate over all stations
      do i = 1,n_sta
!        station coordinates
         pt = sta(i,:)

!        vertices of first sub-triangle, 1st segment
         v1 = fault_top(1,:)
         v2 = fault_top(2,:)
         v3 = fault_base(1,:)
!        closest distance to first sub-triangle
         call ClosDistProjPtTrig3D(v1,v2,v3,pt,r_rup(i),pt_prj(i,:))

!        iterate over all segments (2nd to last)
         do j = 1,(n_flt - 1)

!           vertices of first sub-triangle, segment i
            v1 = fault_top(j,:)
            v2 = fault_top(j+1,:)
            v3 = fault_base(j,:)
!           closest distance to first sub-triangle
            call ClosDistProjPtTrig3D(v1,v2,v3,pt,r_rup_i,pt_prj_i)
            if (r_rup_i < r_rup(i) ) then
               r_rup(i) = r_rup_i
               pt_prj(i,:) = pt_prj_i
            end if

!           vertices of second sub-triangle, segment i
            v1 = fault_top(j+1,:)
            v2 = fault_base(j+1,:)
            v3 = fault_base(j,:)
!           closest distance to second sub-triangle
            call ClosDistProjPtTrig3D(v1,v2,v3,pt,r_rup_i,pt_prj_i)
            if (r_rup_i < r_rup(i) ) then
               r_rup(i) = r_rup_i
               pt_prj(i,:) = pt_prj_i
            end if

         end do
      end do

      end subroutine FaultRrupPtPrj

