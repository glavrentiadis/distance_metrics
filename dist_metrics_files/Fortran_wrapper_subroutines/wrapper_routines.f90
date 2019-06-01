!  *****************************************************************
!  *  Wrapper functions for communication with R                   *
!  *                                                               *
!  *****************************************************************
!     This function uses the dist_metrics subroutine to compute the fualt
!     distance metrics and retruns in R only the GC2 U and T coordinates.
!     Currently it is developed for a single segment fault 
!     For R to be able to communicate with FORTRAN, input and output arrays
!     must be double precision

! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
!     Subroutine for exporting the station GC2 coordinates to R

      subroutine r_wrapper_gc2(flt_xy,n_pt_flt,sta_xy,n_sta, sta_ut)

!     declare variables
      use memory_module
      implicit none
!---- input arguments
      integer, intent(IN) :: n_pt_flt, n_sta !number of fault coordinates and stations
      double precision, intent(IN), dimension(n_pt_flt,2) :: flt_xy !coordinates of fault segments
      double precision, intent(IN), dimension(n_sta,2) :: sta_xy !coordinates of stations
!---- output arguments
      double precision, intent(OUT), dimension(n_sta,2) :: sta_ut
!---- internal variables
      integer :: dm_num_seg, dm_num_sta
      integer, dimension(max_num_seg) :: dm_num_pt_seg
      real, dimension(max_num_pt_seg,3,max_num_seg) :: dm_flt_cor_top
      real, dimension(max_num_pt_seg,3,max_num_seg) :: dm_flt_cor_base 
      real, dimension(max_num_sta,3) :: dm_sta_cor
      real, dimension(max_num_sta) :: dm_r_rup, dm_r_jb
      real, dimension(max_num_sta) :: dm_r_x, dm_r_y, dm_r_y0
      real, dimension(max_num_sta) :: dm_u_pt, dm_t_pt

      interface
        subroutine dist_metrics(num_seg,num_pt_seg,flt_cor_top,flt_cor_base,num_sta,sta_cor, &
                               & r_rup,r_jb,r_x,r_y,r_y0,u_pt,t_pt)
          use memory_module
          integer, intent(IN) :: num_seg, num_sta
          integer, intent(IN), dimension(max_num_seg) :: num_pt_seg
          real, intent(IN), dimension(max_num_pt_seg,3,max_num_seg) :: flt_cor_top
          real, intent(IN), dimension(max_num_pt_seg,3,max_num_seg) :: flt_cor_base 
          real, intent(IN), dimension(max_num_sta,3) :: sta_cor
          real, intent(OUT), dimension(max_num_sta) :: r_rup, r_jb
          real, intent(OUT), dimension(max_num_sta) :: r_x, r_y, r_y0
          real, intent(OUT), dimension(max_num_sta) :: u_pt, t_pt
        end subroutine
      end interface 

!     Pass info to arrays for the dist_metrics subroutine
      dm_num_seg = 1 !number of segments
      dm_num_pt_seg(1) = n_pt_flt !number of fualt coordinates
      dm_flt_cor_top(:,:,:) = 0.0 !initialize fault top coordinates
      dm_flt_cor_top(1:n_pt_flt,1:2,1) = flt_xy !fault top coordinates
      dm_flt_cor_top(1:n_pt_flt,3,1) = 0.0
      dm_flt_cor_base(:,:,:) = 0.0 !initialize fault base coordinates
      dm_flt_cor_base(1:n_pt_flt,1:2,1) = flt_xy !fault base coordinates
      dm_flt_cor_base(1:n_pt_flt,3,1) = -15.0 !default depth
      dm_num_sta = n_sta !number of stations
      dm_sta_cor(:,:) = 0.0 !initialize sation coordinates
      dm_sta_cor(1:n_sta,1:2) = sta_xy !sation coordinates
      dm_sta_cor(1:n_sta,3) = 0.0

!     call dist_metrics subroutine to compute distance metrics
      call dist_metrics(dm_num_seg,dm_num_pt_seg,dm_flt_cor_top,dm_flt_cor_base, &
                       & dm_num_sta,dm_sta_cor, &
                       & dm_r_rup,dm_r_jb,dm_r_x,dm_r_y,dm_r_y0,dm_u_pt,dm_t_pt)     

!     return U and T coordinates
      sta_ut(:,1) = dm_u_pt(1:n_sta)
      sta_ut(:,2) = dm_t_pt(1:n_sta)

      sta_ut(:,1) = dm_u_pt(1:n_sta)
      sta_ut(:,2) = dm_t_pt(1:n_sta) 

      end subroutine r_wrapper_gc2


!  *****************************************************************
!  *   Wrapper functions for communication with Python             *
!  *                                                               *
!  *****************************************************************
!     This functions uses the dist_metrics subroutine to compute the
!     fault distance metrics and returns them to Python
!     Currently it is developed for a single segment fault 
!
! -----------------------------------------------------------
! -----------------------------------------------------------
!     Only distance metrics
! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
!     Subroutine can be used for exporting the distance metrics to Python for multiple stations

      subroutine py_wrapper_dist2stas(flt_xyz_top,flt_xyz_base,n_pt_flt,sta_xyz,n_sta, dist2sta)

!     declare variables
      use memory_module
      implicit none
!---- input arguments
      integer, intent(IN) :: n_pt_flt, n_sta !number of fault coordinates and stations
      double precision, intent(IN), dimension(n_pt_flt,3) :: flt_xyz_top !coordinates of top of fault segments
      double precision, intent(IN), dimension(n_pt_flt,3) :: flt_xyz_base !coordinates of base of fault segments
      double precision, intent(IN), dimension(n_sta,3) :: sta_xyz !coordinates of stations
!---- output arguments
      !different distance metrics between fault and station
      ![R_rup,R_jb,R_x,R_y,R_y0,u_pt,t_pt]
      double precision, intent(OUT), dimension(n_sta,7) :: dist2sta
!---- internal variables
      integer :: dm_num_seg, dm_num_sta
      integer, dimension(max_num_seg) :: dm_num_pt_seg
      real, dimension(max_num_pt_seg,3,max_num_seg) :: dm_flt_cor_top
      real, dimension(max_num_pt_seg,3,max_num_seg) :: dm_flt_cor_base 
      real, dimension(max_num_sta,3) :: dm_sta_cor
      real, dimension(max_num_sta) :: dm_r_rup, dm_r_jb
      real, dimension(max_num_sta) :: dm_r_x, dm_r_y, dm_r_y0
      real, dimension(max_num_sta) :: dm_u_pt, dm_t_pt

      interface
        subroutine dist_metrics(num_seg,num_pt_seg,flt_cor_top,flt_cor_base,num_sta,sta_cor, &
                               & r_rup,r_jb,r_x,r_y,r_y0,u_pt,t_pt)
          use memory_module
          integer, intent(IN) :: num_seg, num_sta
          integer, intent(IN), dimension(max_num_seg) :: num_pt_seg
          real, intent(IN), dimension(max_num_pt_seg,3,max_num_seg) :: flt_cor_top
          real, intent(IN), dimension(max_num_pt_seg,3,max_num_seg) :: flt_cor_base 
          real, intent(IN), dimension(max_num_sta,3) :: sta_cor
          real, intent(OUT), dimension(max_num_sta) :: r_rup, r_jb
          real, intent(OUT), dimension(max_num_sta) :: r_x, r_y, r_y0
          real, intent(OUT), dimension(max_num_sta) :: u_pt, t_pt
        end subroutine
      end interface 

!     Pass info to arrays for the dist_metrics subroutine
      dm_num_seg = 1 !number of segments
      dm_num_pt_seg(1) = n_pt_flt !number of fualt coordinates
      dm_flt_cor_top(:,:,:) = 0.0 !initialize fault top coordinates
      dm_flt_cor_top(1:n_pt_flt,1:3,1) = flt_xyz_top !fault top coordinates
      dm_flt_cor_base(:,:,:) = 0.0 !initialize fault base coordinates
      dm_flt_cor_base(1:n_pt_flt,1:3,1) = flt_xyz_base !fault base coordinates
      dm_num_sta = n_sta !number of stations
      dm_sta_cor(:,:) = 0.0 !initialize sation coordinates
      dm_sta_cor(1:n_sta,1:3) = sta_xyz !sation coordinates

!     call dist_metrics subroutine to compute distance metrics
      call dist_metrics(dm_num_seg,dm_num_pt_seg,dm_flt_cor_top,dm_flt_cor_base, &
                       & dm_num_sta,dm_sta_cor, &
                       & dm_r_rup,dm_r_jb,dm_r_x,dm_r_y,dm_r_y0,dm_u_pt,dm_t_pt)     

!     distance metrics
      dist2sta(:,1) = dm_r_rup(1:n_sta)
      dist2sta(:,2) = dm_r_jb(1:n_sta)
      dist2sta(:,3) = dm_r_x(1:n_sta)
      dist2sta(:,4) = dm_r_y(1:n_sta)
      dist2sta(:,5) = dm_r_y0(1:n_sta)
      dist2sta(:,6) = dm_u_pt(1:n_sta)
      dist2sta(:,7) = dm_t_pt(1:n_sta)

      end subroutine py_wrapper_dist2stas

! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
!     Subroutine can be used for exporting the distance metrics to Python for a single station

      subroutine py_wrapper_dist2sta(flt_xyz_top,flt_xyz_base,n_pt_flt,sta_xyz, dist2sta)

!     declare variables
      use memory_module
      implicit none
!---- input arguments
      integer, intent(IN) :: n_pt_flt !number of fault coordinates
      double precision, intent(IN), dimension(n_pt_flt,3) :: flt_xyz_top !coordinates of top of fault segments
      double precision, intent(IN), dimension(n_pt_flt,3) :: flt_xyz_base !coordinates of base of fault segments
      double precision, intent(IN), dimension(3) :: sta_xyz !coordinates of stations
!---- output arguments
      !different distance metrics between fault and station
      ![R_rup,R_jb,R_x,R_y,R_y0,u_pt,t_pt]
      double precision, intent(OUT), dimension(7) :: dist2sta 
!---- internal variables
      integer :: n_sta
      double precision, dimension(1,7) :: dist2stas !distance metrics, interface between subroutines

      interface
        subroutine Py_wrapper_dist2stas(flt_xyz_top,flt_xyz_base,n_pt_flt,sta_xyz,n_sta,dist2stas)
          integer, intent(IN) :: n_pt_flt, n_sta
          double precision, intent(IN), dimension(n_pt_flt,3) :: flt_xyz_top
          double precision, intent(IN), dimension(n_pt_flt,3) :: flt_xyz_base
          double precision, intent(IN), dimension(n_sta,3) :: sta_xyz
          double precision, intent(OUT), dimension(n_sta,7) :: dist2stas
        end subroutine
      end interface 

      n_sta = 1 !distances are calculated for a single station

      call Py_wrapper_dist2stas(flt_xyz_top,flt_xyz_base,n_pt_flt,sta_xyz,n_sta, dist2stas)
      dist2sta = dist2stas(1,:)

      end subroutine py_wrapper_dist2sta

! -----------------------------------------------------------
! -----------------------------------------------------------
!     Distance metrics and projection points
! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
!     Subroutine can be used for exporting the distance metrics and projection points 
!     to Python for multiple stations

      subroutine py_wrapper_dist_prj2stas(flt_xyz_top,flt_xyz_base,n_pt_flt,sta_xyz,n_sta, dist_prj2sta)

!     declare variables
      use memory_module
      implicit none
!---- input arguments
      integer, intent(IN) :: n_pt_flt, n_sta !number of fault coordinates and stations
      double precision, intent(IN), dimension(n_pt_flt,3) :: flt_xyz_top !coordinates of top of fault segments
      double precision, intent(IN), dimension(n_pt_flt,3) :: flt_xyz_base !coordinates of base of fault segments
      double precision, intent(IN), dimension(n_sta,3) :: sta_xyz !coordinates of stations
!---- output arguments
      !different distance metrics between fault and station
      ![R_rup,R_jb,R_x,R_y,R_y0,u_pt,t_pt]
      double precision, intent(OUT), dimension(n_sta,10) :: dist_prj2sta
!---- internal variables
      integer :: dm_num_seg, dm_num_sta
      integer, dimension(max_num_seg) :: dm_num_pt_seg
      real, dimension(max_num_pt_seg,3,max_num_seg) :: dm_flt_cor_top
      real, dimension(max_num_pt_seg,3,max_num_seg) :: dm_flt_cor_base 
      real, dimension(max_num_sta,3) :: dm_sta_cor
      real, dimension(max_num_sta) :: dm_r_rup, dm_r_jb
      real, dimension(max_num_sta) :: dm_r_x, dm_r_y, dm_r_y0
      real, dimension(max_num_sta) :: dm_u_pt, dm_t_pt
      real, dimension(max_num_sta,3) :: dm_pt_prj

      interface
        subroutine dist_metrics_and_prjpt(num_seg,num_pt_seg,flt_cor_top,flt_cor_base, &
                                         & num_sta,sta_cor, &
                                         & r_rup,r_jb,r_x,r_y,r_y0,u_pt,t_pt,pt_prj)
          use memory_module
          integer, intent(IN) :: num_seg, num_sta
          integer, intent(IN), dimension(max_num_seg) :: num_pt_seg
          real, intent(IN), dimension(max_num_pt_seg,3,max_num_seg) :: flt_cor_top
          real, intent(IN), dimension(max_num_pt_seg,3,max_num_seg) :: flt_cor_base 
          real, intent(IN), dimension(max_num_sta,3) :: sta_cor
          real, intent(OUT), dimension(max_num_sta) :: r_rup, r_jb
          real, intent(OUT), dimension(max_num_sta) :: r_x, r_y, r_y0
          real, intent(OUT), dimension(max_num_sta) :: u_pt, t_pt
          real, intent(OUT), dimension(max_num_sta,3) :: pt_prj
        end subroutine
      end interface 

!     Pass info to arrays for the dist_metrics subroutine
      dm_num_seg = 1 !number of segments
      dm_num_pt_seg(1) = n_pt_flt !number of fualt coordinates
      dm_flt_cor_top(:,:,:) = 0.0 !initialize fault top coordinates
      dm_flt_cor_top(1:n_pt_flt,1:3,1) = flt_xyz_top !fault top coordinates
      dm_flt_cor_base(:,:,:) = 0.0 !initialize fault base coordinates
      dm_flt_cor_base(1:n_pt_flt,1:3,1) = flt_xyz_base !fault base coordinates
      dm_num_sta = n_sta !number of stations
      dm_sta_cor(:,:) = 0.0 !initialize sation coordinates
      dm_sta_cor(1:n_sta,1:3) = sta_xyz !sation coordinates

!     call dist_metrics subroutine to compute distance metrics
      call dist_metrics_and_prjpt(dm_num_seg,dm_num_pt_seg,dm_flt_cor_top,dm_flt_cor_base, &
                                 & dm_num_sta,dm_sta_cor, &
                                 & dm_r_rup,dm_r_jb,dm_r_x,dm_r_y,dm_r_y0,dm_u_pt,dm_t_pt, &
                                 & dm_pt_prj)     

!     distance metrics
      dist_prj2sta(:,1) = dm_r_rup(1:n_sta)
      dist_prj2sta(:,2) = dm_r_jb(1:n_sta)
      dist_prj2sta(:,3) = dm_r_x(1:n_sta)
      dist_prj2sta(:,4) = dm_r_y(1:n_sta)
      dist_prj2sta(:,5) = dm_r_y0(1:n_sta)
      dist_prj2sta(:,6) = dm_u_pt(1:n_sta)
      dist_prj2sta(:,7) = dm_t_pt(1:n_sta)
      dist_prj2sta(:,8:10) = dm_pt_prj(:,:)

      end subroutine py_wrapper_dist_prj2stas

! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
!     Subroutine can be used for exporting the distance metrics to Python for a single station

      subroutine py_wrapper_dist_prj2sta(flt_xyz_top,flt_xyz_base,n_pt_flt,sta_xyz, dist_prj2sta)

!     declare variables
      use memory_module
      implicit none
!---- input arguments
      integer, intent(IN) :: n_pt_flt !number of fault coordinates
      double precision, intent(IN), dimension(n_pt_flt,3) :: flt_xyz_top !coordinates of top of fault segments
      double precision, intent(IN), dimension(n_pt_flt,3) :: flt_xyz_base !coordinates of base of fault segments
      double precision, intent(IN), dimension(3) :: sta_xyz !coordinates of stations
!---- output arguments
      !different distance metrics between fault and station
      ![R_rup,R_jb,R_x,R_y,R_y0,u_pt,t_pt,pt_prj_x,pt_prj_y,pt_prj_z]
      double precision, intent(OUT), dimension(10) :: dist_prj2sta 
!---- internal variables
      integer :: n_sta
      double precision, dimension(1,10) :: dist_prj2stas !distance metrics, interface between subroutines

      interface
        subroutine py_wrapper_dist_prj2stas(flt_xyz_top,flt_xyz_base,n_pt_flt,sta_xyz,n_sta, & 
                                           & dist_prj2stas)
          integer, intent(IN) :: n_pt_flt, n_sta
          double precision, intent(IN), dimension(n_pt_flt,3) :: flt_xyz_top
          double precision, intent(IN), dimension(n_pt_flt,3) :: flt_xyz_base
          double precision, intent(IN), dimension(n_sta,3) :: sta_xyz
          double precision, intent(OUT), dimension(n_sta,10) :: dist_prj2stas
        end subroutine
      end interface

      n_sta = 1 !distances are calculated for a single station

      call py_wrapper_dist_prj2stas(flt_xyz_top,flt_xyz_base,n_pt_flt,sta_xyz,n_sta, dist_prj2stas)
      dist_prj2sta = dist_prj2stas(1,:)

      end subroutine py_wrapper_dist_prj2sta


