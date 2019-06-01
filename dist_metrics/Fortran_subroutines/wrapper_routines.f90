!------------------------------------------------------------------
!     Wrapper function for communication with R
!
!-------------------------------------------------------------------
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
      integer :: n_pt_flt, n_sta !number of fault coordinates and stations
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
      dm_flt_cor_base(1:n_pt_flt,1:2,1) = flt_xy !fault top coordinates
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

