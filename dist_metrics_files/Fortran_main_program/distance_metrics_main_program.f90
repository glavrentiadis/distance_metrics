!------------------------------------------------------------------
!     Distance Metrics Main
!     It is a stand alone program that computes common distance metrics used in GMPES
! 
!-------------------------------------------------------------------

! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
!     Main program for calculating distance metrics

      PROGRAM dist_metrics_main_prog

!     declare variables
      use memory_module
      implicit none
      integer :: i, j, k  !loop indices
      integer :: num_seg_top, num_seg_base, num_sta  !num of segmetns at top, base, num of stations
      integer, dimension(max_num_seg) :: num_pt_seg_top  !num of points per segment at top
      integer, dimension(max_num_seg) :: num_pt_seg_base !num of points per segment at base
      character (len = 200) :: fname_fault_top, fname_fault_base  !file names fault
      character (len = 200) :: fname_sta !file names stations
      !coor. at the top of the fault, 3rd dim for diff. seg.
      real, dimension(max_num_pt_seg,3,max_num_seg) :: flt_cor_top
      !coor. at the base of the fault, 3rd dim for diff. seg.
      real, dimension(max_num_pt_seg,3,max_num_seg) :: flt_cor_base 
      real, dimension(max_num_sta,3) :: sta_cor  !stations' coordinates
!-----distance metrics arrays
      real, dimension(max_num_sta) :: r_rup, r_jb
      real, dimension(max_num_sta) :: r_x, r_y, r_y0
      real, dimension(max_num_sta) :: u_pt, t_pt
!-----error handling
      integer :: f_iostat
      character (256) :: f_iomsg

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

!     I/O section, read fault and station coordinates
! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 

!     read file name of fault coordinates
      read(*,'(a)') fname_fault_top
      read(*,'(a)') fname_fault_base
!     read file name for station coordinates
      read(*,'(a)') fname_sta

!     read fault top coordinates
!     open fault file
      open(unit = 10, file = fname_fault_top, action = 'READ', iostat=f_iostat, iomsg=f_iomsg)
      if(f_iostat /= 0) then
        write(*,*) 'Open file10 failed with iostat = ', f_iostat, ' iomsg = '//trim(f_iomsg) 
      end if

      read(10,*) num_seg_top                        !read number of segments
      do k = 1,num_seg_top
        read(10,*) num_pt_seg_top(k)                !read number of points on segment k
        do i = 1,num_pt_seg_top(k) 
          read(10,*) (flt_cor_top(i,j,k) , j = 1,3) !read fault coordinates line by line (x,y,z)
        end do
      end do
      close(10)                                     !close fualt file

!     read fault base coordinates
!     open fault file
      open(unit = 11, file = fname_fault_base, action = 'READ', iostat=f_iostat, iomsg=f_iomsg)
      if(f_iostat /= 0) then
        write(*,*) 'Open file11 failed with iostat = ', f_iostat, ' iomsg = '//trim(f_iomsg) 
      end if

      read(11,*) num_seg_base                       !read number of segments
      do k = 1,num_seg_base
        read(11,*) num_pt_seg_base(k)
        do i = 1,num_pt_seg_base(k)                 !read number of points on segment k
          read(11,*) (flt_cor_base(i,j,k) , j = 1,3)!read fault coordinates line by line (x,y,z)
        end do
      end do
      close(11)                                     !close fualt file

!     read station coordiantes
!     open fault file
      open(unit = 12, file = fname_sta, action = 'READ', iostat=f_iostat, iomsg=f_iomsg)
      if(f_iostat /= 0) then
        write(*,*) 'Open file12 failed with iostat = ', f_iostat, ' iomsg = '//trim(f_iomsg) 
      end if

      read(12,*) num_sta                            !number of stations to read
      do i = 1,num_sta
        read(12,*) (sta_cor(i,j) , j = 1,3)         !read staion coordinates line by line (x,y,z)
      end do
      close(12)                                     !close station file

!     Input Assertions
! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
!     check number of segments at the top and base
      if (num_seg_top /= num_seg_base) then
        stop 'Error. incompatible number of top and base fault segments'
      end if
!     check number of points per segment
      do k = 1,num_seg_top
        if (num_pt_seg_top(k) /= num_pt_seg_base(k)) then
          stop 
        end if
      end do

!     Compute distance metrics
! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
      call dist_metrics(num_seg_top,num_pt_seg_top,flt_cor_top,flt_cor_base, num_sta,sta_cor, &
                       & r_rup,r_jb,r_x,r_y,r_y0,u_pt,t_pt)

!     Write distance metrics
! --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- 
!     write fault coordinates, top
      open(unit = 21, file = 'distance_metrics.txt', action = 'WRITE', iostat=f_iostat, iomsg=f_iomsg)
      if(f_iostat /= 0) then
        write(*,*) 'Open file21 failed with iostat = ', f_iostat, ' iomsg = '//trim(f_iomsg) 
      end if

!     write file header
      write(21,*) 'Sta ID ','x ','y ','z ','R_rup ','R_jb ','R_x ','R_y ','R_y0 ','u ','t '
!     write distance metrics
      do i = 1,num_sta
         write(21,*) i, (sta_cor(i,j), j = 1,3), r_rup(i), r_jb(i), r_x(i), r_y(i), r_y0(i), &
                                                 u_pt(i), t_pt(i)
      end do


      END PROGRAM dist_metrics_main_prog




