!---------------------------------------------------------------------
!     General functions
!     This file contains functions for general purpose calculations
! 
!     Author: G. Lavrentiadis
!     Date:   04/12/2019
!     Rev:    0.0
! 
!     Contents:	       Matrix inversion functions
! 
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!    MATRIX INVERSION FUNCTIONS
!---------------------------------------------------------------------

      function matinv2(A)
!     Performs a direct calculation of the inverse of a 2×2 matrix.
      real, intent(in), dimension(2,2) :: A   !! Matrix
      real, dimension(2,2)             :: B   !! Inverse matrix
      real                             :: detinv
      real, dimension(2,2)             :: matinv2

!     Calculate the inverse determinant of the matrix
      detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))

!     Calculate the inverse of the matrix
      B(1,1) = +detinv * A(2,2)
      B(2,1) = -detinv * A(2,1)
      B(1,2) = -detinv * A(1,2)
      B(2,2) = +detinv * A(1,1)

      matinv2 = B

      end function

      function matinv3(A)
!     Performs a direct calculation of the inverse of a 3×3 matrix.
      real, intent(in), dimension(3,3) :: A   !! Matrix
      real, dimension(3,3)             :: B   !! Inverse matrix
      real                             :: detinv
      real, dimension(3,3)             :: matinv3

!     Calculate the inverse determinant of the matrix
      detinv = 1/(A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                 - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                 + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))

!     Calculate the inverse of the matrix
      B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
      B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
      B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
      B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
      B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
      B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
      B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
      B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
      B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))

      matinv3 = B

      end function
