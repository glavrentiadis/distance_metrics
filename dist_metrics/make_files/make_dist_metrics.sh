# Creates the dll FORTRAN library for the distance metrics
# usage: ./make_dist_metrics.sh
# Date: March 25, 2019

#remove all compiled files
rm *.o *.so
rm *.mod*

#compile Fortran files
gfortran -fPIC -c distance_metrics_subroutine.f90 closest_dist_trig.f90 GC2_coordinates.f90 R_wrapper_gc2_cor.f90

#compile Fortran files for R
R CMD SHLIB distance_metrics_subroutine.f90 closest_dist_trig.f90 GC2_coordinates.f90 R_wrapper_gc2_cor.f90
#rename dll 
mv distance_metrics_subroutine.so dist_metrics.so

#remove fortran objects
rm *.o
