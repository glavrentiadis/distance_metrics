# Creates the dll FORTRAN library for the distance metrics and the wrapper functions for Python
# usage: ./make_dist_metrics.sh
# Date: March 25, 2019

#create build folder
mkdir build

#copy necessary FORTRAN subroutines into build folder
cd ../Fortran_calc_subroutines/
cp * ./../make_files/build/.
cd ../Fortran_wrapper_subroutines/
cp wrapper_routines.f90 ./../make_files/build/.

#compfile FORTRAN files
cd ../make_files/build

#gfortran -fPIC -c distance_metrics_subroutine.f90 closest_dist_trig.f90 closest_dist_proj_pt_trig.f90 GC2_coordinates.f90 gen_functions.f90

f2py3 -c --f90flags=-ffixed-line-length-none wrapper_routines.f90 distance_metrics_subroutine.f90 closest_dist_trig.f90 closest_dist_proj_pt_trig.f90 GC2_coordinates.f90 gen_functions.f90 -m dist_metrics_prj

#remove build folder
cd ../
rm -r build

