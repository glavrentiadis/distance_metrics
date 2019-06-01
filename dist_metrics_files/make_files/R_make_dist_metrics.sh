# Creates the dll FORTRAN library for the distance metrics
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

gfortran -fPIC -c distance_metrics_subroutine.f90 wrapper_routines.f90 closest_dist_trig.f90 closest_dist_proj_pt_trig.f90 GC2_coordinates.f90 gen_functions.f90

R CMD SHLIB -o R_dist_metrics.so distance_metrics_subroutine.f90 wrapper_routines.f90 closest_dist_trig.f90 closest_dist_proj_pt_trig.f90 GC2_coordinates.f90 gen_functions.f90

#move complide file and delete folder
mv R_dist_metrics.so ../../../dist_metrics_progs/.
cd ..
rm -r build

