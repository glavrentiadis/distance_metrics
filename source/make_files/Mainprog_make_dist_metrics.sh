# Creates the dll FORTRAN library for the distance metrics
# usage: ./make_dist_metrics.sh
# Date: March 25, 2019

#create build folder
mkdir build

#copy necessary FORTRAN subroutines into build folder
cd ../Fortran_calc_subroutines/
cp * ./../make_files/build/.
cd ../Fortran_main_program/
cp * ./../make_files/build/.

#compfile FORTRAN files
cd ../make_files/build

#create objects
gfortran -c -ffixed-line-length-none distance_metrics_subroutine.f90 closest_dist_trig.f90 closest_dist_proj_pt_trig.f90 GC2_coordinates.f90 gen_functions.f90 distance_metrics_main_program.f90 distance_metrics_prj_main_program.f90

#compile main program (w/o prj points)
gfortran -o dist_metrics.bin -ffixed-line-length-none distance_metrics_main_program.f90 distance_metrics_subroutine.f90 closest_dist_trig.f90 closest_dist_proj_pt_trig.f90 GC2_coordinates.f90 gen_functions.f90

#compile main program (w/ projection points)
gfortran -o dist_metrics_prj.bin -ffixed-line-length-none distance_metrics_prj_main_program.f90 distance_metrics_subroutine.f90 closest_dist_trig.f90 closest_dist_proj_pt_trig.f90 GC2_coordinates.f90 gen_functions.f90


#move complide file and delete folder
mv *.bin ../../../dist_metrics_progs/.
cd ..
rm -r build

