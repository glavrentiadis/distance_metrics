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
cp py_dist_metrics_prj.fpy ./../make_files/build/.

#compfile FORTRAN files
cd ../make_files/build

#gfortran -fPIC -c distance_metrics_subroutine.f90 closest_dist_trig.f90 closest_dist_proj_pt_trig.f90 GC2_coordinates.f90 gen_functions.f90

## f2py commands for signature files
#f2py3 -h py_dist_metrics_prj.pyf distance_metrics_subroutine.f90 wrapper_routines.f90 closest_dist_trig.f90 closest_dist_proj_pt_trig.f90 GC2_coordinates.f90 gen_functions.f90 only: py_wrapper_dist2stas py_wrapper_dist2stas py_wrapper_dist_prj2stas py_wrapper_dist_prj2sta:

## f2py compile commands
#f2py3 -c --f90flags=-ffixed-line-length-none py_dist_metrics_prj.pyf distance_metrics_subroutine.f90 wrapper_routines.f90 closest_dist_trig.f90 closest_dist_proj_pt_trig.f90 GC2_coordinates.f90 gen_functions.f90 -m py_dist_metrics_prj only: py_wrapper_dist2stas py_wrapper_dist2stas py_wrapper_dist_prj2stas py_wrapper_dist_prj2sta:

## f2py, compile w/o signature
f2py3 -c -m py_dist_metrics --f90flags=-ffixed-line-length-none distance_metrics_subroutine.f90 wrapper_routines.f90 closest_dist_trig.f90 closest_dist_proj_pt_trig.f90 GC2_coordinates.f90 gen_functions.f90 only: py_wrapper_dist2stas py_wrapper_dist2sta py_wrapper_dist_prj2stas py_wrapper_dist_prj2sta :


#remove build folder
mv py_dist_metrics.*.so ../../../dist_metrics_progs/py_dist_metrics.so
cd ../
rm -r build
#cp py_dist_metrics.so ../../testing/.


