#source files
SRC_DIR=source/calc_subroutines/
SRC=${SRC_DIR}distance_metrics_subroutine.f90 ${SRC_DIR}closest_dist_trig.f90 \
    ${SRC_DIR}closest_dist_proj_pt_trig.f90   ${SRC_DIR}GC2_coordinates.f90 ${SRC_DIR}gen_functions.f90
#main functions
SRC_MDIR=source/main_program/
#wrapper functions source files
SRC_WDIR=source/wrapper_subroutines/

#program output directory
PRGDIR=./progs/

#compile stand alone programs
FC=gfortran
FFLAGS=-ffixed-line-length-none
#compilers for python and R API
#R
FC_R=R CMD SHLIB
#Python
FC_PY=f2py3
FFLAGS_PY=--f90flags=-ffixed-line-length-none

all: clean build build_Py build_R

build: clean
	mkdir -p $(PRGDIR)
	#compile main program (w/o prj points)
	$(FC) -o $(PRGDIR)dist_metrics.bin     $(FFLAGS) $(SRC) $(SRC_MDIR)distance_metrics_main_program.f90 
	#compile main program (w/ projection points)
	$(FC) -o $(PRGDIR)dist_metrics_prj.bin $(FFLAGS) $(SRC) $(SRC_MDIR)distance_metrics_prj_main_program.f90 
	make clean

build_Py: clean
	mkdir -p $(PRGDIR)
	$(FC_PY) -c -m py_dist_metrics $(FFLAGS_PY) $(SRC) $(SRC_WDIR)wrapper_routines.f90 \
		    	only: py_wrapper_dist2stas py_wrapper_dist2sta py_wrapper_dist_prj2stas \
		        py_wrapper_dist_prj2sta :
	mv py_dist_metrics* $(PRGDIR)py_dist_metrics.so
	make clean
	
build_R: clean
	mkdir -p $(PRGDIR)
	$(FC_R) -o $(PRGDIR)R_dist_metrics.so $(SRC) $(SRC_WDIR)wrapper_routines.f90 
	make clean

clean:
	rm -f *.o *.mod *.mod0
	rm -f $(SRC_DIR)*.o
	rm -f $(SRC_WDIR)*.o
