#source /opt/intel/bin/compilervars.sh ia32
source /opt/intel/bin/compilervars.sh intel64
#ifort -o r.x resp_mat.f90 -L/opt/intel/mkl/8.1.1/lib/32/ -lmkl_lapack \
#                                 -lmkl_ia32 -lguide -lpthread
#ifort -o r.x resp_mat.f90 -L/opt/intel/mkl/lib/ia32/ -lmkl_lapack95 -lmkl -lpthread #\
#                                 -lmkl_blacs -lmkl_intel # -lguide -lpthread
#ifort -o r.x resp_mat.f90 -L/opt/intel/mkl/lib/ia32/libmkl_intel.a -L/opt/intel/mkl/lib/ia32/libmkl_lapack95.a -L/opt/intel/mkl/lib/ia32/libmkl_blas95.a # -l/opt/intel/mkl/lib/ia32/libmkl_intel.a #-lmkl_lapack95  -lmkl_blas95.a
ifort -o r.x resp_mat.f90 -L/opt/intel/mkl/lib/intel64/ -I/opt/intel/mkl/include -lmkl_intel_lp64 -lmkl_sequential -lmkl_core #-lguide -lpthread # -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
