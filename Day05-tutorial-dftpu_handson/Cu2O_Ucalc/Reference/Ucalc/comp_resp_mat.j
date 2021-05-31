#!/bin/bash -l

ifort -o r.x resp_mat.f90 -L/usr/global/intel/mkl/10.3.1.107/mkl/lib/intel64 -lmkl_intel_lp64  -lmkl_sequential -lmkl_core
