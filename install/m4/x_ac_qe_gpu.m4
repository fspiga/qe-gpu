# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_GPU], [

AC_ARG_ENABLE(gpu,
   [AS_HELP_STRING([--enable-gpu=<arch|no>],
       [Compile with GPU support (Supported arch: {Kepler, Pascal}, default: no)])],
   [if   test "$enableval" = "yes" ; then
      use_gpu=1
      gpu_arch="$withval"
   else
      use_gpu=0
   fi],
   [use_gpu=0])


AC_ARG_WITH(gpu-profiling,
   [AS_HELP_STRING([--with-gpu-profiling],
       [Enable NVTX profiling information (default: no)])],
   [if  test "$withval" = "no" ; then
      with_nvtx=0
   else
      with_nvtx=1
   fi],
   [with_nvtx=0])

AC_ARG_WITH(cuda_dir,
   [AS_HELP_STRING([--with-cuda-dir=<path>],
    [Specify CUDA installation directory (default is /usr/local/cuda/ - MANDATORY)])],
    [cuda_path="$withval"],
    [cuda_path="/usr/local/cuda/"])

with_nvtx=0
use_gpu=0

  ]
)
