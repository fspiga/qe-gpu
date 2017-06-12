# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_GPU], [

AC_ARG_ENABLE(gpu,
   [AS_HELP_STRING([--enable-gpu=<arch>],
       [Compile with GPU support (Supported arch: {kepler, pascal}, default: pascal)])],
   [if test "$enableval" = "no" ; then
      use_gpu=0
   else
      use_gpu=1
      gpu_arch="$enableval"
   fi],
   [use_gpu=1 gpu_arch="pascal"])

with_nvtx=0
gpu_arch_num=60

if test "$use_gpu" -eq 1
then
  case "$gpu_arch" in
    kepler | Kepler | KEPLER ) 
      gpu_arch_num=35
      ;;
    pascal | Pascal | PASCAL ) 
      gpu_arch_num=60
      ;;
    *)
      AC_MSG_ERROR([ *** GPU arch must be one of: kepler, pascal ])
      ;;
    esac
else
  AC_MSG_ERROR([ *** GPU acceleration must be enable ])
fi


  try_dflags="$try_dflags -DUSE_CUDA"

  try_iflags="$try_iflags -I\$(TOPDIR)/Eigensolver_gpu-0.1/lib_eigsolve"

  f90flags="$f90flags -Mcuda=cc$gpu_arch_num,cuda8.0 -Mlarge_arrays"

  ldflags="$ldflags -Mcuda=cc$gpu_arch_num,cuda8.0 -Mlarge_arrays"
  ld_libs="$ld_libs -Mcudalib=cufft,cublas,cusolver \$(TOPDIR)/Eigensolver_gpu-0.1/lib_eigsolve/lib_eigsolve.a"

  AC_SUBST(gpu_arch_num)
  AC_CONFIG_FILES([install/Makefile.lib_eigsolve:install/Makefile.lib_eigsolve.in])

  ]
)
