# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_SCALAPACK], [

have_scalapack=0

AC_ARG_WITH(scalapack,
   [AS_HELP_STRING([--with-scalapack],
       [(yes|no) Use scalapack if available (default: no)])],
   [if  test "$withval" = "yes" ; then
      with_scalapack=1
   elif  test "$withval" = "no" ; then
      with_scalapack=0
   fi],
   [with_scalapack=0])

  # final check on availability of parallel environment
  for dummy in x # to allow simple 'break'
  do
    test "$have_mpi" -eq 0 && break

    F77=$mpif90
    LIBS="$mpi_libs"

    # look for scalapack if required
    test "$with_scalapack" -eq 0 && break

    if test "$have_mkl" -eq 1
    then
      unset ac_cv_search_pdgemr2d # clear cached value
      LIBS="-lmkl_blacs_lp64 $mpi_libs $math_libs"
      scalapack_libs=-lmkl_blacs_openmpi_lp64
      AC_SEARCH_LIBS(pdgemr2d, "mkl_scalapack_lp64" , have_scalapack=1
                     try_dflags="$try_dflags -D__SCALAPACK"
                     scalapack_libs="-lmkl_scalapack_lp64 $scalapack_libs" )
      test "$have_scalapack" -eq 1 && break
    fi 
  done
   
  # Configuring output message
  if test "$have_scalapack" -eq 1; then
    scalapack_line="SCALAPACK_LIBS=$scalapack_libs"
  else
    scalapack_libs=""
    scalapack_line="@delete@"
  fi
  
  AC_SUBST(scalapack_libs)
  AC_SUBST(scalapack_line)
  ]
)
