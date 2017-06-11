# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation
# Copyright (C) 2017      Filippo Spiga

AC_DEFUN([X_AC_QE_MATHLIB], [

have_mathlib=0
have_mkl=0

ld_library_path=`echo $LD_LIBRARY_PATH | sed 's/:/ /g'`
try_libdirs="none $ld_library_path $libdirs $try_libdirs"

AC_ARG_WITH(netlib,
   [AS_HELP_STRING([--with-netlib],
       [compile with Netlib LAPACK and BLAS (default: no)])],
    [if test "$withval" = "yes" ; then
      use_netlib=1
   else
      use_netlib=0
   fi],
   [use_netlib=0])

if test "$use_netlib" -eq 0
then
   if test "$math_libs" = ""
   then
        case "$arch:$f90" in
        x86_64:pgf* )

                # Check first MKL...
                for dir in $try_libdirs
                do
                        unset ac_cv_search_dgemm # clear cached value
                        if test "$dir" = "none"
                        then
                               try_loption=
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi

                        # Check first MKL...
                        FFLAGS="$test_fflags"
                        LDFLAGS="$MKL_FLAGS $test_ldflags $try_loption"
                        LIBS="-pgf90libs"

                        if test "$use_openmp" -eq 0; then
                              AC_SEARCH_LIBS(dgemm, mkl_intel_lp64,
                                 have_mathlib=1 have_mkl=1
                                 math_libs="$try_loption $LIBS -lmkl_sequential -lmkl_core"
                                 ldflags="$MKL_FLAGS $ldflags",
                                 echo "MKL not found",
                                 -lmkl_sequential -lmkl_core -ldl)
                        else
                              AC_SEARCH_LIBS(dgemm, mkl_intel_lp64,
                                 have_mathlib=1 have_mkl=1
                                 math_libs="$try_loption $LIBS -lmkl_core -lmkl_pgi_thread"
                                 ldflags="$MKL_FLAGS $ldflags",
                                 echo "MKL not found",
                                 -lmkl_sequential -lmkl_core -ldl -lpthread -lm)
                        fi

                        if test "$ac_cv_search_dgemm" != "no"
                        then break ; fi

                done 

                if test "$have_mkl" -eq 0
                then
                   for dir in $try_libdirs
                   do
                        unset ac_cv_search_dgemm # clear cached value
                        if test "$dir" = "none"
                        then
                               try_loption=
                        else
                                echo $ECHO_N "in $dir: " $ECHO_C
                                try_loption="-L$dir"
                        fi

                        FFLAGS="$test_fflags"
                        LDFLAGS=" $test_ldflags $try_loption"
                        LIBS=""

                        # Check PGI BLAS/LAPACK (if BLAS is present, there must be LAPACK as well)
                        unset ac_cv_search_dgemm # clear cached value
                        AC_SEARCH_LIBS(dgemm, blas,                                                                                                                                                                                                                       have_mathlib=1                                                                                                                                                                                                                                math_libs="$try_loption -lblas -llapack")

                        if test "$ac_cv_search_dgemm" != "no"
                        then break ; fi
                   done
                fi

                ;;
        esac
        
   else

        # blas provided in math_libs - not checked!
        have_mathlib=1
   fi
fi

# No blas/lapack library found or NETLIB esplicitly required

if test "$have_mathlib" -eq 0 -o "$use_netlib" -eq 1 ; then
    math_libs="\$(TOPDIR)/LAPACK/liblapack.a \$(TOPDIR)/LAPACK/libblas.a"
    netlib_libs_switch="internal"
else
    netlib_libs_switch="external"
fi

  mathlib_line="MATH_LIBS=$math_libs"

  AC_SUBST(math_libs)
  AC_SUBST(mathlib_line)

  AC_SUBST(netlib_libs_switch)

  AC_CONFIG_FILES(install/make_lapack.inc)
  
  AC_SUBST(math_libs) 
]
)
