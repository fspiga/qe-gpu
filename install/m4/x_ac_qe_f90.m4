# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_F90], [

# check Fortran compiler flags
# have_cpp=0: use external C preprocessing for fortran code -- not supported
# have_cpp=1: use C-like preprocessing in fortran compiler
have_pgi=0
have_cpp=1
xlf_flags=0

echo using F90... $f90

case "$arch:$f90_flavor" in
x86_64:pgf* )
        try_fflags_nomain="-Mnomain"
        try_fflags="-O3 -Mpreprocess"
        try_fflags_openmp="-mp"
        try_f90flags="-O3 -Mpreprocess -Mcache_align"
        try_fflags_noopt="-O0"
        try_ldflags=""
        try_ldflags_openmp="-mp"
        try_dflags="$try_dflags -D__PGI"
        have_cpp=0
        have_pgi=1
        ;;
crayxt*:pgf* )
        try_fflags_nomain="-Mnomain"
        try_fflags="-O3 -Mpreprocess"
        try_fflags_openmp="-mp"
        try_f90flags="-O3 -Mpreprocess -Mcache_align"
        try_fflags_noopt="-O0"
        try_ldflags=""
        try_ldflags_openmp="-mp"
        try_dflags="$try_dflags -D__PGI"
        have_cpp=0
        have_pgi=1
        ;;
ppc64*:xlf* )
        # Need review
        try_f90flags="\$(FFLAGS) -qfree=f90"
        try_fflags_noopt="-q64 -qthreaded -O0"
        try_ldflags="-q64 -qthreaded"
        pre_fdflags="-WF,"
        xlf_flags=1
        ;;
* )
				AC_MSG_ERROR([*** Compiler not properly detected!])
;;

esac

# Constrain true for pre-production release
if test "$have_pgi" -eq 0
then
  AC_MSG_ERROR([*** PGI compiler not detected!])
fi

if test "$use_openmp" -eq 1 ; then
  try_f90flags="$try_f90flags $try_fflags_openmp"
  try_fflags="$try_fflags $try_fflags_openmp"
  try_ldflags="$try_ldflags $try_ldflags_openmp"
fi

if test "$fflags" = ""   ; then fflags=$try_fflags     ; fi
if test "$f90flags" = "" ; then f90flags=$try_f90flags ; fi
if test "$fflags_noopt" = ""   ; then fflags_noopt=$try_fflags_noopt     ; fi
if test "$fflags_nomain" = ""   ; then fflags_nomain=$try_fflags_nomain     ; fi

echo setting FFLAGS... $fflags
echo setting F90FLAGS... $f90flags
echo setting FFLAGS_NOOPT... $fflags_noopt
if test "$fflags_nomain" != "" ; then echo setting FFLAGS_NOMAIN... $fflags_nomain ; fi

if test "$imod" = "" ; then imod="-I" ; fi

# compilation flags for all subsequent tests
# remove all $(...) because at least one compiler doesn't like them
# but if f90flags contains $(FFLAGS), substitute it
if test "`echo $f90flags | grep '$(FFLAGS)'`" != ""
then
        test_fflags="`echo $fflags $f90flags | sed 's/\$([[^)]]*)//g'`"
else
        test_fflags="`echo $f90flags | sed 's/\$([[^)]]*)//g'`"
fi

AC_SUBST(pre_fdflags)
AC_SUBST(f90flags)
AC_SUBST(fflags)
AC_SUBST(fflags_noopt)
AC_SUBST(fflags_nomain)
AC_SUBST(imod)

])
