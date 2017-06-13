# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_MPIF90], [

AC_ARG_ENABLE(parallel,
   [AS_HELP_STRING([--enable-parallel],
       [compile for parallel execution if possible (default: yes)])],
   [set_use_parallel=1
    if   test "$enableval" = "yes" ; then
      use_parallel=1
   else
      use_parallel=0
   fi],
   [set_use_parallel=0 use_parallel=1])


# candidate fortran compilers good for all cases
try_mpif90="mpif90"
try_f90="gfortran f90"

# candidate compilers and flags based on architecture
case $arch in
x86_64 | ppc64 )
        try_f90="pgf90 $try_f90"
        ;;
cray )
        try_f90="ftn"
        try_mpif90="ftn"
        ;;
* )
        AC_MSG_WARN($arch : unsupported architecture?)
        ;;
esac

# check serial Fortran 90 compiler. This must be done before performing
# the check for the parallel compiler (section below) because option
# --disable-parallel will do very strange things otherwise. The reason
# seems to be that autoconf does not repeat all tests for the second
# occurrence of AC_PROG_FC. So the first occurrence is the one that
# must always be performed, the second is optional. PG & CC sep.2006

# use F90 if set
if test "$f90" = "" ; then f90="$try_f90" ; fi
#AC_PROG_FC([pgf90 ftn])
f90=$FC
AC_FC_SRCEXT(f90)

# check parallel Fortran 90 compiler
if test "$use_parallel" -eq 0 ;
then
        mpif90=$f90
else
        unset FC ac_cv_prog_ac_ct_FC ac_cv_fc_compiler_gnu ac_cv_prog_fc_g
        if test "$mpif90" = "" ; then 
	         mpif90="$try_mpif90 $f90"
           AC_PROG_FC($mpif90)
        else
           AC_PROG_FC($mpif90)
           if test "$FC" = "" ; then 
		          AC_MSG_WARN([MPIF90 not found: using MPIF90 anyway])
	  	        FC=$MPIF90
	         fi
        fi
        mpif90=$FC
fi

# check which compiler does mpif90 wrap
case "$arch" in
    x86_64 | ppc64 )
        echo $ECHO_N "checking version of $mpif90... $ECHO_C"
        pgf_version=`$mpif90 -V 2>&1 | grep "^pgf"`
        #
        if test "$pgf_version" != ""
        then
                version=`echo $pgf_version | cut -d ' ' -f2`
                echo "${ECHO_T}pgf90 $version"
                f90_in_mpif90="pgf90"
                # flag to test MKL with PGI
                MKL_FLAGS="-pgf90libs"
        else
                echo "${ECHO_T}unknown, assuming gfortran"
                f90_in_mpif90="gfortran"
        fi
        f90=$f90_in_mpif90
   ;;
esac

if test "$use_parallel" -eq 0
then
  AC_MSG_ERROR([*** MPI in required!])
fi


echo setting F90... $f90
echo setting MPIF90... $mpif90

case "$f90" in
f90 | fc | ftn )
    echo $ECHO_N "checking version wrapped by $f90 command... $ECHO_C"

    if $f90 -V 2>&1 | grep -q "^pgf" ; then
        f90_flavor=pgf
    else
        echo $ECHO_N "unknown, leaving as... $ECHO_C"
        f90_flavor=$f90
    fi
    echo $f90_flavor
    ;;
* )
    f90_flavor=$f90
    ;;
esac

AC_SUBST(f90)
AC_SUBST(mpif90)

])
