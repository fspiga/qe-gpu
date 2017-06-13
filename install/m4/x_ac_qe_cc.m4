# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_CC], [

# candidate C and f77 compilers good for all cases
try_cc="cc gcc"

case "$arch:$f90_flavor" in
*:pgf* )
        try_cc="pgcc $try_cc"
        ;;
esac

# check serial C compiler
if test "$env_cc" = "" ; then cc="$try_cc" ; else cc="$env_cc"; fi
AC_PROG_CC($cc)
cc=$CC

echo setting CC... $cc

AC_SUBST(cc)

# tentative C and loader flags, good for many cases
try_cflags="-O3"
c_ldflags=""
try_cpp="cpp"

case "$arch:$cc" in
*:pgcc )
        try_cflags="-fast -Mpreprocess"
        ;;
esac
if test "$cflags" = "" ; then cflags=$try_cflags ; fi
echo setting CFLAGS... $cflags

# compilation flags for all subsequent tests
test_cflags="`echo $cflags | sed 's/\$([[^)]]*)//g'`"

AC_SUBST(cflags)
])
