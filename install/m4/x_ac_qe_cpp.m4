# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_CPP], [

if test "$cpp" = "" ; then cpp=$try_cpp; fi
case $cpp in
        cpp)  try_cppflags="-P -traditional" ;;
        fpp)  try_cppflags="-P "              ;;
        *)    try_cppflags=""                ;;
esac
if test "$cppflags" = "" ; then cppflags=$try_cppflags ; fi

test_cppflags="$test_cflags"

AC_SUBST(cpp)
AC_SUBST(cppflags)
  
  ]
)
