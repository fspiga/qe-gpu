# Copyright (C) 2001-2016 Quantum ESPRESSO Foundation

AC_DEFUN([X_AC_QE_ARCH], [

  AC_MSG_CHECKING([ARCH])

#cross_compiling=yes

if test "$host" = "" ; then host=$build; fi

# identify host architecture
if test "$arch" = ""
then
        case $host in
                86_64-*-linux-gnu )    arch=x86_64 ;;
                powerpc64le-* ) arch=ppc64  ;;
                * )                     AC_MSG_WARN(Unrecognized build architecture)
        ;;
        esac

        test -d /proc/cray_xt && arch=crayxt
fi

  AC_MSG_RESULT(${arch})
  AC_SUBST(arch)

])
