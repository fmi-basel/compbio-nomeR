## configure.ac is borrowed from https://www.r-bloggers.com/2017/08/setting-up-optional-openmp-support-with-rcpparmadillo/
AC_PREREQ(2.69)

AC_INIT(nomeR, m4_esyscmd_s([awk -e '/^Version:/ {print $2}' DESCRIPTION]))
AC_COPYRIGHT(Copyright (C) 2020 Evgeniy A. Ozonov)


# Set R_HOME, respecting an environment variable if one is set
: ${R_HOME=$(R RHOME)}
if test -z "${R_HOME}"; then
    AC_MSG_ERROR([Could not determine R_HOME.])
fi
# Use R to set CXX and CXXFLAGS
CXX=$(${R_HOME}/bin/R CMD config CXX)
CXXFLAGS=$("${R_HOME}/bin/R" CMD config CXXFLAGS)

# We are using C++
AC_LANG(C++)
AC_REQUIRE_CPP


dnl this the meat of R's m4/openmp.m4
  OPENMP_[]_AC_LANG_PREFIX[]FLAGS=
  AC_ARG_ENABLE([openmp],
    [AS_HELP_STRING([--disable-openmp], [do not use OpenMP])])
  if test "$enable_openmp" != no; then
    AC_CACHE_CHECK([for $[]_AC_CC[] option to support OpenMP],
      [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp],
      [AC_LINK_IFELSE([_AC_LANG_OPENMP],
     [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='none needed'],
     [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp='unsupported'
      for ac_option in -fopenmp -xopenmp -qopenmp \
                           -openmp -mp -omp -qsmp=omp -homp \
               -fopenmp=libomp \
                           -Popenmp --openmp; do
        ac_save_[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
        _AC_LANG_PREFIX[]FLAGS="$[]_AC_LANG_PREFIX[]FLAGS $ac_option"
        AC_LINK_IFELSE([_AC_LANG_OPENMP],
          [ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp=$ac_option])
        _AC_LANG_PREFIX[]FLAGS=$ac_save_[]_AC_LANG_PREFIX[]FLAGS
        if test "$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp" != unsupported; then
          break
        fi
      done])])
    case $ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp in #(
      "none needed" | unsupported)
    ;; #(
      *)
    OPENMP_[]_AC_LANG_PREFIX[]FLAGS=$ac_cv_prog_[]_AC_LANG_ABBREV[]_openmp ;;
    esac
  fi

AC_SUBST(OPENMP_CXXFLAGS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
