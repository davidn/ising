dnl Process this file with autoconf to produce a configure script.
dnl Created by Anjuta application wizard.

AC_INIT(ising, 0.1)

AC_DEFINE(OUTPUT_GNUPLOT,1,Output graph to GNUPLOT)
AC_DEFINE(DEFAULT_SIZE,15,Default lattice size)
AC_DEFINE(DEFAULT_kT,1,Default temperature)
AC_DEFINE(DEFAULT_J,1,Default interaction parameter)
AC_DEFINE(DEFAULT_muH,0,Default external field)

AM_INIT_AUTOMAKE(AC_PACKAGE_NAME, AC_PACKAGE_VERSION)
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
AM_MAINTAINER_MODE

AC_LANG([C++])

AC_ISC_POSIX
AM_PROG_CC_STDC
AC_HEADER_STDC
AX_PROG_CXX_CXX0X(, AC_MSG_ERROR(No C++0x compiler found))

AC_OUTPUT([
Makefile
src/Makefile
tests/Makefile
])
