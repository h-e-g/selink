#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([selink], [1.0], [pierre.boutillier@laposte.net, epatin@pasteur.fr, glaval@pasteur.fr])
AC_CONFIG_SRCDIR([src/selink.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC_C99
AC_USE_SYSTEM_EXTENSIONS

# Checks for libraries.
AC_SEARCH_LIBS([log], [m])
AC_SEARCH_LIBS([pthread_create], [pthread])

# Checks for typedefs, structures, and compiler characteristics
AC_C_INLINE
AC_TYPE_SIZE_T
AC_SYS_LARGEFILE

# Checks for header files.
AC_CHECK_HEADERS([ctype.h err.h fcntl.h getopt.h libgen.h stdbool.h stdio.h sys/mman.h])
AC_CHECK_HEADERS([pthread.h])

# Checks for library functions.
AC_FUNC_STAT
AC_FUNC_MMAP
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([strtol strdup strchr pow asprintf])

AM_INIT_AUTOMAKE
AC_CONFIG_FILES([Makefile src/Makefile test/Makefile])
AC_OUTPUT
