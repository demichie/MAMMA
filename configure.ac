#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(MAMMA,1.0,mattia.demichielivitturi@ingv.it)
AM_INIT_AUTOMAKE(1.0)

AC_CONFIG_SRCDIR(src/c_main.f90)

# Checks for programs.
AC_PROG_FC

# Checks for libraries.



# files to generate via autotools (prepare .am or .in source files)
AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])


AC_PREFIX_DEFAULT([..])

# finally this generates the Makefiles etc. for the build
AC_OUTPUT
