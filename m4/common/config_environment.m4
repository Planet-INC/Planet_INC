# SYNOPSIS
#
#   Queries configuration environment.
#
#   AX_SUMMARIZE_ENV([, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
#
# DESCRIPTION
#
#   Queries compile environment and SVN revision for use in configure summary 
#   and pre-processing macros.
#
# LAST MODIFICATION
#
#   19 Feb. 2014
#
# COPYLEFT
#
#   Copyright (c) 2014 Sylvain A. Plessis <splessis@ices.utexas.edu>
#   Copyright (c) 2010 Karl W. Schulz <karl@ices.utexas.edu>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_SUMMARIZE_ENV],
[

AC_CANONICAL_HOST

BUILD_USER=${USER}
BUILD_ARCH=${host}
BUILD_HOST=${ac_hostname}
BUILD_DATE=`date +'%F %H:%M'`

# Determine method for querying Source code revisioning (assumes SVN)

AC_PATH_PROG(svnquery,svnversion)
AC_PATH_PROG(gitquery,git)

# svnversion 1.7.x changed the return value for unversioned directories
if test "x${svnquery}" = "x" || test "`${svnquery} -n $srcdir`" = "exported" || test "`${svnquery} -n $srcdir`" = "Unversioned directory"; then
   SVN_REVISION="cat $srcdir/dist_version"
   SVN_CHECKOUT=false
   BUILD_DEVSTATUS="External Release"
else
   SVN_REVISION="${svnquery} -n $srcdir"
   SVN_CHECKOUT=true
   BUILD_DEVSTATUS="Development Build"
fi
if [test -d $srcdir/.git] ; then
   GIT_HASH="${gitquery} --git-dir=$srcdir/.git rev-parse HEAD"
   GIT_TAG="${gitquery} --git-dir=$srcdir/.git describe --abbrev=0"
   BUILD_DEVSTATUS="Development Build"
   GIT_CHECKOUT=true
else
   GIT_HASH="cat $srcdir/dist_version"
   GIT_TAG="External Release"
   BUILD_DEVSTATUS="External Release"
   GIT_CHECKOUT=false
fi

# Query current version.

if test "${SVN_CHECKOUT}" == "true"; then
  BUILD_VERSION=`${SVN_REVISION}`
else
  BUILD_VERSION=`${GIT_HASH}`
fi

AC_SUBST(SVN_REVISION)
AC_SUBST(BUILD_DEVSTATUS)
AM_CONDITIONAL(SVN_CHECKOUT,test x${SVN_CHECKOUT} = xtrue )
AC_SUBST(GIT_HASH)
AC_SUBST(GIT_TAG)
AM_CONDITIONAL(GIT_CHECKOUT,test x${GIT_CHECKOUT} = xtrue )


# Versioning info - check local developer version (if checked out)

AC_DEFINE_UNQUOTED([BUILD_USER],     "${BUILD_USER}",     [The fine user who built the package])
AC_DEFINE_UNQUOTED([BUILD_ARCH],     "${BUILD_ARCH}",     [Architecture of the build host])
AC_DEFINE_UNQUOTED([BUILD_HOST],     "${BUILD_HOST}",     [Build host name])
AC_DEFINE_UNQUOTED([BUILD_VERSION],  "${BUILD_VERSION}",  [SVN revision/git hash])
AC_DEFINE_UNQUOTED([BUILD_DEVSTATUS],"${BUILD_DEVSTATUS}",[Dev/Release build])
AC_DEFINE(         [BUILD_DATE],     __DATE__" "__TIME__, [Build date])

AC_SUBST(BUILD_USER)
AC_SUBST(BUILD_ARCH)
AC_SUBST(BUILD_HOST)
AC_SUBST(BUILD_DATE)
AC_SUBST(BUILD_VERSION)

])
