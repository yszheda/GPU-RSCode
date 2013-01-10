# ------------------------------------------------------------------------------
AC_PREREQ([2.65])
AC_INIT([cuAutotools], [5.0.0], [zarnick@geekvault.org])
AM_INIT_AUTOMAKE([cuAutotools], [5.0.0])
# ------------------------------------------------------------------------------
modified by Shuai YUAN
# ------------------------------------------------------------------------------
AC_CONFIG_SRCDIR([src/matrix.cu])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AM_PROG_CC_C_O

# Checks for libraries.
AC_CHECK_LIB([m], [pow])
# ------------------------------------------------------------------------------
modified by Shuai YUAN
add libstdc++ lib to the linker command line
# ------------------------------------------------------------------------------
AC_CHECK_LIB(stdc++, main,,AC_MSG_ERROR(gdu requires libstdc++))

# Checks for header files.
AC_CHECK_HEADERS([limits.h stdlib.h string.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([pow])

AC_CONFIG_FILES([Makefile
                 src/Makefile])

#find out what version we are running
ARCH=`uname -m`
if [[ $ARCH == "x86_64" ]];
then
  SUFFIX="64"
else
  SUFFIX=""
fi

# Setup CUDA paths
# ------------------------------------------------------------------------------
AC_ARG_WITH([cuda],
   [  --with-cuda=PATH  		    prefix where cuda is installed [default=auto]])
if test -n "$with_cuda"
then
   CUDA_CFLAGS="-I$with_cuda/include"
   CUDA_LIBS="-L$with_cuda/lib$SUFFIX"
   CUDA_LDFLAGS="-L$with_cuda/lib$SUFFIX"
   NVCC="$with_cuda/bin/nvcc"
else
   CUDA_CFLAGS="-I/usr/local/cuda/include"
   CUDA_LIBS="-L/usr/local/cuda/lib$SUFFIX"
   CUDA_LDFLAGS="-L/usr/local/cuda/lib$SUFFIX"
   NVCC="nvcc"
fi
AC_ARG_ENABLE([debug],
              [--enable-debug               enable nVidia CUDA debug mode [default=no]],
              [NVCCFLAGS="-G -g"],
              [NVCCFLAGS=""]
              )
AC_SUBST(CUDA_CFLAGS)
AC_SUBST(CUDA_LIBS)
AC_SUBST(NVCC)
AC_SUBST(NVCCFLAGS)

#Check for CUDA libraries
save_LDFLAGS="$LDFLAGS"
LDFLAGS="$LDFLAGS $CUDA_LDFLAGS"
AC_CHECK_LIB([cudart], [cudaMalloc])
LDFLAGS="$save_LDFLAGS"

AC_OUTPUT
