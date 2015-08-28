#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([coronos], [0.01], [tdennis@gi.alaska.edu])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([src/coronos.cpp])
AC_CONFIG_HEADERS([config.h])


# Get the architecture

AC_MSG_CHECKING([architecture])
ARCH=`uname -m`

if test x$ARCH == xx86_64
then
  SUFFIX="64"
else
  SUFFIX=""
fi

AC_MSG_RESULT([$ARCH])

# Check for CUDA

AH_TEMPLATE([HAVE_CUDA_H], [Define to one if using cuda])
AH_TEMPLATE([DONT_HAVE_CUDA_H], [Define to one if not using cuda])

AC_ARG_WITH([cuda], AS_HELP_STRING([--with-cuda=PATH],
                                   [use gpu acceleration (default is no)])],
                                   [],[])    

# Start CUDA vs. NOT-CUDA CONDITIONAL
# -> IF CUDA ==>

if test -n "$with_cuda" 
then 

  # -> first check if the request makes sense:

  # ->-> try nvcc --version

    NVCC=""
    AC_CHECK_PROG([NVCC],[nvcc],[nvcc],[""],[],[])

    echo "NVCC = " $NVCC
  # ->-> try whereis

    WHERECUDA=""
    if test ! -n "$NVCC"
    then
      WHERECUDA=`whereis cuda`
      WHERECUDA=${WHERECUDA#cuda:}
      echo "WHERECUDA" =$WHERECUDA
    fi

    USR_LOCAL_CUDA="/usr/local/cuda"

    if test -n "$WHERECUDA"
    then
      echo "CHAKRAH!"
      W_LEN=${#WHERECUDA}
      POS=`expr $W_LEN - 15`
      TEST_CUDA_DEFAULT=${WHERECUDA:$POS}
      echo "TEST_CUDA_DEFAULT = " $TEST_CUDA_DEFAULT
    else
      TEST_CUDA_DEFAULT=""
    fi
    
    if test x$TEST_CUDA_DEFAULT == x$USR_LOCAL_CUDA
    then
      TEST_CUDA_DEFAULT="$TEST_CUDA_DEFAULT/bin"
      echo "TEST_CUDA_DEFAULT = " $TEST_CUDA_DEFAULT
      NVCC=""
      AC_CHECK_PROG([NVCC],[nvcc],[nvcc],[""],[$TEST_CUDA_DEFAULT],[])
    fi
    
    if test ! -n "$NVCC"
    then
      AC_MSG_ERROR([cuda not found, check modules and try again])
    fi
    
  AC_DEFINE([HAVE_CUDA_H])


  # ->-> Establish library search variables for headers and paths

  # ->-> CUDA:

  CUDA_H="cuda.h"
  CUDA_LIBS=" -lcudart"

  # ->-> FFT:

  FFTW_H="cufft.h"
  FFTW_LIBS=" -lcufft"


  if test -n "$CRAY_CUDATOOLKIT_DIR"
  then
    CUDA_PATH_DEFAULT="$CRAY_CUDATOOLKIT_DIR"
    FFTW_PATH_DEFAULT="$CRAY_CUDATOOLKIT_DIR"
    BLAS_PATH_DEFAULT="$CRAY_CUDATOOLKIT_DIR"
  else
    CUDA_PATH_DEFAULT="$CRAY_CUDATOOLKIT_DIR"
    FFTW_PATH_DEFAULT="$CRAY_CUDATOOLKIT_DIR"
    BLAS_PATH_DEFAULT="$CRAY_CUDATOOLKIT_DIR"
  fi  

 # ->-> BLAS:

  BLAS_H="cublas.h"
  BLAS_LIBS=" -lcublas"
 
  CRS_CXX="mpiCC CC"
  CRS_CC=" mpicc cc "

else # -> IF NOT-CUDA ==>

  AC_DEFINE([DONT_HAVE_CUDA_H])


  # ->-> Establish library search variables for headers and paths

    # ->->-> CUDA:

    CUDA_H=""
    CUDA_LIBS=""
    CUDA_PATH_DEFAULT=""

    # ->->-> FFT:

    FFTW_H="fftw3.h"
    FFTW_LIBS="-lfftw3 -lm"

    if test -n "$CRAY_FFTW_DIR"
    then
      FFTW_PATH_DEFAULT="$CRAY_FFTW_DIR/interlagos"
    else
      FFTW_PATH_DEFAULT="/usr/local/pkg/fftw/fftw-3.3.3.pgi-13.4"
    fi

    # ->->-> BLAS 

    BLAS_H="gsl/gsl_blas.h"
    if test -n "$GSL_ROOT" 
    then
      BLAS_PATH_DEFAULT="$GSL_ROOT"
    else
      BLAS_PATH_DEFAULT="/usr/local/pkg/gsl/gsl-1.15.gnu-4.7.3"
    fi

    BLAS_LIBS="-lgsl -lgslcblas"

    CRS_CXX="mpiCC CC"
    CRS_CC="mpicc cc"

fi # -> END CUDA vs. NOT-CUDA CONDITIONAL

# Checks for programs.

# Checks for libraries.


# Checks for header files.


# -> FFT:

  AC_MSG_CHECKING([whether FFTW_PATH is set])

  if test -n "$FFTW_PATH"
  then
    AC_MSG_RESULT([yes])
    AC_MSG_NOTICE([FFTW_PATH = $FFTW_PATH])
  else
    AC_MSG_RESULT([no])
    FFTW_PATH=$FFTW_PATH_DEFAULT
    AC_MSG_WARN([FFTW_PATH not set, using default: $FFTW_PATH])
  fi

  AC_MSG_CHECKING([header $FFTW_H ])
  AC_CHECK_HEADERS([$FFTW_PATH/include/$FFTW_H],[],[],[])

# -> BLAS:

  AC_MSG_CHECKING([whether BLAS_PATH is set])

# -> -> BLAS headers


  if test -n "$BLAS_PATH"
  then
    AC_MSG_RESULT([yes])
    AC_MSG_NOTICE([BLAS_PATH = $BLAS_PATH])
  else
    AC_MSG_RESULT([no])
    BLAS_PATH=$BLAS_PATH_DEFAULT
    AC_MSG_WARN([BLAS_PATH not set, using default: $BLAS_PATH])
  fi

  AC_CHECK_HEADERS([$BLAS_PATH/include/$BLAS_H],[],[],[])

# -> -> BLAS libs

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# -> linear algebra:

# -> fft's


# Wrap up

FFTW_INC="$FFTW_PATH/include"
BLAS_INC="$BLAS_PATH/include"

CRS_INC_FLAGS=" -I$FFTW_INC -I$BLAS_INC"

if test ! -n "$CUDA_PATH"
then
  CUDA_PATH=$CUDA_PATH_DEFAULT
fi

if test  -n "$CUDA_PATH"
then
  CUDA_INC="$CUDA_PATH/include"
  CRS_INC_FLAGS="-I$CUDA_INC $CRS_INC_FLAGS"
  CUDA_LDFLAGS="-L$CUDA_PATH/lib$SUFFIX"
  BLAS_LDFLAGS="-L$BLAS_PATH/lib$SUFFIX"
  FFTW_LDFLAGS="-L$FFTW_PATH/lib$SUFFIX"
else
  CUDA_INC=""
  CUDA_LDFLAGS=""
  BLAS_LDFLAGS="-L$BLAS_PATH/lib"
  FFTW_LDFLAGS="-L$FFTW_PATH/lib"
fi

echo "CUDA_PATH:     " $CUDA_PATH
echo "FFTW_PATH:     " $FFTW_PATH
echo "BLAS_PATH:     " $BLAS_PATH

echo "CUDA_LDFLAGS:  " $CUDA_LDFLAGS
echo "FFTW_LDFLAGS:  " $FFTW_LDFLAGS
echo "BLAS_LDFLAGS:  " $BLAS_LDFLAGS

CRS_LIBS="$CUDA_LIBS $FFTW_LIBS $BLAS_LIBS"
echo "CRS_LIBS: " $CRS_LIBS

CRS_CPPFLAGS="$CRS_INC_FLAGS"
echo "CRS_CPPFLAGS: " $CRS_CPPFLAGS

CRS_LDFLAGS="$CUDA_LDFLAGS $FFTW_LDFLAGS $BLAS_LDFLAGS"

echo "CRS_LDFLAGS: " $CRS_LDFLAGS

NVCC="nvcc"

INCMPI="/opt/cray/mpt/5.6.5/gni/mpich2-cray/74/include"

AC_SUBST(INCMPI)

AC_SUBST(CRS_CPPFLAGS)
AC_SUBST(CRS_LIBS)
AC_SUBST(CRS_LDFLAGS)
AC_SUBST(CRS_NVCC)

AC_PROG_CXX([$CRS_CXX])
AC_PROG_CC([$CRS_CC])

AM_CONDITIONAL([CUDA],[test -n "$with_cuda"])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
