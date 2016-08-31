#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.63])
AC_INIT([coronos], [0.01], [tdennis10@alaska.edu])
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([src/coronos.cpp])
AC_CONFIG_HEADERS([config.h])

# Get the architecture
# -> if ARCH is 64 are lib64's gauranteed to exist? 

AC_MSG_CHECKING([architecture])

                   ARCH=`uname -m`

if test x$ARCH==xx86_64
then
                 SUFFIX="64"
else
                 SUFFIX=""
fi

AC_MSG_RESULT([$ARCH])

# Establish precision to be used 

AH_TEMPLATE([LD_PRECISION_H], [Define to one if using long double precision])
AH_TEMPLATE([OD_PRECISION_H], [Define to one if using ordinary double precision])

AC_ARG_WITH([long_doubles], AS_HELP_STRING([--with-long_doubles],
                                   [use long double precision (default is no)])],
                                   [],[])

if test -n "$with_long_doubles"
then
  AC_DEFINE([LD_PRECISION_H])
else
  AC_DEFINE([OD_PRECISION_H])
fi

    # ->->-> compiler info
    
  CRS_CXX="mpiCC CC gcc"
  AC_PROG_CXX([$CRS_CXX])
  AC_MSG_RESULT([$CRS_CXX])
  CRS_CC=" mpicc cc  gcc"
  AC_PROG_CC([$CRS_CC])
  AC_MSG_RESULT([$CRS_CC])

# Check for CUDA

AH_TEMPLATE([     HAVE_CUDA_H], [Define to one if using cuda])
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

                 FFTW_H="cufft.h"
               CUFFTW_H="cufftw.h"
              FFTW_LIBS=" -lcufft -lcufftw"

 # ->-> BLAS:

                 BLAS_H="cublas.h"
              BLAS_LIBS=" -lcublas"
           
                CRS_CXX="mpiCC CC"
                 CRS_CC=" mpicc cc "

  AC_MSG_CHECKING([header $CUFFTW_H ])
  AC_CHECK_HEADERS([$FFTW_PATH/include/$CUFFTW_H],[],[],[])

else # -> IF NOT-CUDA ==>

  AC_DEFINE([DONT_HAVE_CUDA_H])

  # ->-> Establish library search variables for headers and paths

    # ->->-> CUDA:

               CUDA_H=""
            CUDA_LIBS=""
    CUDA_PATH_DEFAULT=""

   
    # ->->-> FFT:

                    FFTW_H="fftw3.h"
                 FFTW_LIBS="-lfftw3"

         FFTW_PATH_DEFAULT="/usr/local/pkg/fftw"

fi # -> END CUDA vs. NOT-CUDA CONDITIONAL

# Checks for programs.

# Checks for libraries.

# Checks for header files.

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
    AC_MSG_NOTICE([BLAS_PATH=$BLAS_PATH])
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

               FFTW_INC="$FFTW_PATH"

          CRS_INC_FLAGS=" -I.. -I$FFTW_INC "

if test ! -n "$CUDA_PATH"
then
              CUDA_PATH=$CUDA_PATH_DEFAULT
fi

if test  -n "$CUDA_PATH"
then
                   NVCC="nvcc"
                 INCMPI="/opt/cray/mpt/5.6.5/gni/mpich2-cray/74/include"

               CUDA_INC="$CUDA_PATH/include"
          CRS_INC_FLAGS="-I$CUDA_INC -I$INCMPI $CRS_INC_FLAGS"
           CUDA_LDFLAGS="-L$CUDA_PATH/lib$SUFFIX"
           BLAS_LDFLAGS="-L$BLAS_PATH/lib$SUFFIX"
           FFTW_LDFLAGS="-L$FFTW_PATH/lib$SUFFIX"
else
               CUDA_INC=""
           CUDA_LDFLAGS=""

#-> something in here about finding where mpi.h is

           FFTW_LDFLAGS="-L$FFTW_PATH"
fi

               CRS_LIBS="$CUDA_LIBS $FFTW_LIBS "
           CRS_CPPFLAGS="$CRS_INC_FLAGS"
            CRS_LDFLAGS="$CUDA_LDFLAGS $FFTW_LDFLAGS"

echo "CUDA_PATH:    " $CUDA_PATH
echo "FFTW_PATH:    " $FFTW_PATH
echo "BLAS_PATH:    " $BLAS_PATH

echo "CUDA_LDFLAGS: " $CUDA_LDFLAGS
echo "FFTW_LDFLAGS: " $FFTW_LDFLAGS
echo "BLAS_LDFLAGS: " $BLAS_LDFLAGS

echo "CRS_LIBS:     " $CRS_LIBS
echo "CRS_CPPFLAGS: " $CRS_CPPFLAGS
echo "CRS_LDFLAGS:  " $CRS_LDFLAGS


AC_SUBST(INCMPI)
AC_SUBST(CRS_CPPFLAGS)
AC_SUBST(CRS_LIBS)
AC_SUBST(CRS_LDFLAGS)
AC_SUBST(CRS_NVCC)

AM_CONDITIONAL([CUDA],[test -n "$with_cuda"])

AC_CONFIG_FILES([Makefile
                 src/Makefile])
AC_OUTPUT
