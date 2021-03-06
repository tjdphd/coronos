#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
#

AC_PREREQ([2.63])
AC_INIT([coronos], [0.01], [tdennis10@alaska.edu])
AC_LANG(C++)
AM_INIT_AUTOMAKE([foreign])
AC_CONFIG_SRCDIR([src/coronos.cpp])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

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
AC_MSG_NOTICE([setting suffix to $SUFFIX])

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

                  USR="/usr"
            USR_LOCAL="/usr/local"
                  OPT="/opt"

  # ->->-> compiler info
    
  CRS_CXX="mpic++ mpiCC CC "
  AC_PROG_CXX([$CRS_CXX])
  CRS_CC="mpicc mpiCC cc "
  AC_PROG_CC([$CRS_CC])

  AX_COMPILER_VERSION

  AC_MSG_NOTICE([found compiler vendor $ax_cv_cxx_compiler_vendor])
  AC_MSG_NOTICE([found compiler version $ax_cv_cxx_compiler_version])

  PORTLAND=`expr match "$ax_cv_cxx_compiler_vendor" '\(portland\)'`
     INTEL=`expr match "$ax_cv_cxx_compiler_vendor" '\(intel\)'`
       GNU=`expr match "$ax_cv_cxx_compiler_vendor" '\(gnu\)'`

  if  test -n "$PORTLAND"
  then VENDOR=$PORTLAND
  fi
  if  test -n "$INTEL"
  then VENDOR=$INTEL
  fi
  if  test -n "$GNU"
  then VENDOR=$GNU
  fi

  echo "PORTLAND = " $PORTLAND
  echo "INTEL    = " $INTEL
  echo "GNU      = " $GNU
  echo "VENDOR   = " $VENDOR

# Find mpi.h

   if  test -n "$MPIHOME"
   then
     MPI_DIR=$MPIHOME
   fi

   if test -n "MPI_HOME"
   then
     MPI_DIR=$MPI_HOME
   fi
     MPI_INC="$MPI_DIR/include"
     MPI_LIB="$MPI_DIR/lib"
   
   if test -n "$MPI_DIR"
   then
     AC_CHECK_HEADERS([$MPI_INC/"mpi.h"],[],[],[])
   fi
   

  if [[ $VENDOR == $INTEL ]]
  then
    echo "VENDOR really is $VENDOR"
  fi

# Check for CUDA

AH_TEMPLATE([DONT_HAVE_CUDA_H], [Define to one if not using cuda])
AH_TEMPLATE(     [HAVE_CUDA_H], [Define to one if using cuda])

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

    FFTW_PATH_DEFAULT="$USR_LOCAL/fftw"
               FFTW_H="fftw3.h"
            FFTW_LIBS="-lm -lfftw3"


  if test -n "$FFTWHOME"
  then
  FFTW_PATH=$FFTWHOME
  fi
  if test -n "$FFTW_ROOT"
  then
  FFTW_PATH=$FFTW_ROOT
  fi

  AC_MSG_CHECKING([whether FFTW_PATH is set])

  if test -n "$FFTW_PATH"
  then
    AC_MSG_RESULT([yes])
    AC_MSG_NOTICE([FFTW_PATH = $FFTW_PATH])
    FFTW_PATH_INC=$FFTW_PATH/include
    FFTW_PATH_LIB_SFX="$FFTW_PATH/lib$SUFFIX"
    FFTW_PATH_LIB="$FFTW_PATH/lib"
    if test -e $FFTW_PATH_LIB_SFX
    then
    FFTW_PATH_LIB="$FFTW_PATH_SFX"
    fi
    echo "FFTW_PATH_LIB = " $FFTW_PATH_LIB
  else
    AC_MSG_RESULT([no])
    AC_MSG_WARN([FFTW_PATH not set])

    AC_MSG_CHECKING([if FFTW is available])
    FFTW_DIR_NUM=`ls -1d $OPT/fftw* 2>/dev/null | wc -l`
    echo "FFTW_DIR_NUM = " $FFTW_DIR_NUM
    if [[ $FFTW_DIR_NUM -eq '1' ]]
    then
    FFTW_ROOT_DIR=`ls -1d $OPT/fftw* 2>/dev/null`
    elif [[ $FFTW_DIR_NUM -eq '0' ]]
    then
      FFTW_DIR_NUM=`ls -1d $USR_LOCAL/fftw* 2>/dev/null | wc -l`
      if [[ $FFTW_DIR_NUM -eq '0' ]]
      then
        USR_LOC_PKG=`ls -1d /usr/local/pkg 2>/dev/null`
        if test -n $USR_LOC_PKG
        then
          echo "USR_LOC_PKG = " $USR_LOC_PKG
          FFTW_DIR_NUM=`ls -1d $USR_LOC_PKG/fftw* 2>/dev/null | wc -l`
          if [[ "$FFTW_DIR_NUM" -eq "0" ]]
          then
            echo "FFTW_DIR_NUM = " $FFTW_DIR_NUM
            AC_MSG_ERROR([after valient effort, giving up seeking FFTW_PATH. Please set manually],[1])
          else 
            if [[ "$FFTW_DIR_NUM" -eq "1" ]]
            then
              FFTW_ROOT_DIR=`ls -1d $USR_LOC_PKG/fftw*`
            else
              AC_MSG_WARN([found multiple FFTW directories, proceeding with caution, but you may want to set FFTW_PATH manually ])
              FFTW_ROOT_DIR=`ls -1d $USR_LOCAL 2>/dev/null`
            fi
          fi
        fi
      fi
    else
      AC_MSG_WARN([found multiple FFTW directories, proceeding with caution, but you may want to set FFTW_PATH manually ])
      FFTW_ROOT_DIR=`ls -1d $_OPT 2>/dev/null`
    fi    # should have FFTW_ROOT_DIR by now still need to set FFTW_PATH

    echo "FFTW_ROOT_DIR = " $FFTW_ROOT_DIR

    if test -r $FFTW_ROOT_DIR
    then
      echo "looking for FFTW_PATH"
      FFTW_INC_FOUND=`ls -1d $FFTW_ROOT_DIR/include 2>/dev/null | wc -l`
      if [[ "$FFTW_INC_FOUND" -eq "1" ]] 
      then
        echo "found an include file"
        FFTW_PATH=$FFTW_ROOT_DIR
      else
        echo "looking for subdirs starting with fftw"
        FFTW_DIR_NUM=`ls -1d $FFTW_ROOT_DIR/fftw* 2>/dev/null | wc -l`
        if [[ "$FFTW_DIR_NUM" -eq "0" ]] 
        then
          echo "did not find fftw-prefixed files looking more broadly"
          FFTW_ANYTHING_IN_HERE=`ls -1d $FFTW_ROOT_DIR/* 2>/dev/null | wc -l`
          if test -n $FFTW_ANYTHING_IN_HERE
          then
            echo "found something, looking for default"
            FFTW_DEFAULT_FOUND=`ls -1d $FFTW_ROOT_DIR/default 2>/dev/null | wc -l`
            if [[ "$FFTW_DEFAULT_FOUND" -eq "0" ]]
            then
              echo "giving up, nothing inside default"
              AC_MSG_ERROR([I am very sorry, but I do not know what is going on. You will have to set FFTW_PATH manually and try again.],[1])
            else 
              echo "found default directory resetting fftw root directory looking for interlagos"
              FFTW_ROOT_DIR=$FFTW_ROOT_DIR/default
              if test -n "$LOADEDMODULES"
              then
                INTERLAGOS_LOADED=`echo $LOADEDMODULES | grep -i interlagos`
                echo "INTERLAGOS_LOADED = " $INTERLAGOS_LOADED
                if test -n $INTERLAGOS_LOADED
                then
                  echo "found interlagos setting FFTW_PATH"
                  FFTW_PATH=$FFTW_ROOT_DIR/interlagos
                else
                  echo "did not find interagos"
                  AC_MSG_ERROR([I am very sorry, but I do not know what is going on. You will have to set FFTW_PATH manually and try again.],[1])
                fi
              else
                AC_MSG_ERROR([I am very sorry, but I do not know what is going on. You will have to set FFTW_PATH manually and try again.],[1])
              fi
            fi
          else
                echo "at line 261 should not quit here"
            AC_MSG_ERROR([FFTW_ROOT_DIR appears to be empty. Consider setting FFTW_PATH manually],[1])
          fi
        elif [[ "$FFTW_DIR_NUM" -eq "1" ]]
        then
          FFTW_ROOT_DIR=`ls -1d $FFTW_ROOT_DIR/fftw*`
          if test -r $FFTW_ROOT_DIR
          then
            FFTW_INC_FOUND=`ls -1d $FFTW_ROOT_DIR/include`
            if [[ "$FFTW_INC_FOUND" -eq "1" ]] 
            then
              FFTW_PATH=$FFTW_ROOT_DIR
            else
              AC_MSG_ERROR([after valient effort, giving up seeking FFTW_PATH. Please set manually],[1])
            fi
          else
            FFTW_PATH=$FFTW_ROOT_DIR
            AC_MSG_WARN([Cannot read $FFTW_ROOT_DIR. Assigning to FFTW_PATH anyway, you may want to manually set FFTW_PATH])
          fi
        else
          if test -n "$LOADEDMODULES" 
          then
             PGI_LOADED=`echo $LOADEDMODULES | grep -i pgi`
           INTEL_LOADED=`echo $LOADEDMODULES | grep -i intel`
             GNU_LOADED=`echo $LOADEDMODULES | grep -i gnu`
          fi

# ~       FOR intel compilers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

          if test -n "$INTEL"
          then

            echo "INTEL_LOADED=" $INTEL_LOADED

            if test -n "$PGI_LOADED" || test -n "$GNU_LOADED"
            then
               AC_MSG_WARN([multiple loaded programming enviroment modules detected. 
                            If configure fails unload redundant modules and try again.])
            fi

            
            INTEL_VER=$ax_cv_c_compiler_version

            FFTW_INTEL_VER_NUM=`ls -1d $FFTW_ROOT_DIR/fftw-*intel* 2>/dev/null |wc -l`
                        

            echo "INTEL_VER          = " $INTEL_VER
            echo "FFTW_INTEL_VER_NUM = " $FFTW_INTEL_VER_NUM
            echo "FFTW_INTEL_VER     = " $FFTW_INTEL_VER
            echo "FFTW_ROOT_DIR      = " $FFTW_ROOT_DIR

            if test -n $FFTW_INTEL_VER
            then
              FFTW_PATH=$FFTW_INTEL_VER
            fi

          fi

# ~       FOR portland compilers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

          if test -n "$PORTLAND"
          then

            echo "PGI_LOADED = " $PGI_LOADED

            if test -n "$INTEL_LOADED" || test -n "$GNU_LOADED"
            then
               AC_MSG_WARN([multiple loaded programming enviroment modules detected. 
                            If configure fails unload redundant modules and try again.])
            fi

            if test -n "$PGI"
            then
              echo "PGI        = " $PGI
              PGI_VER=${PGI##"$USR_LOC_PKG/pgi/"}
               echo "PGI_VER = " $PGI_VER

              if test -n "$with_long_doubles"
              then
                PGI_VER=".precision."$PGI_VER
              fi

              if test -n "$PGI_VER"
              then
                FFTW_PGI_VER=`ls -1d $FFTW_ROOT_DIR/fftw-*.?.$PGI_VER* 2>/dev/null`
                
                if test -n "$FFTW_PGI_VER" 
                then
                  if test -e $FFTW_PGI_VER
                  then
                    FFTW_PATH=$FFTW_PGI_VER
                  else
                    FFTW_PATH=$FFTW_ROOT_DIR
                    AC_MSG_WARN([the directory $FFTW_PGI_VER does not exist!])
                  fi
                fi
              fi
            fi
          fi

# ~       FOR GNU compilers ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

          if test -n "$GNU"
          then

            echo "GNU_LOADED = " $GNU_LOADED

            if test -n "$INTEL_LOADED" || test -n "$PGI_LOADED"
            then
               AC_MSG_WARN([multiple loaded programming enviroment modules detected. 
                            If configure fails unload redundant modules and try again.])
            fi

#                GNU_VER="gnu-"$ax_cv_c_compiler_version
                 GNU_VER="gnu"
            FFTW_GNU_VER=`ls -1d $FFTW_ROOT_DIR/fftw-*$GNU_VER* 2>/dev/null`
            echo "FFTW_GNU_VER" $FFTW_GNU_VER
            if test -n $FFTW_GNU_VER 
            then
              if test -e $FFTW_GNU_VER
              then
                FFTW_PATH=$FFTW_GNU_VER
              else
                AC_MSG_WARN([the directory $FFTW_GNU_VER does not exist!])
              fi
            fi


          fi

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

        fi
      fi
    else
      FFTW_PATH=$FFTW_ROOT_DIR
      AC_MSG_WARN([Cannot read $FFTW_ROOT_DIR. Assigning to FFTW_PATH anyway, you may want to manually set FFTW_PATH])
    fi

    AC_MSG_NOTICE([FFTW_PATH set to $FFTW_PATH])

    FFTW_PATH_INC="$FFTW_PATH/include"
    FFTW_PATH_LIB_SFX="$FFTW_PATH/lib$SUFFIX"
    FFTW_PATH_LIB="$FFTW_PATH/lib"
    echo "FFTW_PATH_LIB = " $FFTW_PATH_LIB
    if test -e $FFTW_PATH_LIB_SFX
    then
    FFTW_PATH_LIB="$FFTW_PATH_SFX"
    fi

  fi

fi # -> END CUDA vs. NOT-CUDA CONDITIONAL

# Checks for programs.

# Checks for libraries.

# Checks for header files.

AC_CHECK_HEADERS([$FFTW_PATH_INC/$FFTW_H],[
AC_MSG_NOTICE([please disregard the superfluous warning above regarding rejection by the pre-processor.])
],[],[])

# -> BLAS:

# AC_MSG_CHECKING([whether BLAS_PATH is set])

# -> -> BLAS headers

# if test -n "$BLAS_PATH"
# then
#   AC_MSG_RESULT([yes])
#   AC_MSG_NOTICE([BLAS_PATH=$BLAS_PATH])
# else
#   AC_MSG_RESULT([no])
#             BLAS_PATH=$BLAS_PATH_DEFAULT
#   AC_MSG_WARN([BLAS_PATH not set, using default: $BLAS_PATH])
# fi

# AC_CHECK_HEADERS([$BLAS_PATH/include/$BLAS_H],[],[],[])

# -> -> BLAS libs

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

# -> linear algebra:

# -> fft's

# Wrap up


          CRS_INC_FLAGS=" -I.. -I$FFTW_PATH_INC -I$MPI_INC"

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
           FFTW_LDFLAGS="-L$FFTW_PATH_LIB"
else
               CUDA_INC=""
           CUDA_LDFLAGS=""

#-> something in here about finding where mpi.h is

           FFTW_LDFLAGS="-L$FFTW_PATH_LIB"
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
