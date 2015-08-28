#!/bin/bash
 
################################################################################
#                                                                              #
# bash-shell script: cts_idl.s - Timothy J. Dennis 12 April 2013               #
#                             tdennis@gi.alaska.edu                            #
#                                                                              #
# Developed for use with the family of Heliospheric Physics codes originating  #
# with the work of C.S. Ng and including 'rmct_np' and 'fourfields.'           #
#                                                                              #
# Description: This script is designed to be executed from within a pbs script #
#              though it can be run from the command line if the user knows    #
#              which arguments to pass to it, and their proper ordering.       #
#              This is, however, not recommended for jobs requiring large      #
#              numbers of plots - particularly when the data is of high        #
#              resolution.                                                     #
#                                                                              #
#              It's function is to convert the eps files generated from the    #
#              coordinated execution of the c-shell script 'cts_idl.s', the    #
#              IDL batch script 'cts_batch', and the IDL procedure script -    #
#              and its dependencies - referenced therein to a set of gif       #
#              files. These gif files are then used to construct an animation  #
#              showing the time-evolution of a "horizontal" layer of the       #
#              solution space of the code from which the data was obtained.    #
#                                                                              #
#              As currently written the code can handle output for the fields  #
#              listed below:                                                   #
#                                                                              #
#              Field                                 cont_fld  Model Eq.'s     #
#                                                                              #
#              phi - 'Stream Function'               'p'       RMHD/FF         #
#              A   - 'Flux Function'                 'a'       RMHD/FF         #
#              J   - 'Current Density (z-comp)       'j'       RMHD/FF         #
#                                                                              #
#              Z   - 'B-field fluctations (z-comp)   'bz'      FF              #
#              V   - 'Velocity fluctuations (z-comp) 'vz'      FF              #
#                                                                              #
#                                                                              #
#              Also as currently written, this script should reside in a sub-  #
#              directory of the code run-directory called 'cts' though this    #
#              will likely change.                                             #
#                                                                              #
################################################################################


       pfx=$1 # field/quantity to be contoured - see above
   res_str=$2 # string giving the x (or y) X z-resolution of the data
desc_label=$3 # descriptive label for run
     n_slc=$4 # contoured layer

cur_dir=`pwd` # directory originating script invocation

cd cts

cts_eps=`ls $pfx/eps/$pfx_contour*stp-*ff-spec.eps`

if [ ! -e $pfx/gif ]   # create gif directory if needed
  then 
    mkdir $pfx/gif
fi

if [ ! -e $pfx/jpg ]   # create jpg directory if needed
  then 
    mkdir $pfx/jpg
fi

trunc=$pfx'/eps/'

for plot in $cts_eps; do
        plot=${plot#$trunc}
    len_plot=${#plot}                                    # Get the file's basename and create a gif version
           n=`expr $len_plot - 4`                        # with convert
   base_name=${plot:0:$n}
   echo "creating = " $base_name'.gif'
   convert -density 144x144 $pfx/eps/$base_name.eps $pfx/gif/$base_name.gif
done

cts_gif=`ls $pfx/gif/$1_contour*stp-*ff-spec.gif`

cts_anim=$1'_cts_'$res_str'_'$desc_label'-anim.gif'
convert -delay 10 -loop 0 -coalesce -background white -dispose 1 $cts_gif $cts_anim

for plot in $cts_eps; do
        plot=${plot#$trunc}
    len_plot=${#plot}                                    # Get the file's basename and create a gif version
           n=`expr $len_plot - 4`                        # with convert
   base_name=${plot:0:$n}
   echo "creating = " $base_name'.jpg'
   convert -density 144x144 $pfx/eps/$base_name.eps $pfx/jpg/$base_name.jpg
done


sep=`expr index $res_str '_'`
echo 'sep       = ' $sep

z_res=${res_str:$sep}
echo 'z_res     = ' $z_res

len_z_res=`expr length $z_res`
echo 'len_z_res = ' $len_z_res

len_n_slc=`expr length $n_slc`
echo 'len_n_slc = ' $len_n_slc

len_dif=`expr $len_z_res - $len_n_slc`

echo 'len_dif   = ' $len_dif

if [$len_dif -eq 0]
then
  zero_str=''
else
  zero_str='000000000000000000000000000000'
  pos=`expr 30 - $len_dif`
  echo 'pos = ' $pos
  zero_str=${zero_str:$pos}
fi

echo 'zero_str  = ' $zero_str

#if [ $len_n_slc -eq 1 ]
#then zero_str='00'
#fi
#if [ $len_n_slc -eq 2 ]
#then zero_str='0'
#fi
#if [ $len_n_slc -eq 2 ]
#then zero_str=''
#fi

str_n_slc=$zero_str$n_slc

echo "str_n_slc = " $str_n_slc

cts_jpg=$pfx/jpg/$1_contour_slc-${str_n_slc}_stp-%03d-ff-spec.jpg  # going to have to pass the slice as argument

cts_anim=$1'_cts_'$res_str'_'$desc_label'-anim.mp4'

if [ -e $cts_anim ]                                      # ffmpeg won't overwrite
then 
  rm $cts_anim
fi

ffmpeg -sameq -i $cts_jpg $cts_anim

exit 0
