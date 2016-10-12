#!/bin/bash
#
# function evaluate removes parameter name from line and reads parameter value #
# from remainder
#
evaluate() {                                     # e.g. - if $seek=p1...
         line=${line#$seek}                      # chop p1 off of line
          val=`echo $line | awk '{print $1}'`    # the line now starts with the value of p1
}
#
################################################################################
#
# create provisional archive sub-directory name based on name of current       #
# working directory:
#
# -> first obtain some necessary information from coronos.in
#
FILE="coronos.in"
   k=1
while read line                                  # Loop over lines of coronos.in
do
  for seek in "nprofile" "bdrys" "prefix" "run_label" "p1" "p2" "p3" "np"  \
              "nnodes" "ppn" "srun" "calcqvz" "calcsvz" "qout_pref" "spout_pref"
  do
    if [ `expr "$line" : $seek` -ne 0 ]
    then
      evaluate
       case "$seek" in
         "prefix"    )    file_pref=$val ;;
         "run_label" )    run_label=$val ;;
         "p1"        )           p1=$val ;;
         "p2"        )           p2=$val ;;
         "p3"        )           p3=$val ;;
         "np"        )           np=$val ;;
         "nnodes"    )       nnodes=$val ;;
         "bdrys"     )        bdrys=$val ;;
         "ppn"       )          ppn=$val ;;
         "srun"      )         srun=$val ;;
         "calcqvz"   )      calcqvz=$val ;;
         "calcsvz"   )      calcsvz=$val ;;
         "qout_pref" )    qout_pref=$val ;;
         "spout_pref")   spout_pref=$val ;;
       esac
     fi
   done
   ((k++))
done < $FILE
#
# -> next construct a resolution label
#
xres=$((2**$p1))                                 # calculate resolution in x
yres=$((2**$p2))                                 # calculate resolution in y
zres=$(($p3*$np))                                # calculate resolution in z
#

if [[ "$xres" -eq "$yres" ]]
then
  res_label="$xres"_"$zres"
else
  res_label="$xres"_"$yres"_"$zres"
fi
#
# -> now construct the archive directory name

           sep='/'
#         prefix='/import/c'                                       # Needed as a prefix to the archive pathname when
        prefix='/center'                                           # Needed as a prefix to the archive pathname when
  archive_root=$ARCHIVE                                            # $ARCHIVE is the environment variable specifying the
        whoami=`whoami`                                            # So is this. Find user's userid.
        center=$prefix'/w/'$whoami
        whoami=$sep$whoami$sep                                     # This too. Find the current directory.
current_wrkdir=$(echo $PWD)

  center_match=$(echo `expr match "$current_wrkdir" "$center"`)    # should be further developed and perhaps made
                                                                   # configurable so that the addition of known 
                                                                   # HPC file systems can be automated.

if [ $center_match -gt 0 ]                                         # If the HPC filesystem is $CENTER
then
  this_work_space=${current_wrkdir:0:$center_match}                # This should equal $CENTER.
     workdir_root=$center                                          # So we'll set workdir_root to $CENTER
                                                                   # and compare these strings later.
                                                                       
else                                                               # Or maybe we're not on a known HPC filesystem at all,

  echo "crs-archive.s - ERROR: the working directory " $workdir_root " is not recognized "
        exit 1
fi

if [[ ! "$this_work_space" == "$workdir_root" ]]                     # The call to 'match' can give an erroneous
  then                                                               # result because partial matches are possible.
                                                                     # So here is where we do a final check to make 
                                                                     # sure these two strings are actually equal.

    echo "crs-archive.s - ERROR: Your current directory does not appear to be working directory."
    echo "Exiting..."
    exit 1
fi

position=`expr match $current_wrkdir $workdir_root`                  # Determine the current directory's
  subdir=${current_wrkdir:position}                                  # path relative to root work directory
 
subdir_abspath=$archive_root$subdir                                  # Define the absolute path to the 
                                                                     # current directory's archive directory.
                                                                     # Note: "/import" has not been prefixed yet!

# -> we have the name, now check if it exists already and create if not
  
  if [ -e "$subdir_abspath" ]
  then
     if [ -d "$subdir_abspath" ]
     then
       if [ -w "$subdir_abspath" ]
       then
         echo "found writeable archive..." 
         echo $subdir_abspath
       fi
     else
       echo "crs-archive: ERROR - file $subdir_abspath exists but is not a writeable archive" 
       exit 1
     fi
  else
    mkdir -p $subdir_abspath
  fi 
#
# Now determine subruns to archive
#
if [ -z $2 ]
then
  args_abs=1
  if [ -z $1 ]
  then
        stop=`expr $srun - 2`
        if [ "$stop" -gt 0 ]
        then
          start=`expr $stop - 1`
        else
          start=$stop
        fi
  else
    start=$1
     stop=$start
  fi
else
  args_abs=0
     start=$1
      stop=$2
fi
#
# now gather file names
#
      field_data_files=$file_pref'_'$res_label'.???.o'$run_label
        sr_input_files=$file_pref'_'$res_label'.00.o'$run_label
  energy_tracking_file=$file_pref'_'$res_label'.o'$run_label
rand_num_tracking_file=$file_pref'_'$res_label'r'
  q_vs_z_tracking_file=$qout_pref'_'$res_label'.o'$run_label
 sp_vs_z_tracking_file=$spout_pref'_'$res_label'.???_???.o'$run_label

#
# now commence archiving unless there has been only one subrun
#
for j in `seq $start $stop`
do
  if [ "$j" -le "$srun" ]
  then

    if  [ "$j" -ge 0 ]
    then
      if    [[ "$args_abs" -eq 0 && "$j" -le "$srun" ]] || [[ "$args_abs" -eq 1 && "$j" -le ` expr $srun - 1 ` ]]
      then
#
# where the actual archiving gets done
#
        echo "archiving for subrun = " $j "..."

        cp coronos.in            $subdir_abspath
        cp $energy_tracking_file $subdir_abspath
        cp $sr_input_files$j     $subdir_abspath
        cp $field_data_files$j   $subdir_abspath
        if [ "$j" -gt 0 ]
        then
          if [ "$calcqvz" -eq 1 ]
          then
            cp $q_vs_z_tracking_file$j $subdir_abspath
          fi
#
          if [ "$calcsvz" -eq 1 ]
          then
            cp $sp_vs_z_tracking_file$j $subdir_abspath
          fi
#
          if [ "$bdrys" -ne 0 ]
          then
            cp $rand_num_tracking_file$j'.tgz' $subdir_abspath
          fi
        fi
        
      else
        echo "crs-archive: WARNING-refusing to archive data needed for next subrun"
        echo "                     if run is complete pass subrun numbers as arguments."
      fi
    fi

  else
    echo "crs-archive: WARNING-loop index exceeds number of subruns. Ignoring..."
    break
  fi
done

exit 0
