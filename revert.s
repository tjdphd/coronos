#!/bin/bash

evaluate() { 

         line=${line#$seek}                      # chop p1 off of line
          val=`echo $line | awk '{print $1}'`    # the line now starts with the value of p1
}

clean=$1

if [ -z $clean ]
then
echo "reverting without clean..."
elif  [ "$clean" == "clean" ]
then
make clean
fi

FILE="coronos.in"

k=1
while read line                                  # Loop over lines of coronos.in
do
  for seek in "nprofile" "bdrys" "prefix" "run_label" "p1" "p2" "p3" "np" "data_dir" \
              "nnodes" "ppn" "srun" "calcqvz" "calcsvz" "qout_pref" "spout_pref"
  do
    if [ `expr "$line" : $seek` -ne 0 ]          # find line containing current seek string
    then
      evaluate                                   # val now contains the value of current seek string
       case "$seek" in                           # set parameters appropriately
         "prefix"    )       prefix=$val ;;
         "run_label" )    run_label=$val ;;
         "p1"        )           p1=$val ;;
         "p2"        )           p2=$val ;;
         "p3"        )           p3=$val ;;
         "np"        )           np=$val ;;
         "nnodes"    )       nnodes=$val ;;
         "ppn"       )          ppn=$val ;;
         "calcqvz"   )      calcqvz=$val ;;
         "calcsvz"   )      calcsvz=$val ;;
         "qout_pref" )    qout_pref=$val ;;
         "spout_pref")   spout_pref=$val ;;
         "data_dir"  )     data_dir=$val
       esac
    fi
  done
  ((k++))
done < $FILE
echo p3 $p3
echo np $np

xres=$((2**$p1))                                 # calculate resolution in x
yres=$((2**$p2))                                 # calculate resolution in y
zres=$(($p3*$np))                                # calculate resolution in z
echo zres $zres

if [[ "$xres" -eq "$yres" ]]
then
  res_label="$xres"_"$zres"
else
  res_label="$xres"_"$yres"_"$zres"
fi


        pbs_job_prefix=$prefix"-"$res_label"-"
        slm_job_prefix=$prefix"-"$run_label"-"$res_label"-sr-"
           data_prefix="./"$data_dir"/"$prefix"_"$res_label"."
          edata_prefix="./"$data_dir"/"$prefix"_"$res_label".o"$run_label
           rand_prefix=$prefix"_"$res_label"r"
  q_vs_z_tracking_file="./"$data_dir"/"$qout_pref'_'$res_label'.o'$run_label
 sp_vs_z_tracking_file="./"$data_dir"/"$spout_pref'_'$res_label'.???_???.o'$run_label

#rm  $job_prefix*.pbs
rm  $pbs_job_prefix*.pbs.[oe]*
rm  $pbs_job_prefix*.pbs
rm  $slm_job_prefix*.slurm
rm  $data_prefix???.o$run_label?
rm  $data_prefix???.o$run_label??
rm  $data_prefix???.o$run_label???
rm  $data_prefix???.o$run_label*.gz
rm  $data_prefix??.o$run_label?
rm  $data_prefix??.o$run_label??
rm  $data_prefix??.o$run_label???
rm  $edata_prefix
rm  $rand_prefix?
rm  $rand_prefix??
rm  $rand_prefix???
rm  $q_vs_z_tracking_file?
rm  $q_vs_z_tracking_file??
rm  $q_vs_z_tracking_file???
rm  $sp_vs_z_tracking_file?
rm  $sp_vs_z_tracking_file??
rm  $sp_vs_z_tracking_file???
