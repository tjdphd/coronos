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
         "data_dir"  )     data_dir=$val ;;
         "srun"      )         srun=$val
       esac
    fi
  done
  ((k++))
done < $FILE

xres=$((2**$p1))                                 # calculate resolution in x
yres=$((2**$p2))                                 # calculate resolution in y
zres=$(($p3*$np))                                # calculate resolution in z

if [[ "$xres" -eq "$yres" ]]
then
  res_label="$xres"_"$zres"
else
  res_label="$xres"_"$yres"_"$zres"
fi


        pbs_job_prefix=$prefix"-"$res_label"-"
        slm_job_prefix=$prefix"-"$run_label"-"$res_label"-sr-"
        slm_arc_prefix=$prefix"-"$run_label"-"$res_label"-ar-"
           data_prefix="./"$data_dir"/"$prefix"_"$res_label"."
          edata_prefix="./"$data_dir"/"$prefix"_"$res_label".o"$run_label
           rand_prefix="./"$data_dir"/"$prefix"_"$res_label"r"
  q_vs_z_tracking_file="./"$data_dir"/"$qout_pref'_'$res_label'.o'$run_label
 sp_vs_z_tracking_file="./"$data_dir"/"$spout_pref'_'$res_label'.???_???.o'$run_label

      pbs_in=$pbs_job_prefix*.pbs
     pbs_out=$pbs_job_prefix*.pbs.[oe]*
      slm_in=$slm_job_prefix*.slurm
      arc_in=$slm_arc_prefix*.slurm
     slm_out='tjd.'
     nds_out='nodes.'


fld_data_out=$data_prefix???.o$run_label
par_data_out=$data_prefix??.o$run_label
erg_data_out=$data_prefix??.o$run_label
rnd_data_out=$rand_prefix
qvz_data_out=$q_vs_z_tracking_file
svz_data_out=$sp_vs_z_tracking_file

if [ -e $pbs_in       ]
then 
  rm $pbs_in
  rm $pbs_out
fi

if [ -e $fld_data_out0 ]; then rm $fld_data_out?; fi
if [ -e $par_data_out0 ]; then rm $par_data_out?; fi
if [ -e $erg_data_out0 ]; then rm $erg_data_out?; fi
if [ -e $rnd_data_out0 ]; then rm $rnd_data_out?; fi
if [ -e $qvz_data_out0 ]; then rm $qvz_data_out?; fi
if [ -e $svz_data_out0 ]; then rm $svz_data_out?; fi

rm $pbs_in
rm $pbs_out
rm $slm_in
rm $arc_in
rm $slm_out*
rm $nds_out*

if (( "$srun" >= "10" ))
then

  if [ -e $fld_data_out10 ]; then rm $fld_data_out??; fi
  if [ -e $par_data_out10 ]; then rm $par_data_out??; fi
  if [ -e $erg_data_out10 ]; then rm $erg_data_out??; fi
  if [ -e $rnd_data_out10 ]; then rm $rnd_data_out??; fi
  if [ -e $qvz_data_out10 ]; then rm $qvz_data_out??; fi
  if [ -e $svz_data_out10 ]; then rm $svz_data_out??; fi

  if (( "$srun" >= "100" ))
  then

   if [ -e $fld_data_out100 ]; then rm $fld_data_out???; fi
   if [ -e $par_data_out100 ]; then rm $par_data_out???; fi
   if [ -e $erg_data_out100 ]; then rm $erg_data_out???; fi
   if [ -e $rnd_data_out100 ]; then rm $rnd_data_out???; fi
   if [ -e $qvz_data_out100 ]; then rm $qvz_data_out???; fi
   if [ -e $svz_data_out100 ]; then rm $svz_data_out???; fi

  fi
fi

cmp_data_out=$data_prefix???.o$run_label*.gz

if [ -e $edata_prefix ]; then rm  $edata_prefix; fi
if [ -e $cmp_data_out ]; then rm  $cmp_data_out; fi
if [ -e "crs_init.in" ]; then cp "crs_init.in" "coronos.in"; fi
