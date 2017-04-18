#!/bin/bash
#
evaluate() {                                     # e.g. - if $seek = p1...
         line=${line#$seek}                      # chop p1 off of line
          val=`echo $line | awk '{print $1}'`    # the line now starts with the value of p1
}
#
if [ -z $3 ]
then
 rsc_man=".pbs"
else
 rsc_man="."$3
fi
#
if [ -z  $2 ]
then 
    start=1
  if [ -z $1 ]
  then
     stop=1
  else
     stop=$1
  fi
else
    start=$1
     stop=$2
fi
#
    last=$(($start-1))
  penult=$(($stop-1))
#
FILE="coronos.in"
#
k=1
while read line                                  # Loop over lines of coronos.in
do
  for seek in "nprofile" "prefix" "run_label" "p1" "p2" "p3" "np" "nnodes" "ppn"
  do
    if [ `expr "$line" : $seek` -ne 0 ]          # find line containing current seek string
    then
      evaluate                                   # val now contains the value of current seek string
       case "$seek" in                           # set parameters appropriately
         "prefix"    )    prefix=$val ;;
         "run_label" ) run_label=$val ;;
         "p1"        )        p1=$val ;;
         "p2"        )        p2=$val ;;
         "p3"        )        p3=$val ;;
         "np"       )         np=$val ;;
         "nnodes"    )    nnodes=$val ;;
         "ppn"       )       ppn=$val ;;
       esac
    fi
  done
  ((k++))
done < $FILE
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
         job_name=$prefix"-"$run_label"-"$res_label"-sr-"
         arc_name=$prefix"-"$run_label"-"$res_label"-ar-"
#
echo " "
echo " job_name = " $job_name
echo " arc_name = " $arc_name
echo "Job Specifications:"
echo " "
echo "   First subrun:    $start"
echo "   Final subrun:    $stop"
echo "   X-resolution:    $xres"
echo "   Y-resolution:    $yres"
echo "   Z-resolution:    $zres"
echo "number of nodes:    $nnodes"
echo "processes per node: $ppn"
echo "total processes:    $np"
echo " "
echo -n "Does this look correct? (y/n):"
read ans
echo -n "specify run wall time (nn:nn:nn):"
read run_wall
echo -n "specify archive wall time (nn:nn:nn):"
read arch_wall
echo -n "specify partition:"
read partition
echo -n "specify project:"
read project
#
if  [ "$ppn" -le  16 ] 
then
  cores=16
elif [ "$ppn" -gt 16 ]
then
  cores=$ppn
fi
#
   str_zeros="00000000000000000000000000000000"
str_stop_len=${#stop}
if [[ "$str_stop_len" < 3 ]]
then
     min_len=3 
else
     min_len=$str_stop_len
fi
   add_zeros=${str_zeros:0:$min_len-1}
#
  cur_dir=`pwd`
#
     subr=$job_name$j_sr$rsc_man
  subarch=$arc_name$j_ar$rsc_man
#
if [[ "$rsc_man" == ".pbs" ]]
then
for j in `seq $start $stop`
do
  str_j_len=${#j}
  if [[ "$str_j_len"  < "$min_len" ]]
  then
    j_sr=$add_zeros$j
    j_ar=`expr $j - 1`
    j_ar=$add_zeros$j_ar
    j_nx=`expr $j + 1`
    j_nx=$add_zeros$j_nx
  fi
     subr=$job_name$j_sr$rsc_man
  subarch=$arc_name$j_ar$rsc_man
#
  echo "#!/bin/bash"                      > $subr
  echo "#PBS -W group_list=mhdturb"      >> $subr
  echo "#PBS -q standard_$cores"         >> $subr
  echo "#PBS -l nodes=$nnodes:ppn=$ppn"  >> $subr
  echo "#PBS -l walltime=$run_wall"      >> $subr
  echo "#PBS -r n"                       >> $subr
  echo " "                               >> $subr
  echo "cd \$PBS_O_WORKDIR"              >> $subr
  echo "mpirun -np $np src/coronos"      >> $subr
# echo "qsub $subarch"                   >> $subr
#
  nextsr=$job_name$j_nx$rsc_man
  echo "nextsr = " $nextsr
  echo "#!/bin/bash"                     > $subarch
  echo "#PBS -W group_list=mhdturb"     >> $subarch
  echo "#PBS -q transfer"               >> $subarch
  echo "#PBS -l nodes=1:ppn=1"          >> $subarch
  echo "#PBS -l walltime=$arch_wall"    >> $subarch
  echo "#PBS -r n"                      >> $subarch
  echo " "                              >> $subarch
  echo "cd \$PBS_O_WORKDIR"             >> $subarch
  echo "#./crs-archive.s"               >> $subarch
  if [[ "$j" -ne "$stop" ]]
  then
    echo "qsub $nextsr"                 >> $subarch
  fi
done
#
elif [[ "$rsc_man" == ".slurm" ]]
then
#
  jobr=$job_name$j_sr_$start-$stop$rsc_man
  joba=$arc_name$j_ar_$start-$stop$rsc_man
#
  echo "#!/bin/sh"                                                                            > $jobr
  echo " "                                                                                   >> $jobr
# echo "#SBATCH -A $project"                                                                 >> $jobr
  echo "#SBATCH --job-name=coronos"                                                          >> $jobr
  echo "#SBATCH --output=tjd.%j.%N.out"                                                      >> $jobr
  echo "#SBATCH --partition=$partition"                                                      >> $jobr
  echo "#SBATCH --nodes=$nnodes"                                                             >> $jobr
  echo "#SBATCH --tasks-per-node=$ppn"                                                       >> $jobr
  echo "#SBATCH --export=ALL"                                                                >> $jobr
  echo "#SBATCH -t $run_wall"                                                                >> $jobr
  echo "#SBATCH  --array=$start-$stop%1"                                                     >> $jobr       
  echo "#SBATCH --mail-user=tdennis10@alaska.edu"                                            >> $jobr
  echo "#SBATCH --mail-type=FAIL"                                                            >> $jobr
  echo " "                                                                                   >> $jobr
  echo ". /usr/share/Modules/init/sh"                                                        >> $jobr
  echo " "                                                                                   >> $jobr
  echo "cd $cur_dir"                                                                         >> $jobr
  echo "srun -l /bin/hostname | sort -n | awk '{print \$2}' > ./nodes.\$SLURM_ARRAY_TASK_ID" >> $jobr
  echo "time mpirun -np $np -machinefile ./nodes.\$SLURM_ARRAY_TASK_ID ./src/coronos"        >> $jobr
# echo "time ibrun -v ./src/coronos"                                                         >> $jobr
  echo " "                                                                                   >> $jobr
  echo "#EOF "                                                                               >> $jobr
  chmod u+x $jobr
#
  echo "#!/bin/sh"                                                                            > $joba
  echo " "                                                                                   >> $joba
  echo "#SBATCH --partition=transfer"                                                        >> $joba
  echo "#SBATCH --ntasks=1"                                                                  >> $joba
  echo "#SBATCH --tasks-per-node=1"                                                          >> $joba
  echo "#SBATCH --mail-user=tdennis10@alaska.edu"                                            >> $joba
  echo "#SBATCH --mail-type=FAIL"                                                            >> $joba
  echo "#SBATCH --time=$arch_wall"                                                           >> $joba
  echo "#SBATCH --output=tjd.%j"                                                             >> $joba
  echo " "                                                                                   >> $joba
  echo ". /usr/share/Modules/init/sh"                                                        >> $joba
  echo " "                                                                                   >> $joba
  echo "cd $cur_dir"                                                                         >> $joba
  echo "srun -l /bin/hostname | sort -n | awk '{print \$2}' > ./nodes.\$SLURM_JOB_ID"        >> $joba
  echo "time $cur_dir/crs-archive.s $start $stop"                                            >> $joba
  echo " "                                                                                   >> $joba
  echo "#EOF "                                                                               >> $joba
  chmod u+x $joba
else
  echo "init-coronos: ERROR - the resource management option " $rsc_man " is unknown."
fi

exit 0
