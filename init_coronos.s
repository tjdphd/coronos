#!/bin/bash

evaluate() { 
         line=${line#$seek}                      # chop p1 off of line
          val=`echo $line | awk '{print $1}'`    # the line now starts with the value of p1
}

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

    last=$(($start-1))
  penult=$(($stop-1))

FILE="coronos.in"

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

xres=$((2**$p1))                                 # calculate resolution in x
yres=$((2**$p2))                                 # calculate resolution in y
zres=$(($p3*$np))                                # calculate resolution in z

echo "zres = " $zres
echo "p3   = " $p3
echo "np   = " $n3

if [[ "$xres" -eq "$yres" ]]
then
  res_label="$xres"_"$zres"
else
  res_label="$xres"_"$yres"_"$zres"
fi

         job_name=$prefix"-"$res_label"-"
        arch_name=$job_name$"archive-"

echo " "
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

if  [ "$ppn" -le 12 ] 
then
  cores=12
elif [ "$ppn" -gt 12 ]
then
  cores=16
fi

for j in `seq $start $stop`
do
#
   subr=$job_name$j".pbs"
subarch=$arch_name$(($j-1))".pbs"
#
echo "#!/bin/bash"                        > $subr
echo "#PBS -W group_list=mhdturb"        >> $subr
echo "#PBS -q standard_$cores"           >> $subr
echo "#PBS -l nodes=$nnodes:ppn=$ppn"    >> $subr
echo "#PBS -l walltime=$run_wall"        >> $subr
echo "#PBS -r n"                         >> $subr
echo " "                                 >> $subr
echo "cd \$PBS_O_WORKDIR"                >> $subr
echo "mpirun -np $np src/coronos"        >> $subr
#echo "qsub $subarch"                     >> $subr
#
 nextsr=$job_name$(($j+1))".pbs"
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
echo "qsub $nextsr"                   >> $subarch
fi
#
done

exit 0
