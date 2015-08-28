#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#      script:  cts_gen.s 
#     
#     version: 0.0.1 
#     
#        date: 13 April 2013
#     
#      author:  Timothy J. Dennis
#     
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 
# Description:
#
# This script is designed to help the user prepare pbs batch script for
# post-processing. As of this writing it's an early verion which 
# puts something together for the script 'cts_idl.s' which is for making
# contour plots and animations of the fields
#
# The script is interactive and asks the user to input the following
# information:
#
# 1.) the field to be plot. Allowed values are: 'p', 'a', 'j', 'bz', 'vz'
#
# 2.) The descriptive label
#
# 3.) The nominal number of contours to use
#
# 4.) the first step 
#
# 5.) the last step
#
# 6.) resolution string
#
# 7.) total number of steps to base contours levels on
#
# 
# If the script's book-keeping efforts have been sufficiently successful 
# it will then execute the script "rungen-ff.s" which produces the batch
# files for the job chain. It then exits.
# 
# Note: at any time, the directory can be cleared of the "history"  
#       produced by run_init.s if one runs its companion script 
#       "revert.s."
#    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#!/bin/bash

this_host=`echo $HOSTNAME`
this_host=${this_host:0:4}
echo "this_host = " $this_host

   pp_queue="shared"                            # queue for post-processing
   echo -n "Which field?:"
   read fld
   echo -n "Run label?:"
   read desc_label 
   echo -n "How many contours?:"
   read  n_cnts
   echo -n "first step:"
   read first_step
   echo -n "last step:"
   read last_step
   echo -n "x-y resolution:"
   read x_y_res
   echo -n "z-resolution:"
   read z_res
   res_label=$x_y_res'_'$z_res
   echo "res_label = " $res_label
   echo -n "total steps:"
   read tot_steps
   echo -n "slice:"
   read slc
   echo -n "run wall-time:"
   read max_run_wall

   n_nodes=1
   ppn=1

   echo "#!/bin/bash" > gen.tmp
   echo "#PBS -W group_list=mhdturb" >> gen.tmp
   echo "#PBS -q " $pp_queue >> gen.tmp 
   echo "#PBS -l nodes="$n_nodes":ppn="$ppn >> gen.tmp 
   echo "#PBS -l walltime="$max_run_wall >> gen.tmp 
   echo "#PBS -r n" >> gen.tmp 
   echo " " >> gen.tmp
   echo "cd \$PBS_O_WORKDIR" >> gen.tmp
   echo "./cts_idl.s" $fld $desc_label $n_cnts $first_step $last_step $res_label  $tot_steps $slc >> gen.tmp 

  mv gen.tmp cts.pbs

exit 0 
