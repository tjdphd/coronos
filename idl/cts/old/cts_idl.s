#!/bin/csh

################################################################################
#                                                                              #
# c-shell script: cts_idl.s - Timothy J. Dennis 12 April 2013                  #
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
#              It's function is to prepare a compute node's environment for    #
#              submission of  the batch script 'cts_batch' to an instantiation #
#              of idl on a compute node and allows for the generation of large #
#              numbers of contour plots based on high-resolution data to be    #
#              produced 'automatedly' so as to avoid the need to stress the    #
#              resources of a login node, or the patience of the code-jockey   #
#              (that would me or you) who needs said contour plots to move his #
#              research forward.                                               #
#                                                                              #
#              Note: though currently designed specifically for use with       #
#                    IDL contour-plotting scripts, it can and should be        #
#                    generalized to work with other scripts needed by the      #
#                    family of Heliospheric codes for which it was developed,  #
#                    but these things take time...                             #
#                                                                              #
# Usage Example:                                                               #
#                                                                              #
#              ./cts_idl.s p fa 31 0 99                                        #
#                                                                              #
#              Note: the argument values given here must be consistent with    #
#                    similar arguments chosen in the scripts used when         #
#                    initiating the run that produced the data to be post-     #
#                    processed.                                                #
#                                                                              #
################################################################################
#
# Arguments to 'cts_idl.s' are ultimately to be passed to 'cts_batch'. This is
# accomplished by defining corresponding environment variables which may
# subsequently be accessed by IDL via the function 'GETENV'
#
###########################
                          #
setenv CTS_CONT_FLD   $1  # - contoured field - 'p', 'a', 'bz', 'vz', or 'j'
setenv CTS_DESC_LABEL $2  # - string coding for information about the data
setenv CTS_NUM_CONT   $3  # - Nominal number of contours to include in plots
setenv CTS_FIRST_STEP $4  # - First data-set in series to be contoured
setenv CTS_LAST_STEP  $5  # - Last data-set in series to be contoured
setenv CTS_RES_STR    $6  # - String indicating run-resolution
setenv CTS_TOT_STPS   $7  # - Total number of steps upon which to base contour
                          #   levels
setenv CTS_SLC        $8  # - Which layer or "slice" to be contoured
                          #
###########################
#
# Next we need to be sure that the appropriate sub-directory
# that is assumed to exist by the batch script will in fact
# exist.
#
############################################
                                           #
if( ! -d cts) then                         # upon first run the output
  mkdir cts                                # directory will not yet exist
  mkdir cts/$CTS_CONT_FLD                  # so we create it and whichever
  mkdir cts/$CTS_CONT_FLD/eps              # sub-directories not-yet 
else if( ! -d cts/$CTS_CONT_FLD) then      # created so we have place for
  mkdir cts/$CTS_CONT_FLD                  # for eps output
  mkdir cts/$CTS_CONT_FLD/eps              # 
else if( ! -d cts/$CTS_CONT_FLD/eps ) then # 
  mkdir cts/$CTS_CONT_FLD/eps              # 
endif                                      #
                                           #
############################################
#
# Next we load the idl module, tell it where to find IDL sources and
# tell it the name of the batch script we want it to execute
#
######################################
                                     #
#source /etc/profile.d/modules.csh    #
module load idl                      # load idl onto the compute node
                                     #
set SRCDIR=$PWD/idl/ff_contour_ts    # idl script source directory
echo "SRCDIR = " $SRCDIR             #
setenv IDL_PATH $IDL_DIR/lib:$SRCDIR #
echo "IDL_PATH = " $IDL_PATH         #
setenv IDL_STARTUP $PWD/cts_batch    # idl batch script to execute
                                     #
idl                                  # invoke idl
                                     #
######################################
#
# Next initiate animation script
#
##################################################################
                                                                 #
./cts_anim.s $CTS_CONT_FLD $CTS_RES_STR $CTS_DESC_LABEL $CTS_SLC #
                                                                 #
##################################################################
#
# Now clean things up.
#
####################################
                                   #
module unload idl                  #
unsetenv IDL_STARTUP               #
unsetenv CTS_CONT_FLD              #
unsetenv CTS_DESC_LABEL            #
unsetenv CTS_NUM_CONT              #
unsetenv CTS_FIRST_STEP            #
unsetenv CTS_LAST_STEP             #
unsetenv CTS_RES_STR               #
unsetenv CTS_TOT_STPS              #
unsetenv CTS_SLC                   #
                                   #
exit 0                             #
                                   #
####################################
##
###
#### The End ####
