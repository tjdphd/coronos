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
#                                                                              #
################################################################################
#
###########################
                          #
setenv CC_QTY        $1   # - quantity to be compared 
setenv CC_FIRST_STEP $2   # - first output step to compare
setenv CC_LAST_STEP  $3   # - last  output step to compare
setenv CC_FIRST_SLC  $4   # - first output slice to compare
setenv CC_LAST_SLC   $5   # - last output slice to compare
setenv CC_RES_STR    $6   # - string indicating resolution
setenv CC_LABEL_ONE  $7   # - label for first data set
setenv CC_LABEL_TWO  $8   # - label for second data set
setenv CC_TOL        $9   # - tolerance
                          #   
###########################
#
######################################
                                     #
#source /etc/profile.d/modules.csh   #
module load idl                      # load idl onto the compute node
                                     #
set SRCDIR=$PWD/idl/pro/ff_cc        # idl script source directory
echo "SRCDIR   = " $SRCDIR           #
setenv IDL_PATH $IDL_DIR/lib:$SRCDIR #
echo "IDL_PATH = " $IDL_PATH         #
setenv IDL_STARTUP $PWD/idl/bin/cc_batch     # idl batch script to execute
                                     #
idl                                  # invoke idl
                                     #
######################################
