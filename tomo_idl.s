#!/bin/csh

################################################################################
#                                                                              #
# c-shell script: tomo_idl.s - Timothy J. Dennis 19 March 2014                 #
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
#              submission of the batch script 'tomo_batch' to an instantiation #
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
#              ./tomo_idl.s p fa 31 0 99                                       #
#                                                                              #
#              Note: the argument values given here must be consistent with    #
#                    similar arguments chosen in the scripts used when         #
#                    initiating the run that produced the data to be post-     #
#                    processed.                                                #
#                                                                              #
################################################################################
#
# Arguments to 'tomo_idl.s' are ultimately to be passed to 'tomo_batch'. This is
# accomplished by defining corresponding environment variables which may
# subsequently be accessed by IDL via the function 'GETENV'
#
#############################
                            #
setenv TOMO_CONT_FLD    $1  # - contoured field - 'p', 'a', 'bz', 'vz', or 'j'
setenv TOMO_DESC_LABEL  $2  # - string coding for information about the data
setenv TOMO_NUM_CONT    $3  # - Nominal number of contours to include in plots
setenv TOMO_FIRST_SLICE $4  # - First data-set in series to be contoured
setenv TOMO_LAST_SLICE  $5  # - Last data-set in series to be contoured
setenv TOMO_RES_STR     $6  # - String indicating run-resolution
setenv TOMO_TOT_SLCS    $7  # - Total number of steps upon which to base contour
                            #   levels
setenv TOMO_STP         $8  # - Which time-step is to be "tomographically" contoured
                            #
#############################
#
# Next we need to be sure that the appropriate sub-directory
# that is assumed to exist by the batch script will in fact
# exist.
#
##############################################
                                             #
if( ! -d tomo) then                          # upon first run the output
  mkdir tomo                                 # directory will not yet exist
  mkdir tomo/$TOMO_CONT_FLD                  # so we create it and whichever
  mkdir tomo/$TOMO_CONT_FLD/eps              # sub-directories not-yet 
else if( ! -d tomo/$TOMO_CONT_FLD) then      # created so we have place for
  mkdir tomo/$TOMO_CONT_FLD                  # for eps output
  mkdir tomo/$TOMO_CONT_FLD/eps              #
else if( ! -d tomo/$TOMO_CONT_FLD/eps ) then #
  mkdir tomo/$TOMO_CONT_FLD/eps              #
endif                                        #
                                             #
##############################################
#
# Next we load the idl module, tell it where to find IDL sources and
# tell it the name of the batch script we want it to execute
#
######################################
                                     #
#source /etc/profile.d/modules.csh    #
module load idl                      # load idl onto the compute node
                                     #
set SRCDIR=$PWD/idl/ff_tomo          # idl script source directory
echo "SRCDIR = " $SRCDIR             #
setenv IDL_PATH $IDL_DIR/lib:$SRCDIR #
echo "IDL_PATH = " $IDL_PATH         #
setenv IDL_STARTUP $PWD/tomo_batch   # idl batch script to execute
                                     #
idl                                  # invoke idl
                                     #
######################################
#
# Next initiate animation script
#
########################################################################
                                                                       #
./tomo_anim.s $TOMO_CONT_FLD $TOMO_RES_STR $TOMO_DESC_LABEL $TOMO_STP  #
                                                                       #
########################################################################
#
# Now clean things up.
#
####################################
                                   #
module unload idl                  #
unsetenv IDL_STARTUP               #
unsetenv TOMO_CONT_FLD             #
unsetenv TOMO_DESC_LABEL           #
unsetenv TOMO_NUM_CONT             #
unsetenv TOMO_FIRST_SLICE          #
unsetenv TOMO_LAST_SLICE           #
unsetenv TOMO_RES_STR              #
unsetenv TOMO_TOT_SLCS             #
unsetenv TOMO_STP                  #
                                   #
exit 0                             #
                                   #
####################################
##
###
#### The End ####
