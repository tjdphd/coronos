#!/bin/csh

##################################################################################
#                                                                                #
# c-shell script: glb_ext_idl.s - Timothy J. Dennis 13 June 2016                 #
#                             tdennis@gi.alaska.edu                              #
#                                                                                #
# Developed for use with the family of Heliospheric Physics codes originating    #
# with the work of C.S. Ng and including 'rmct_np' and 'fourfields.'             #
#                                                                                #
# Description: This script is designed to be executed from within a pbs script   #
#              though it can be run from the command line if the user knows      #
#              which arguments to pass to it, and their proper ordering.         #
#              This is, however, not recommended for jobs requiring large        #
#              numbers of plots - particularly when the data is of high          #
#              resolution.                                                       #
#                                                                                #
#              It's function is to prepare a compute node's environment for      #
#              submission of  the batch script 'sheet_batch' to an instantiation #
#              of idl on a compute node and allows for the generation of...      #
#                                                                                #
# Usage Example:                                                                 #
#                                                                                #
#                                                                                #
#                                                                                #
##################################################################################
#
# Arguments to 'glb_ext.s' are ultimately to be passed to 'glb_ext_batch'. This is
# accomplished by defining corresponding environment variables which may
# subsequently be accessed by IDL via the function 'GETENV'
#
############################
                           #
setenv GLEX_DESC_LABEL $1  #
setenv GLEX_RES_STR    $2  #
setenv GLEX_QTY        $3  #
setenv GLEX_FIRST_STP  $4  # - First time-step
setenv GLEX_LAST_STP   $5  # - Last  time-step
setenv GLEX_FIRST_SLC  $6  # - First slice
setenv GLEX_LAST_SLC   $7  # - Last slice
                           #
############################
#
# Next we load the idl module, tell it where to find IDL sources and
# tell it the name of the batch script we want it to execute
#
#######################################
                                      #
module load idl                       # load idl onto the compute node
                                      #
set       SRCDIR=$PWD/idl/glb_ext     # idl script source directory
echo   "SRCDIR   = " $SRCDIR          #
setenv IDL_PATH $IDL_DIR/lib:$SRCDIR  #
echo "IDL_PATH   = " $IDL_PATH        #
setenv IDL_STARTUP $PWD/glb_ext_batch # idl batch script to execute
                                      #
idl                                   # invoke idl
                                      #
#######################################
#
# Now clean things up.
#
###########################
                          #
module unload idl         #
unsetenv IDL_STARTUP      #
unsetenv GLEX_DESC_LABEL  #
unsetenv GLEX_RES_STR     #
unsetenv GLEX_QTY         #
unsetenv GLEX_FIRST_STP   #
unsetenv GLEX_LAST_STP    #
unsetenv GLEX_FIRST_SLC   #
unsetenv GLEX_LAST_SLC    #
                          #
exit 0                    #
                          #
###########################
##
###
#### The End ####
