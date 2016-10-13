#!/bin/csh

##################################################################################
#                                                                                #
# c-shell script: jsheet_idl.s - Timothy J. Dennis 24 March 2014                 #
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
#              ./jstruct_idl.s ts, 1, 32, 1, 100                                 #
#                                                                                #
#                                                                                #
##################################################################################
#
# Arguments to 'jsheet_idl.s' are ultimately to be passed to 'jsheet_batch'. This is
# accomplished by defining corresponding environment variables which may
# subsequently be accessed by IDL via the function 'GETENV'
#
#######################################
                                      #
setenv JSHT_QTY        $1             #
setenv JSHT_DESC_LABEL $2             #
setenv JSHT_ACTION     $3             #
setenv JSHT_RES_STR    $4             #
setenv JSHT_THRSH      $5             #
setenv JSHT_FIRST_STP  $6             # - First time-step
setenv JSHT_LAST_STP   $7             # - Last  time-step
setenv JSHT_FIRST_SLC  $8             # - First slice
setenv JSHT_LAST_SLC   $9             # - Last slice
setenv JSHT_Y_ROT      $10            # - angle of rotation about y
setenv JSHT_N_CNTRS    $11            # - number of contours and/or isosurfaces
                                      #
if( -e glb_ext.out) then              # - if global extrema data is available
  setenv JSHT_GLB_EXT 'y'             #
endif                                 #
                                      #
#######################################
#
# Next we need to be sure that the appropriate sub-directory
# that is assumed to exist by the batch script will in fact
# exist.
#
############################################
                                           #
if( ! -d jstruct) then                     # upon first run the output
  mkdir jstruct                            # directory will not yet exist
  mkdir jstruct/eps                        # directory will not yet exist
  mkdir jstruct/gif                        # directory will not yet exist
  mkdir jstruct/jpg                        # directory will not yet exist
endif                                      #
                                           #
if( ! -d ostruct) then                     # upon first run the output
  mkdir ostruct                            # directory will not yet exist
  mkdir ostruct/eps                        # directory will not yet exist
  mkdir ostruct/gif                        # directory will not yet exist
  mkdir ostruct/jpg                        # directory will not yet exist
endif                                      #
                                           #
############################################
#
# Next we load the idl module, tell it where to find IDL sources and
# tell it the name of the batch script we want it to execute
#
#######################################
                                      #
module load idl                       # load idl onto the compute node
                                      #
set SRCDIR=$PWD/idl/jstruct           # idl script source directory
echo "SRCDIR = " $SRCDIR              #
setenv IDL_PATH $IDL_DIR/lib:$SRCDIR  #
echo "IDL_PATH = " $IDL_PATH          #
setenv IDL_STARTUP $PWD/jstruct_batch # idl batch script to execute
                                      #
idl                                   # invoke idl
                                      #
#######################################
#
# Next initiate animation script
#
########################################################################################################################################
                                                                                                                                       #
#./jstruct_anim.s  $JSHT_DESC_LABEL $JSHT_ACTION $JSHT_RES_STR $JSHT_THRSH $JSHT_FIRST_STP $JSHT_LAST_STP $JSHT_FIRST_SLC $JSHT_LAST_SLC #
                                                                                                                                       #
########################################################################################################################################
#
# Now clean things up.
#
####################################
                                   #
module unload idl                  #
unsetenv IDL_STARTUP               #
unsetenv JSHT_QTY                  #
unsetenv JSHT_DESC_LABEL           #
unsetenv JSHT_ACTION               #
unsetenv JSHT_RES_STR              #
unsetenv JSHT_THRSH                #
unsetenv JSHT_FIRST_STP            #
unsetenv JSHT_LAST_STP             #
unsetenv JSHT_FIRST_SLC            #
unsetenv JSHT_LAST_SLC             #
unsetenv JSHT_Y_ROT                #
unsetenv JSHT_N_CNTRS              # 
                                   #
exit 0                             #
                                   #
####################################
##
###
#### The End ####
