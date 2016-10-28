#!/bin/csh

#################################################################################
#                                                                               #
# c-shell script: evsz_idl.s - Timothy J. Dennis 10 April 2014                  #
#                             tdennis@gi.alaska.edu                             #
#                                                                               #
# Developed for use with the family of Heliospheric Physics codes originating   #
# with the work of C.S. Ng and including 'rmct_np' and 'fourfields.'            #
#                                                                               #
# Description: This script is designed to be executed from within a pbs script  #
#              though it can be run from the command line if the user knows     #
#              which arguments to pass to it, and their proper ordering.        #
#              This is, however, not recommended for jobs requiring large       #
#              numbers of plots - particularly when the data is of high         #
#              resolution.                                                      #
#                                                                               #
#              It's function is to prepare a compute node's environment for     #
#              submission of  the batch script 'evsz_batch' to an instantiation #
#              of idl on a compute node and allows for the generation of large  #
#              numbers of plots based on high-resolution data to be             #
#              produced 'automatedly' so as to avoid the need to stress the     #
#              resources of a login node, or the patience of the code-jockey    #
#              (that would me or you) who needs said contour plots to move his  #
#              research forward.                                                #
#                                                                               #
# Usage Example:                                                                #
#                                                                               #
#              ./evsz_idl.s ts 0 99                                             #
#                                                                               #
#              Note: the argument values given here must be consistent with     #
#                    similar arguments chosen in the scripts used when          #
#                    initiating the run that produced the data to be post-      #
#                    processed.                                                 #
#                                                                               #
#################################################################################
#
# Arguments to 'cts_idl.s' are ultimately to be passed to 'cts_batch'. This is
# accomplished by defining corresponding environment variables which may
# subsequently be accessed by IDL via the function 'GETENV'
#
############################
                           #
setenv EVSZ_DSC_LAB     $1 # - string coding for information about the data
setenv EVSZ_RES_STR     $2 #
setenv EVSZ_EFLD        $3 # - string coding for information about the data
setenv EVSZ_FIRST_STEP  $4 # - First data-set in series to be plotted
setenv EVSZ_LAST_STEP   $5 # - Last data-set in series to be plotted
                           #
############################
#
#
# Next we need to be sure that the appropriate sub-directory
# that is assumed to exist by the batch script will in fact
# exist.
#
############################################
                                           #
if( ! -d ra_evsz) then                     # upon first run the output
  mkdir ra_evsz                            # directory will not yet exist
  mkdir -p ra_evsz/$EVSZ_EFLD/eps          # directory will not yet exist
else if( ! -d ra_evsz/$EVSZ_EFLD) then     #
  mkdir -p ra_evsz/$EVSZ_EFLD/eps          # directory will not yet exist
else if( ! -d ra_evsz/$EVSZ_EFLD/eps) then #
  mkdir ra_evsz/$EVSZ_EFLD/eps             # directory will not yet exist
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
set SRCDIR=$PWD/idl/pro/evsz          #
echo "SRCDIR = " $SRCDIR              #
setenv IDL_PATH $IDL_DIR/lib:$SRCDIR  #
echo "IDL_PATH = " $IDL_PATH          #
setenv IDL_STARTUP $PWD/idl/bin/ra_evsz_batch # idl batch script to execute
                                      #
idl                                   # invoke idl
                                      #
#######################################
#
# Now clean things up.
#
####################################
                                   #
module unload idl                  #
unsetenv IDL_STARTUP               #
unsetenv EVSZ_DSC_LAB              #
unsetenv EVSZ_RES_STR              #
unsetenv EVSZ_EFLD                 #
unsetenv EVSZ_FIRST_STEP           #
unsetenv EVSZ_LAST_STEP            #
                                   #
exit 0                             #
                                   #
####################################
##
###
#### The End ####
