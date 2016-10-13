#!/bin/csh
#
############################
                           #
setenv SPEC_DSC_LAB     $1 # - string coding for information about the data
setenv SPEC_RES_STR     $2 #
setenv SPEC_SFLD        $3 # - string coding for information about the data
setenv SPEC_LAYER       $4 # - layer to average over
setenv SPEC_FIRST_STEP  $5 # - First data-set in series to be plotted
setenv SPEC_LAST_STEP   $6 # - Last data-set in series to be plotted
setenv SPEC_MINMAX_MODE $7 # - either 'global' or 'local'
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
if( ! -d ra_spec) then                     # upon first run the output
  mkdir ra_spec                            # directory will not yet exist
  mkdir -p ra_spec/$SPEC_SFLD/eps          # directory will not yet exist
else if( ! -d ra_spec/$SPEC_SFLD) then     #
  mkdir -p ra_spec/$SPEC_SFLD/eps          # directory will not yet exist
else if( ! -d ra_spec/$SPEC_SFLD/eps) then #
  mkdir ra_spec/$SPEC_SFLD/eps             # directory will not yet exist
endif                                      #
                                           #
############################################
#
# Next we load the idl module, tell it where to find IDL sources and
# tell it the name of the batch script we want it to execute
#
######################################
                                     #
source /etc/profile.d/modules.csh    #
module load idl                      # load idl onto the compute node
                                     #
set SRCDIR=$PWD/idl/pro/spec             #
setenv IDL_PATH $IDL_DIR/lib:$SRCDIR #
echo "IDL_PATH = " $IDL_PATH         #
setenv IDL_STARTUP $PWD/idl/bin/ra_spec_batch   # idl batch script to execute
                                     #
idl                                  # invoke idl
                                     #
######################################
#
# Now clean things up.
#
####################################
                                   #
module unload idl                  #
unsetenv IDL_STARTUP               #
unsetenv SPEC_DSC_LAB              #
unsetenv SPEC_RES_STR              #
unsetenv SPEC_SFLD                 #
unsetenv SPEC_LAYER                #
unsetenv SPEC_FIRST_STEP           #
unsetenv SPEC_LAST_STEP            #
unsetenv SPEC_MINMAX_MODE          #
                                   #
exit 0                             #
                                   #
####################################
##
###
#### The End ####
