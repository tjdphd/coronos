#!/bin/csh

####################################################################################
#                                                                                  #
# c-shell script: ffts_idl.s - Timothy J. Dennis 18 October 2013                   #
#                             tdennis@gi.alaska.edu                                #
#                                                                                  #
# Developed for use with the family of Heliospheric Physics codes originating      #
# with the work of C.S. Ng and including 'rmct_np' and 'fourfields.'               #
#                                                                                  #
# Description: This script is designed to be executed from within a pbs script     #
#              though it can be run from the command line if the user knows        #
#              which arguments to pass to it, and their proper ordering.           #
#              This is, however, not recommended for jobs requiring large          #
#              numbers of plots - particularly when the data is of high            #
#              resolution.                                                         #
#                                                                                  #
#              It's function is to prepare a compute node's environment for        #
#              submission of  the batch script 'ffts_batch' to an instantiation    #
#              of idl on a compute node and allows for the generation of large     #
#              numbers of time-series plots based on high-resolution data to be    #
#              produced 'automatedly' so as to avoid the need to stress the        #
#              resources of a login node, or the patience of the code-jockey       #
#              (that would me or you) who needs said time-series plots to move his #
#              research forward.                                                   #
#                                                                                  #
#              Note: though currently designed specifically for use with           #
#                    IDL plotting scripts, it can and should be                    #
#                    generalized to work with other scripts needed by the          #
#                    family of Heliospheric codes for which it was developed,      #
#                    but these things take time...                                 #
#                                                                                  #
# Usage Examples:                                                                  #
#                                                                                  #
#              ./ffts_idl.s       Will prompt for descriptive label and field      #
#              ./ffts_idl.s ts    Will promt for field to plot                     #
#              ./ffts_idl.s ts 2  Plots field 2 for run labeled 'ts'               #
#                                                                                  #
#              Note: the argument values given here must be consistent with        #
#                    similar arguments chosen in the scripts used when             #
#                    initiating the run that produced the data to be post-         #
#                    processed.                                                    #
#                                                                                  #
####################################################################################
#
# Arguments to 'ffts_idl.s' are ultimately to be passed to 'ffts_batch'. This is
# accomplished by defining corresponding environment variables which may
# subsequently be accessed by IDL via the function 'GETENV'
#
###################################################################
                                                                  #
set dsc_lab='empty'                                               #
                                                                  #
if ($#argv < 1) then                                              #
                                                                  #
  echo -n "Descriptive Label: "                                   #
  set dsc_lab=$<                                                  #
                                                                  #
endif                                                             #
                                                                  #
if ($#argv < 2) then                                              #
                                                                  #
  if ($#argv > 0) then                                            #
                                                                  #
      set dsc_lab=$1                                              #
      echo "dsc_lab = " $dsc_lab                                  #
                                                                  #
  endif                                                           #
                                                                  #
  echo " "                                                        #
  echo "Time-Series Plotting Options For Run " $dsc_lab " :"      #
  echo " "                                                        #
  echo "  1 - Perpendicular Kinetic Energy                     "  #
  echo "  2 - Perpendicular Magnetic Energy                    "  #
  echo "  3 - Maximum Vorticity                                "  #
  echo "  4 - Location in z for maximum vorticity              "  #
  echo "  5 - Location in plain maximum vorticity              "  #
  echo "  6 - Square magnitude of vorticity Laplacian          "  #
  echo "  7 - Maximum Current                                  "  #
  echo "  8 - Location in z for maximum current                "  #
  echo "  9 - Location in plain for maximum current            "  #
  echo " 10 - Square magnitude of current Laplacian            "  #
  echo " 11 - Perpendicular viscous dissipation                "  #
  echo " 12 - Perpendicular resistive dissipation              "  #
  echo " 13 - Poynting flux                                    "  #
  echo " 14 - Average injected foot-point energy per unit time "  #
  echo " 15 - Average dissipated energy per unit time          "  #
  echo " 16 - Average rate-of-change of total energy           "  #
  echo " 17 - Instantaneous rate-of-change of total energy     "  #
  echo " 18 - Average rate of energy gain minus energy loss    "  #
  echo " 19 - consE                                            "  #
  echo " 20 - Current time in units of correlation time        "  #
  echo " 21 - dt                                               "  #
  echo " 22 - dtvb                                             "  #
  echo " 23 - sqrt(ce/cee) - whatever that is                  "  #
  echo " 24 - aveK2/t                                          "  #
  echo " 25 - sqrt(aveME/t)                                    "  #
  echo " 26 - sqrt(2.*AVEpe/t)                                 "  #
  echo " 27 - Footpoint Energy (fp)                            "  #
  echo " 28 - Internal Energy (he)                             "  #
  echo " 29 - Parallel Kinetic Energy                          "  #
  echo " 30 - Parallel Magnetic Energy                         "  #
  echo " 31 - Internal Energy                                  "  #
  echo " 32 - Square magnitude of Z Laplacian                  "  #
  echo " 33 - Square magnitude of V_z Laplacian                "  #
  echo " 34 - Parallel Conductive heat loss                    "  #
  echo " 35 - Parallel viscous dissipation                     "  #
  echo " "                                                        #
  echo -n "field: "                                               # - interactive bits
                                                                  #
  set fld_one=$<                                                  #
                                                                  #
  echo -n "Second Field?[y/n]:"                                   #
  set fld_two=$<                                                  #
                                                                  #
  if ($fld_two == 'y') then                                       #
                                                                  #
    echo -n "field two: "                                         #
    set fld_two=$<                                                #
                                                                  #
  else                                                            #
                                                                  #
    set fld_two='0'                                               #
                                                                  #
  endif                                                           #
                                                                  #
else if ($#argv >= 2) then                                        #
                                                                  #
  set dsc_lab=$1                                                  #
  set fld_one=$2                                                  #
                                                                  #
  if ($#argv > 2) then                                            #
                                                                  #
    set fld_two=$3                                                #
  else                                                            #
    set fld_two='0'                                               #
  endif                                                           #
                                                                  #
endif                                                             #
                                                                  #
if ($fld_two == '0') then                                         #
                                                                  #
  set lmt_fld=$fld_one                                            #
                                                                  #
else                                                              #
                                                                  #
  echo " "                                                        #
  echo "Limit Options: "                                          #
  echo " "                                                        #
  echo "  1 - use a field "                                       #
  echo "  2 - manual "                                            #
  echo " "                                                        #
  echo -n "Select limit option: "                                 #
  set lim_opt=$<                                                  #
                                                                  #
  if ($#argv > 3) then                                            #
                                                                  #
    set lmt_fld=$4                                                #
                                                                  #
  else                                                            #
                                                                  #
    echo -n "Limiting field: "                                    #
    set lmt_fld=$<                                                #
                                                                  #
  endif                                                           #
                                                                  #
endif                                                             #
                                                                  #
echo " "                                                          #
echo "dsc_lab: " $dsc_lab                                         #
echo "fld_one: " $fld_one                                         #
echo "fld_two: " $fld_two                                         #
echo "lmt_fld: " $lmt_fld                                         #
echo " "                                                          #
                                                                  #
###################################################################
#
####################################
                                   #
  setenv FFTS_DSC_LAB    $dsc_lab  # - descriptive label for run
  setenv FFTS_FLD_ONE    $fld_one  # - First field to plot
  setenv FFTS_FLD_TWO    $fld_two  # - Second field to plot if present and greater than 0
  setenv FFTS_LMT_FLD    $lmt_fld  # - Which field determines plot limits if present and greater than 0
# setenv FFTS_LOW_YLT    $5        # - Manually chosen lower limit in vertical (Y) axis if present
# setenv FFTS_UPP_YLT    $6        # - Manually chosen upper limit in vertical (Y) axis if present
                                   #
####################################
#
# Next we need to be sure that the appropriate sub-directory
# that is assumed to exist by the batch script will in fact
# exist.
#
###############################################
                                              #
if( ! -d ffts) then                           # upon first run the output
  mkdir ffts                                  # directory will not yet exist
endif                                         # so we create it.
                                              #
###############################################
#
# Next we load the idl module, tell it where to find IDL sources and
# tell it the name of the batch script we want it to execute
#
######################################
                                     #
source /etc/profile.d/modules.csh    #
module load idl                      # load idl onto the compute node
                                     #
set SRCDIR=$PWD/idl/pro/ff_time_series   # idl script source directory
setenv IDL_PATH $IDL_DIR/lib:$SRCDIR #
setenv IDL_STARTUP $PWD/idl/bin/ffts_batch   # idl batch script to execute
                                     #
idl                                  # invoke idl
                                     #
######################################
#
# Now clean things up.
#
#####################################
                                    #
module unload idl                   #
unsetenv IDL_STARTUP                #
                                    #
unsetenv FFTS_DSC_LAB               #
unsetenv FFTS_FLD_ONE               #
#unsetenv FFTS_FLD_TWO              #
#unsetenv FFTS_LTG_FLD              #
#unsetenv FFTS_LOW_YLT              #
#unsetenv FFTS_UPP_YLT              #
                                    #
exit 0                              #
                                    #
#####################################
#
##
###
#####
###### The End ####
