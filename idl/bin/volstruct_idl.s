# !/bin/bash

##################################################################################
#                                                                                #
# c-shell script: volstruct_idl.s - Timothy J. Dennis 08 July 2016               #
#                             tdennis10@alaska.edu                               #
#                                                                                #
# Developed for use with the family of Heliospheric Physics codes originating    #
# with the work of C.S. Ng and including 'rmct_np', 'fourfields,' and 'coronos'  #
#                                                                                #
# Description: This script is designed to be executed from within a pbs script   #
#              though it can be run from the command line if the user knows      #
#              which arguments to pass to it, and their proper ordering.         #
#              This is, however, not recommended for jobs requiring large        #
#              numbers of plots - particularly when the data is of high          #
#              resolution.                                                       #
#                                                                                #
#              It's purpose is to prepare a node's environment for               #
#              submission of  the batch script 'volstruct_batch' to an           #
#              instantiation of idl on a compute node and allows for the         #
#              generation of...                                                  #
#                                                                                #
# Usage Example:                                                                 #
#                                                                                #
#              ./volstruct_idl.s j, ts, vol, tseries, y                          #
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
       export VSTR_QTY=$1             #
export VSTR_DESC_LABEL=$2             #
      export VSTR_TYPE=$3             #
    export VSTR_ACTION=$4             #
                                      #
if [ -e glb_ext.out ]
then
  export VSTR_GLB_EXT='y'
fi
                                      #
                                      #
#######################################
#
# Now we begin an interactive portion of the script
# for the purposes of specifying optional keywords
# to the base idl script
#


                  prefix="rmct2"
export       VSTR_PREFIX=$prefix

echo -n "pls_min: "
read pls_min

export     VSTR_PLS_MIN=$pls_min
echo -n "slice_range: "
read slice_range
export VSTR_SLICE_RANGE=$slice_range
echo -n "gap: "
read gap
export          VSTR_GAP=$gap
echo -n "window_size: "
read window_size
export VSTR_WINDOW_SIZE=$window_size

if [ "$VSTR_TYPE" == "vol" ]
then
  echo -n "threshold: "
  read threshold
  export VSTR_THRESHOLD=$threshold
elif [ "$VSTR_TYPE" == "isf" ]
then
  echo -n "width: "
  read width
  export      VSTR_WIDTH=$width
  echo -n "n_cntrs: "
  read n_cntrs
  export    VSTR_N_CNTRS=$n_cntrs
  echo -n "inc_cont: "
  read inc_cont
  export   VSTR_INC_CONT=$inc_cont
fi

if [ "$VSTR_ACTION" == "tseries" ]
then
  echo -n "first_step: "
  read first_step
  export VSTR_FIRST_STEP=$first_step
  echo -n "last_step: "
  read last_step
  export VSTR_LAST_STEP=$last_step
  echo -n "yrot: "
  read yrot
  export      VSTR_YROT=$yrot
elif [ "$VSTR_ACTION" == "rotate" ]
then
  echo -n "n_step: "
  read n_step
  export VSTR_N_STEP=$n_step
  echo -n "n_frames: "
  read n_frames
  export VSTR_N_FRAMES=$n_frames
  echo -n "start_angle: "
  read strt_ang
  export VSTR_START_ANGLE=$strt_ang
  echo -n "stop_angle: "
  read stp_ang
  export VSTR_STOP_ANGLE=$stp_ang
elif [ "$VSTR_ACTION" == "scan" ]
then
  echo -n "n_step: "
  read n_step
  export    VSTR_N_STEP=$n_step
  echo -n "n_frames: "
  read n_frames
  export  VSTR_N_FRAMES=$n_frames
  echo -n "yrot: "
  read yrot
  export      VSTR_YROT=$yrot
fi
                                       #
#########################################
#
## Next we need to be sure that the appropriate sub-directory
## that is assumed to exist by the batch script will in fact
## exist.
##
#############################################
#                                           #
OUTDIR=$1'struct'

if [ ! -d $OUTDIR ] 
then
  mkdir -p $OUTDIR/eps
  mkdir    $OUTDIR/pdf
  mkdir    $OUTDIR/jpg
fi

#############################################
##
## Next we load the idl module, tell it where to find IDL sources and
## tell it the name of the batch script we want it to execute
##
##########################################

source /etc/profile.d/modules.csh
module load idl

            SRCDIR=$PWD/idl/volstruct
     echo "SRCDIR = " $SRCDIR
   export IDL_PATH=$IDL_DIR/lib:$SRCDIR
   echo "IDL_PATH = " $IDL_PATH
export IDL_STARTUP=$PWD/volstruct_batch

$IDL_DIR/bin/idl

##########################################
##
## Next initiate animation script
##
#############################################################################################################################
##                                                                                                                          #
##./jsheet_anim.s  $JSHT_DESC_LABEL $JSHT_RES_STR $JSHT_THRSH $JSHT_FIRST_STP $JSHT_LAST_STP $JSHT_FIRST_SLC $JSHT_LAST_SLC #
##                                                                                                                          #
#############################################################################################################################
##
## Now clean things up.
##
#####################################

module unload idl

  export      IDL_STARTUP=""
  export         VSTR_QTY=""
  export  VSTR_DESC_LABEL=""
  export        VSTR_TYPE=""
  export      VSTR_ACTION=""
  export     VSTR_GLB_EXT=""
  export      VSTR_PREFIX=""
  export     VSTR_PLS_MIN=""

  export VSTR_SLICE_RANGE=""
  export   VSTR_THRESHOLD=""
  export       VSTR_WIDTH=""
  export    VSTR_INC_CONT=""
  export  VSTR_FIRST_STEP=""
  export   VSTR_LAST_STEP=""
  export      VSTR_N_STEP=""
  export    VSTR_N_FRAMES=""
  export VSTR_START_ANGLE=""
  export  VSTR_STOP_ANGLE=""
  export         VSTR_GAP=""
  export VSTR_WINDOW_SIZE=""
  export        VSTR_YROT=""

exit 0

####################################
##
###
#### The End ####
