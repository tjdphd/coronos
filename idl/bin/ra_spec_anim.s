#!/bin/bash

 desc_label=$1
    res_str=$2
        qty=$3
first_layer=$4
 last_layer=$5
 first_step=$6
  last_step=$7
  anim_mode=$8

zero_str='000000000000000000000000000000'
wild_str='??????????????????????????????'
     sep='_'
     hph='-'
 cur_dir=`pwd` # directory originating script invocation

cd ra_spec

if [[ "$anim_mode" == "vsz" ]]
then
  if [ $first_step -eq $last_step ]
  then
#
    len_first_step=`expr length $first_step`
#
    if [ $len_first_step -lt 3 ]
    then
          len_dif=`expr 3 - $len_first_step`
              pos=`expr 30 - $len_dif`
    trnc_zero_str=${zero_str:$pos}
       first_step=$trnc_zero_str$first_step
    fi
#
    len_last_layer=`expr length $last_layer`
#
    if [ $len_last_layer -lt 3 ]
    then
      trnc_wild_str='???'
    else
                pos=`expr 30 - $len_last_layer`
      trnc_wild_str=${wild_str:$pos}
    fi
#
    spec_eps=`ls -1 $qty/eps/$qty$sep'ra_spec'$sep$res_str$sep'zl_'*$sep'proc'$hph$trnc_wild_str$sep'layer'$hph$trnc_wild_str$sep'srun'$sep$first_step'.eps'`
    spec_gif_anim=$qty'_ra_spec_step'$sep$first_step$sep'layers'$sep$first_layer$hph$last_layer$sep$res_str$hph$desc_label'-anim.gif'
    spec_mpf_anim=$qty'_ra_spec_step'$sep$first_step$sep'layers'$sep$first_layer$hph$last_layer$sep$res_str$hph$desc_label'-anim.mp4'
  else
    echo "ra_spec_anim: error - cannot simultaneously animate for varying z and t. Exiting..."
    exit 1
  fi
else
  if [[ "$anim_mode" == "vst" ]]
  then
    echo "why hello there!"
    if [ $first_layer -eq $last_layer ]
    then
      len_first_layer=`expr length $first_layer`
#
      if [ $len_first_layer -lt 3 ]
      then
             len_dif=`expr 3 - $len_first_layer`
                 pos=`expr 30 - $len_dif`
       trnc_zero_str=${zero_str:$pos}
         first_layer=$trnc_zero_str$first_layer
      fi
#
      len_last_step=`expr length $last_step`
#
      if [ $len_last_step -lt 3 ]
      then
        trnc_wild_str='???'
      else
                  pos=`expr 30 - $len_last_step`
        trnc_wild_str=${wild_str:$pos}
      fi
#
      spec_eps=`ls -1 $qty/eps/$qty$sep'ra_spec'$sep$res_str$sep'zl_'*$sep'proc'$hph*$sep'layer'$hph$first_layer$sep'srun'$sep$trnc_wild_str'.eps'`
      spec_gif_anim=$qty'_ra_spec$hph'layer'$sep$first_layer'$sep'steps'$sep$first_step$hph$last_step$sep$res_str$hph$desc_label'-anim.gif'
      spec_mpf_anim=$qty'_ra_spec$hph'layer'$sep$first_layer'$sep'steps'$sep$first_step$hph$last_step$sep$res_str$hph$desc_label'-anim.mp4'
    fi
  else
    echo "ra_spec_anim: error - cannot simultaneously animate for varying z and t. Exiting..."
    exit 1
  fi
fi

if  ! [   -e $qty/png ]
then
  mkdir $qty/png
fi

trunc=$qty'/eps/'
echo "trunc = " $trunc
for plot in $spec_eps; do
        plot=${plot#$trunc}
    len_plot=${#plot}                                    # Get the file's basename and create a gif version
           n=`expr $len_plot - 4`                        # with convert
   base_name=${plot:0:$n}
   echo "creating = " $base_name'.png'
   convert -density 300x300 $qty/eps/$base_name.eps -background white -flatten -resize 1024x1024 $qty/png/$base_name.png 
done

spec_png=`ls -1 $qty/png/$qty$sep'ra_spec'$sep$res_str$sep'zl_'*$sep'proc'$hph$trnc_wild_str$sep'layer'$hph$trnc_wild_str$sep'srun'$sep$first_step'.png'`
convert -delay 20 -dispose previous -loop 0 $spec_png $spec_gif_anim

$HOME/local/bin/ffmpeg -framerate 24 -pattern_type glob -i $qty/png/"*.png" -pix_fmt yuv420p $spec_mpf_anim

exit 0
