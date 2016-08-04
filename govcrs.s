#!/bin/bash

 job_1=`qsub rmct2-64_16-1.pbs`
 job_2=`qsub -W depend=afterok:$job_1 rmct2-64_16-2.pbs`
 job_3=`qsub -W depend=afterok:$job_2 rmct2-64_16-3.pbs`
 job_4=`qsub -W depend=afterok:$job_3 rmct2-64_16-4.pbs`
 job_5=`qsub -W depend=afterok:$job_4 rmct2-64_16-5.pbs`
 job_6=`qsub -W depend=afterok:$job_5 rmct2-64_16-6.pbs`
 job_7=`qsub -W depend=afterok:$job_6 rmct2-64_16-7.pbs`
 job_8=`qsub -W depend=afterok:$job_7 rmct2-64_16-8.pbs`
 job_9=`qsub -W depend=afterok:$job_8 rmct2-64_16-9.pbs`

job_10=`qsub -W depend=afterok:$job_9  rmct2-64_16-10.pbs`
job_11=`qsub -W depend=afterok:$job_10 rmct2-64_16-11.pbs`
job_12=`qsub -W depend=afterok:$job_11 rmct2-64_16-12.pbs`
job_13=`qsub -W depend=afterok:$job_12 rmct2-64_16-13.pbs`
job_14=`qsub -W depend=afterok:$job_13 rmct2-64_16-14.pbs`
job_15=`qsub -W depend=afterok:$job_14 rmct2-64_16-15.pbs`
job_16=`qsub -W depend=afterok:$job_15 rmct2-64_16-16.pbs`
job_17=`qsub -W depend=afterok:$job_16 rmct2-64_16-17.pbs`
job_18=`qsub -W depend=afterok:$job_17 rmct2-64_16-18.pbs`
job_19=`qsub -W depend=afterok:$job_18 rmct2-64_16-19.pbs`

job_20=`qsub -W depend=afterok:$job_19 rmct2-64_16-20.pbs`
job_21=`qsub -W depend=afterok:$job_20 rmct2-64_16-21.pbs`
job_22=`qsub -W depend=afterok:$job_21 rmct2-64_16-22.pbs`
job_23=`qsub -W depend=afterok:$job_22 rmct2-64_16-23.pbs`
job_24=`qsub -W depend=afterok:$job_23 rmct2-64_16-24.pbs`
job_25=`qsub -W depend=afterok:$job_24 rmct2-64_16-25.pbs`
job_26=`qsub -W depend=afterok:$job_25 rmct2-64_16-26.pbs`
job_27=`qsub -W depend=afterok:$job_26 rmct2-64_16-27.pbs`
job_28=`qsub -W depend=afterok:$job_27 rmct2-64_16-28.pbs`
job_29=`qsub -W depend=afterok:$job_28 rmct2-64_16-29.pbs`

job_30=`qsub -W depend=afterok:$job_29 rmct2-64_16-30.pbs`
job_31=`qsub -W depend=afterok:$job_30 rmct2-64_16-31.pbs`
job_32=`qsub -W depend=afterok:$job_31 rmct2-64_16-32.pbs`
job_33=`qsub -W depend=afterok:$job_32 rmct2-64_16-33.pbs`
job_34=`qsub -W depend=afterok:$job_33 rmct2-64_16-34.pbs`
job_35=`qsub -W depend=afterok:$job_34 rmct2-64_16-35.pbs`
job_36=`qsub -W depend=afterok:$job_35 rmct2-64_16-36.pbs`
job_37=`qsub -W depend=afterok:$job_36 rmct2-64_16-37.pbs`
job_38=`qsub -W depend=afterok:$job_37 rmct2-64_16-38.pbs`
job_39=`qsub -W depend=afterok:$job_38 rmct2-64_16-39.pbs`

job_40=`qsub -W depend=afterok:$job_39 rmct2-64_16-40.pbs`
job_41=`qsub -W depend=afterok:$job_40 rmct2-64_16-41.pbs`
job_42=`qsub -W depend=afterok:$job_41 rmct2-64_16-42.pbs`
job_43=`qsub -W depend=afterok:$job_42 rmct2-64_16-43.pbs`
job_44=`qsub -W depend=afterok:$job_43 rmct2-64_16-44.pbs`
job_45=`qsub -W depend=afterok:$job_44 rmct2-64_16-45.pbs`
job_46=`qsub -W depend=afterok:$job_45 rmct2-64_16-46.pbs`
job_47=`qsub -W depend=afterok:$job_46 rmct2-64_16-47.pbs`
job_48=`qsub -W depend=afterok:$job_47 rmct2-64_16-48.pbs`
job_49=`qsub -W depend=afterok:$job_48 rmct2-64_16-49.pbs`

job_50=`qsub -W depend=afterok:$job_49 rmct2-64_16-50.pbs`
job_51=`qsub -W depend=afterok:$job_50 rmct2-64_16-51.pbs`
job_52=`qsub -W depend=afterok:$job_51 rmct2-64_16-52.pbs`
job_53=`qsub -W depend=afterok:$job_52 rmct2-64_16-53.pbs`
job_54=`qsub -W depend=afterok:$job_53 rmct2-64_16-54.pbs`
job_55=`qsub -W depend=afterok:$job_54 rmct2-64_16-55.pbs`
job_56=`qsub -W depend=afterok:$job_55 rmct2-64_16-56.pbs`
job_57=`qsub -W depend=afterok:$job_56 rmct2-64_16-57.pbs`
job_58=`qsub -W depend=afterok:$job_57 rmct2-64_16-58.pbs`
job_59=`qsub -W depend=afterok:$job_58 rmct2-64_16-59.pbs`

job_60=`qsub -W depend=afterok:$job_59 rmct2-64_16-60.pbs`
job_61=`qsub -W depend=afterok:$job_60 rmct2-64_16-61.pbs`
job_62=`qsub -W depend=afterok:$job_61 rmct2-64_16-62.pbs`
job_63=`qsub -W depend=afterok:$job_62 rmct2-64_16-63.pbs`
job_64=`qsub -W depend=afterok:$job_63 rmct2-64_16-64.pbs`
job_65=`qsub -W depend=afterok:$job_64 rmct2-64_16-65.pbs`
job_66=`qsub -W depend=afterok:$job_65 rmct2-64_16-66.pbs`
job_67=`qsub -W depend=afterok:$job_66 rmct2-64_16-67.pbs`
job_68=`qsub -W depend=afterok:$job_67 rmct2-64_16-68.pbs`
job_69=`qsub -W depend=afterok:$job_68 rmct2-64_16-69.pbs`

job_70=`qsub -W depend=afterok:$job_69 rmct2-64_16-70.pbs`
job_71=`qsub -W depend=afterok:$job_70 rmct2-64_16-71.pbs`
job_72=`qsub -W depend=afterok:$job_71 rmct2-64_16-72.pbs`
job_73=`qsub -W depend=afterok:$job_72 rmct2-64_16-73.pbs`
job_74=`qsub -W depend=afterok:$job_73 rmct2-64_16-74.pbs`
job_75=`qsub -W depend=afterok:$job_74 rmct2-64_16-75.pbs`
job_76=`qsub -W depend=afterok:$job_75 rmct2-64_16-76.pbs`
job_77=`qsub -W depend=afterok:$job_76 rmct2-64_16-77.pbs`
job_78=`qsub -W depend=afterok:$job_77 rmct2-64_16-78.pbs`
job_79=`qsub -W depend=afterok:$job_78 rmct2-64_16-79.pbs`

job_80=`qsub -W depend=afterok:$job_79 rmct2-64_16-80.pbs`
job_81=`qsub -W depend=afterok:$job_80 rmct2-64_16-81.pbs`
job_82=`qsub -W depend=afterok:$job_81 rmct2-64_16-82.pbs`
job_83=`qsub -W depend=afterok:$job_82 rmct2-64_16-83.pbs`
job_84=`qsub -W depend=afterok:$job_83 rmct2-64_16-84.pbs`
job_85=`qsub -W depend=afterok:$job_84 rmct2-64_16-85.pbs`
job_86=`qsub -W depend=afterok:$job_85 rmct2-64_16-86.pbs`
job_87=`qsub -W depend=afterok:$job_86 rmct2-64_16-87.pbs`
job_88=`qsub -W depend=afterok:$job_87 rmct2-64_16-88.pbs`
job_89=`qsub -W depend=afterok:$job_88 rmct2-64_16-89.pbs`

job_90=`qsub -W depend=afterok:$job_89 rmct2-64_16-90.pbs`
job_91=`qsub -W depend=afterok:$job_90 rmct2-64_16-91.pbs`
job_92=`qsub -W depend=afterok:$job_91 rmct2-64_16-92.pbs`
job_93=`qsub -W depend=afterok:$job_92 rmct2-64_16-93.pbs`
job_94=`qsub -W depend=afterok:$job_93 rmct2-64_16-94.pbs`
job_95=`qsub -W depend=afterok:$job_94 rmct2-64_16-95.pbs`
job_96=`qsub -W depend=afterok:$job_95 rmct2-64_16-96.pbs`
job_97=`qsub -W depend=afterok:$job_96 rmct2-64_16-97.pbs`
job_98=`qsub -W depend=afterok:$job_97 rmct2-64_16-98.pbs`
job_99=`qsub -W depend=afterok:$job_98 rmct2-64_16-99.pbs`

job_100=`qsub -W depend=afterok:$job_99  rmct2-64_16-100.pbs`

exit 0
