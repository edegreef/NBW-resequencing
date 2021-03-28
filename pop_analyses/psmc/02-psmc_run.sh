#!/bin/bash

cd /home/degreefe/NBW/psmc

# AR_O-BKW1313-718-505
# IC_O-BKW1330-716-505
# LB_HamLB03-2
# NF_NBW-2017-03
# SS_NBW-2016-09

sample=SS_NBW-2016-09

/home/degreefe/programs/psmc/psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o $sample.N25.t15.r5.psmc $sample.psmcfa
