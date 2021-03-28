#!/bin/bash

# plot psmc results (-R prints intermediate outputs which help for plotting in R, -u mutation rate -g generation time)
# super quick to run so just doing on command line

/home/degreefe/programs/psmc/utils/psmc_plot.pl -R -u 1.53e-08 -g 17.8 -p AR_1.53e-8_g17.8_plot AR_O-BKW1313-718-505.N25.t15.r5.psmc

/home/degreefe/programs/psmc/utils/psmc_plot.pl -R -u 1.53e-08 -g 17.8 -p IC_1.53e-8_g17.8_plot IC_O-BKW1330-716-505.N25.t15.r5.psmc

/home/degreefe/programs/psmc/utils/psmc_plot.pl -R -u 1.53e-08 -g 17.8 -p LB_1.53e-8_g17.8_plot LB_HamLB03-2.N25.t15.r5.psmc

/home/degreefe/programs/psmc/utils/psmc_plot.pl -R -u 1.53e-08 -g 17.8 -p NF_1.53e-8_g17.8_plot NF_NBW-2017-03.N25.t15.r5.psmc

/home/degreefe/programs/psmc/utils/psmc_plot.pl -R -u 1.53e-08 -g 17.8 -p SS_1.53e-8_g17.8_plot SS_NBW-2016-09.N25.t15.r5.psmc
