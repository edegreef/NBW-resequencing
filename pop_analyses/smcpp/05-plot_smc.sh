#!/bin/bash

# plot smc++ models
# super nice that this program outputs .pdf and .csv nicely labeled so can plot in R easily
# here is example for plotting the 5group (AR,IC,LB,NF,SS) runs

source /home/edegreef/smcpp/bin/activate

smc++ plot -c /scratch/edegreef/SMC_min100kb/smc_AR_plot.pdf /scratch/edegreef/SMC_min100kb/smc_out_timepoints/AR/run_*/model.final.json
smc++ plot -c /scratch/edegreef/SMC_min100kb/smc_IC_plot.pdf /scratch/edegreef/SMC_min100kb/smc_out_timepoints/IC/run_*/model.final.json
smc++ plot -c /scratch/edegreef/SMC_min100kb/smc_LB_plot.pdf /scratch/edegreef/SMC_min100kb/smc_out_timepoints/LB/run_*/model.final.json
smc++ plot -c /scratch/edegreef/SMC_min100kb/smc_NF_plot.pdf /scratch/edegreef/SMC_min100kb/smc_out_timepoints/NF/run_*/model.final.json
smc++ plot -c /scratch/edegreef/SMC_min100kb/smc_SS_plot.pdf /scratch/edegreef/SMC_min100kb/smc_out_timepoints/SS/run_*/model.final.json

