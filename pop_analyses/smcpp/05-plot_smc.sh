#!/bin/bash

# plot smc++ models
# super nice that this program outputs .pdf and .csv nicely labeled so can plot in R easily
#salloc --ntasks=1 --cpus-per-task=1 --mem-per-cpu=10GB --time=2:00:00 --account=def-coling

source /home/edegreef/smcpp/bin/activate

smc++ plot -c /scratch/edegreef/smcpp/plot/smc_ALLAR_plot.pdf /scratch/edegreef/smcpp/smc_out_timepoints/ALL.AR-NBW-12-710-509/run_*/model.final.json
smc++ plot -c /scratch/edegreef/smcpp/plot/smc_ALLIC_plot.pdf /scratch/edegreef/smcpp/smc_out_timepoints/ALL.IC-Ha14-01-711-503/run_*/model.final.json
smc++ plot -c /scratch/edegreef/smcpp/plot/smc_ALLLB_plot.pdf /scratch/edegreef/smcpp/smc_out_timepoints/ALL.LB-HamLB03-1/run_*/model.final.json
smc++ plot -c /scratch/edegreef/smcpp/plot/smc_ALLNF_plot.pdf /scratch/edegreef/smcpp/smc_out_timepoints/ALL.NF-Hyam-2016-06-pooled/run_*/model.final.json
smc++ plot -c /scratch/edegreef/smcpp/plot/smc_ALLSS_plot.pdf /scratch/edegreef/smcpp/smc_out_timepoints/ALL.SS-NBW-2016-09/run_*/model.final.json

smc++ plot -c /scratch/edegreef/smcpp/plot/smc_ARICAR_plot.pdf /scratch/edegreef/smcpp/smc_out_timepoints/ARIC.AR/run_*/model.final.json
smc++ plot -c /scratch/edegreef/smcpp/plot/smc_ARICIC_plot.pdf /scratch/edegreef/smcpp/smc_out_timepoints/ARIC.IC/run_*/model.final.json
smc++ plot -c /scratch/edegreef/smcpp/plot/smc_LBNFLB_plot.pdf /scratch/edegreef/smcpp/smc_out_timepoints/LBNF.LB/run_*/model.final.json
smc++ plot -c /scratch/edegreef/smcpp/plot/smc_LBNFNF_plot.pdf /scratch/edegreef/smcpp/smc_out_timepoints/LBNF.NF/run_*/model.final.json

smc++ plot -c /scratch/edegreef/smcpp/plot/smc_SS_plot.pdf /scratch/edegreef/smcpp/smc_out_timepoints/SS/run_*/model.final.json
smc++ plot -c /scratch/edegreef/smcpp/plot/smc_LB_plot.pdf /scratch/edegreef/smcpp/smc_out_timepoints/LB/run_*/model.final.json
smc++ plot -c /scratch/edegreef/smcpp/plot/smc_NF_plot.pdf /scratch/edegreef/smcpp/smc_out_timepoints/NF/run_*/model.final.json
smc++ plot -c /scratch/edegreef/smcpp/plot/smc_AR_plot.pdf /scratch/edegreef/smcpp/smc_out_timepoints/AR/run_*/model.final.json
smc++ plot -c /scratch/edegreef/smcpp/plot/smc_IC_plot.pdf /scratch/edegreef/smcpp/smc_out_timepoints/IC/run_*/model.final.json
