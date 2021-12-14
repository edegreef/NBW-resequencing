#!/bin/bash


source /home/edegreef/smcpp/bin/activate

smc++ plot -c /scratch/edegreef/smcpp/plot/smc_ICSS_split_plot.pdf /scratch/edegreef/smcpp/split_estimate/ICSS/combined_run*/model.final.json
