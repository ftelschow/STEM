#!/bin/bash
########################################################################
####  script for automatically plotting all figures of simulations  ####
########################################################################
# create a file containing the simulations for the T plots
ls simulations/*_T_*_prethresh[!-]_transform* > Tsimulations.txt
ls simulations/FieldTYPE_Z_Msim1000_* > Histsimulations.txt
ls simulations/FieldTYPE_Z_N20Msim10000_* > Zsimulations.txt

# plot the T-field simulations using matlab
matlab -nodesktop -nosplash -r "run scripts/plotscript.m"
rm Tsimulations.txt Histsimulations.txt Zsimulations.txt

