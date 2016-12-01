#!/bin/bash -l
#PBS -l walltime=96:00:00,nodes=1:ppn=13,mem=2580mb
#PBS -m abe
#PBS -M perre035@umn.edu

cd ~/kitaev
module load matlab
matlab -nodisplay -r "maxNumCompThreads(13)" < do_some_mcmc.m