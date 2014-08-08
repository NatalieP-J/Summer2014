#!/bin/csh
#PBS -l nodes=1:ppn=8
#PBS -q workq
#PBS -l walltime=00:01:00
#PBS -N batchtest
cd $PBS_O_WORKDIR
module load python/2.7.6
python motherscript.py
