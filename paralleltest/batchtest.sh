#!/bin/csh
#PBS -l nodes=5:ppn=8
#PBS -q workq
#PBS -l walltime = 00:00:01
#PBS -N batchtest
cd /mnt/scratch-lustre/njones/Summer2014/paralleltest
module load python/2.7.6
python motherscript.py