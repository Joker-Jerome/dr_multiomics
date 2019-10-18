#!/bin/bash
#SBATCH --output dsq-gsea_task-%A_%1a-%N.out
#SBATCH --array 0-9
#SBATCH --job-name gsea
#SBATCH --mem-per-cpu=32G -t 10:00:00 -p general,bigmem,pi_zhao

# DO NOT EDIT LINE BELOW
/ysm-gpfs/apps/software/dSQ/0.96/dSQBatch.py /gpfs/ysm/project/zhao/zy92/NASH/dr_multiomics/code/gsea_task.txt /gpfs/ysm/project/zhao/zy92/NASH/dr_multiomics/code

