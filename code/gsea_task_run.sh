#!/bin/bash
#SBATCH --output dsq-gsea_task-%A_%1a-%N.out
#SBATCH --array 0-9
#SBATCH --job-name gsea
#SBATCH --mem-per-cpu=125G -t 10:00:00 -p bigmem,week,day

# DO NOT EDIT LINE BELOW
/gpfs/loomis/apps/avx/software/dSQ/0.96/dSQBatch.py /gpfs/ysm/project/zhao/zy92/NASH/dr_multiomics/code/gsea_task.txt /gpfs/ysm/project/zhao/zy92/NASH/dr_multiomics/code

