#! /bin/bash
###
###Script used to run the pipeline
###
#SBATCH -t 30:00:00
#SBATCH --partition=normal_q
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tnchau@vt.edu
#SBATCH --output=spmarker_Rice_111022.out


echo "Starting"
date

time python SPmarker.py \
            -d 111022_Rice_work/ -o 111022_Rice_out/ \
            -mtx 111022_Rice_expression.csv \
            -meta 111022_Rice_cellType.csv


echo "Finished"
date

exit;
