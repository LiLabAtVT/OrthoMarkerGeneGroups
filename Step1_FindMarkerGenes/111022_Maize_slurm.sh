#! /bin/bash
###
###Script used to run the pipeline
###
#SBATCH -t 4:00:00
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH -p p100_dev_q
#SBATCH -A introtogds
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tnchau@vt.edu
#SBATCH --output=111022_spmarker_Maize.out


echo "Starting"
date

time python SPmarker.py \
            -d 111022_Maize_work/ -o 111022_Maize_out/ \
            -mtx 111022_Maize_expression.csv \
            -meta 111022_Maize_cellType.csv


echo "Finished"
date

exit;
