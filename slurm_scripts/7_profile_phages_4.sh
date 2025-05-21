#!/bin/bash
#
#SBATCH --job-name=s7_propagate
#SBATCH --output=log/s7_s4_%a.log
#SBATCH --error=err/s7_s4_%a.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --array=48-55
#SBATCH --time=04:15:00 # set to 3 hours when all samples are run # 8h
#SBATCH --cpus-per-task=8
#SBATCH --mem=130G #150
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

new="/scratch/project_2007362"
n=$SLURM_ARRAY_TASK_ID
source global_env.sh
 
sname_mag=`sed -n "${n} p" $sepsafile`
sep=$(echo $sname_mag | awk -F"_" '{print $1}')
sa=$(echo $sname_mag | sed "s/${sep}_//g" )

set -ex 
awk -F"\t" -v h=$sep -v sa=$sa -v p="/scratch/project_2007362/sandro/HeP_samples/1_hostremoval/" '{if($1==h && $4==sa){print p$2".fq.gz"}}' $sampletable |
	xargs -P 8 -n 1 -I {} bash ./getdormantactive.sh {} ${sname_mag}

wait
#P 10 parallel processes
#n 1 command per independent seq 1,2,3,4...
#I {} the value of seq command, 1,2,3,4...

