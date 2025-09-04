#!/bin/bash
#
#SBATCH --job-name=merge
#SBATCH --output=log/merge_%a.txt
#SBATCH --error=err/merge_%a.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --time=00:30:00
#SBATCH --array=401-481 #698 699
#SBATCH --cpus-per-task=8
#SBATCH --mem=15GB
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

#baby_microbiome /scratch/project_2003252/
#new microbiome /scratch/project_2007362/
old="/scratch/project_2003252"
new="/scratch/project_2007362"

n=$SLURM_ARRAY_TASK_ID
pair=`sed -n "${n} p" pairs.txt`
cd $new/sandro/sep_samples

r1=$(echo "$pair" | awk '{print $1}')
r2=$(echo "$pair" | awk '{print $2}')
sname=$(basename $r1 | awk -F"_" '{print $1}')

fastp -i $r1 -I $r2 --merge --include_unmerged --merged_out "${sname}.fq.gz" -h ${sname}.html -j ${sname}.json -l 50 -q 20 -n 0 -w 8 -f 2 -t 2 -F 2 -T 2 --dedup --dup_calc_accuracy 3 --overlap_len_require 20 --fix_mgi_id --detect_adapter_for_pe --correction --overlap_diff_limit 3 

mv -f ${sname}.fq.gz 0_merge/.
mv -f ${sname}.html ${sname}.json 0_merge/logs
