#!/bin/bash
#
#SBATCH --job-name=host_removal
#SBATCH --output=log/hr_%a.log
#SBATCH --error=err/hr_%a.err
#SBATCH --partition=small
#SBATCH --account=Project_2007362
#SBATCH --time=03:00:00
#SBATCH --array=1-200 # 397
#SBATCH --cpus-per-task=8
#SBATCH --mem=15G
#SBATCH --mail-type=END
#SBATCH --mail-user=sandro.valenzuela@helsinki.fi

#baby_microbiome /scratch/project_2003252/
#new microbiome /scratch/project_2007362/
old="/scratch/project_2003252"
new="/scratch/project_2007362"

n=$SLURM_ARRAY_TASK_ID
r3=`sed -n "${n} p" merged.txt` #qc samples
cd $new/sandro/sep_samples

sname=$(basename $r3 | awk -F"/" '{gsub(".fq.gz","",$NF);print $NF}') #we use $2 because there is qc_ as prefix.

ml minimap2/2.24

minimap2 -t 8 -u both --split-prefix=${sname}.tmp -a $new/software/minimap2DB/human_shortread.mmi $r3 | samtools view -h -f 4 | samtools sort -n | samtools view -bh > ${sname}_nohost.bam
samtools fastq -n ${sname}_nohost.bam > ${sname}.fq
gzip ${sname}.fq
mv ${sname}.fq.gz 1_hostremoval
rm -f ${sname}_nohost.bam

