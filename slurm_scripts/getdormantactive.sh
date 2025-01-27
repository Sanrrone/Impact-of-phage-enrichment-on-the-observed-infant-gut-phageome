#!/bin/bash
#
new="/scratch/project_2007362"

source global_env.sh

export PATH=$new/software/mambaforge/bin:$PATH
eval "$(conda shell.bash hook)"
conda activate $PROP_ENV

r3=$1
sname_mag=$2
sname=$(basename $r3 | sed "s/.fq.gz//g")
bfolder="${sname_mag}_bams"
c=3
assembly=${home_project}/2_assembly/${sname_mag}.fna

ml minimap2
##############################################################################################
cd ${home_project}/3_phage_annotation
set -ex
bfile=${sname}.bam

if [ ! -f "$bfolder/$bfile" ];then
minimap2 -t $c --eqx --split-prefix=${sname}_propagate_tmp -a $assembly $r3 | samtools view -h -F 2308 | samtools sort | samtools view -bh > ${home_project}/2_assembly/$bfolder/$bfile
fi
echo -e "scaffold\tfragment\tstart\tstop" > ${sname}_procoords.tsv
awk -F"\t" -v s=${sname_mag} '{if($1==s)print $3"\t"$2"\t"$4"\t"$5}' ~/phage_approach/supp_files/prophages_coordinates.tsv >> ${sname}_procoords.tsv
rm -rf ${sname}_propagate
Propagate -f $assembly -b ${home_project}/2_assembly/$bfolder/$bfile -v ${sname}_procoords.tsv -o ${sname}_propagate -c 1.5 -t $c --clean
mv ${sname}_propagate/${sname}_propagate.tsv .
rm -rf ${sname}_propagate ${sname}_procoords.tsv

