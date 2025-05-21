new="/scratch/project_2007362"

n=$1
curd=$2

source global_env.sh

r3=`sed -n "${n} p" $sfile`
sname=$(basename $r3 | sed "s/.fq.gz//g")

cd ${home_project}/6_tax_profile/bact

set -ex

bfile=${sname}_map.bam
samtools view $bfile | cut -f1,3 | uniq | awk '{print $1"\t"$2}' > ${sname}_map.massign
awk -F'\t' '{n[$2]+=1}END{for(k in n)print k"\t"n[k]}' ${sname}_map.massign > ${sname}_map.rawcounts
samtools depth $bfile | awk -F'\t' '{if($3>0){n[$1]+=1;depth[$1]+=$3};if($3>=3){n3[$1]+=1;depth3[$1]+=$3};if($3>=5){n5[$1]+=1;depth5[$1]+=$3};if($3>=10){n10[$1]+=1;depth10[$1]+=$3}}END{for(contig in n){printf("%s\t%d\t%d\t%d\t%d\t",contig,n[contig],n3[contig],n5[contig],n10[contig]);if(n3[contig]==0){n3[contig]=1};if(n5[contig]==0){n5[contig]=1};if(n10[contig]==0){n10[contig]=1}print depth[contig]/n[contig]"\t"depth3[contig]/n3[contig]"\t"depth5[contig]/n5[contig]"\t"depth10[contig]/n10[contig] }}' > ${sname}_map.basescov
awk -F'\t' '{if(NR==FNR){n[$1]=$2}else{if($1 in n){print $1"\t"n[$1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}}}' ${sname}_map.rawcounts ${sname}_map.basescov > ${sname}.depths
#kraken:taxid|3018277|HumGut_18277_1
awk -F'\t' 'BEGIN{print "ID\tlength\traw_abundance\tbpcov1\tbpcov3\tbpcov5\tbpcov10\tmdepth1\tmdepth3\tmdepth5\tmdepth10"}{if(NR==FNR){n[$1]=$16}else{split($1,a,"|");$1=a[3];gsub("_1$","",$1);if($1 in n)print $1"\t"n[$1]"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}}' $HUMGUTTSV ${sname}.depths > ${sname}.bactstats
rm -f ${sname}_map.massign ${sname}_map.rawcounts ${sname}_map.basescov ${sname}.depths ${sname}_map.bam
cd $curd
