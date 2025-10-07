##########
# check for reads in intergenic loci
# count read pairs with featurecounts 
# (this is the version that is going into the manuscript)

#intergenic regions
a=($(gunzip -c /home/projects/rth/co2capture/subprojects/RNA-seq_1v15-pct-CO2/PCC7002-genome.gff.gz | awk -F "\t" '$3=="region"{print $1}'))
for i in "${a[@]}"; do gunzip -c /home/projects/rth/co2capture/subprojects/RNA-seq_1v15-pct-CO2/PCC7002-genome.gff.gz | awk -F "\t" -v chr=$i -v offset=100 'BEGIN{OFS="\t"; s=0}$3!="exon" && $3!="CDS" && $3!="" && $1==chr{split($9,n,";");sub("ID=","",n[1]);if(s!=0 && e+offset<$4-offset){print $1,$2,"intergenic",e+offset,$4-offset,$6,$7,$8,"ID="gid".."n[1]}; e=$5; gid=n[1]; s=s+1}'; done > orphan-rna/PCC7002-genome-intergenic-offset100-nostrand.gff

#feature counts with both paired-end reads mapping inside the intergenic region
module load subread/2.0.3
for bam in /home/projects/rth/co2capture/subprojects/RNA-seq_1v15-pct-CO2/RNA-browser-Dataset2/analysis/21_sorted/*.bam; do
 bamid=`basename ${bam%_CO2-adaptation.bam*}`
 featureCounts --fracOverlap 0.5 -B -p --countReadPairs -P -D 700 -t intergenic -g ID -a orphan-rna/PCC7002-genome-intergenic-offset100-nostrand.gff -o orphan-rna/featurecounts-readpair/${bamid}-featurecounts-readpair.txt $bam 
done

#combine feature counts from all samples
tmp=$(mktemp -d)
# Extract the gene names etc
input=($(ls orphan-rna/featurecounts-readpair/*-featurecounts-readpair.txt))
tail +2 ${input[0]} | cut -f 1-6 > $tmp/00_annot
# for each file extract only the counts
for i in "${input[@]}"; do
 bsn=$(basename $i .featureCounts)
 echo $bsn > $tmp/$bsn
 tail +3 $i | cut -f 7 >> $tmp/$bsn
done
# collect columns together
paste $tmp/* > orphan-rna/featurecounts-readpair/intergenic-offset100-nostrand-featurecounts-readpairs.txt
rm -rf $tmp
#filter intergenic loci length > 100nts and total raw counts per samples greater than 5000
awk -F "\t" -v cutoff=5000 '$6>100 && ($7+$8+$9+$10>cutoff || $11+$12+$13+$14>cutoff || $15+$16+$17+$18>cutoff || $19+$20+$21+$22>cutoff){print $0 "\t" $7+$8+$9+$10+$11+$12+$13+$14+$15+$16+$17+$18+$19+$20+$21+$22}' orphan-rna/featurecounts-readpair/intergenic-offset100-nostrand-featurecounts-readpairs.txt > orphan-rna/featurecounts-readpair/intergenic-offset100-nostrand-featurecounts-readpairs-l100-c5kpercond.txt
mergeid.pl -a orphan-rna/featurecounts-readpair/intergenic-offset100-nostrand-featurecounts-readpairs-l100-c5kpercond.txt -b orphan-rna/featurecounts/intergenic-offset100-nostrand-featurecounts-l100-c32k-annotated.tsv -i 1,1 -j -o 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,45,46 -t > orphan-rna/featurecounts-readpair/intergenic-offset100-nostrand-featurecounts-readpairs-l100-c5kpercond-annotated.txt
