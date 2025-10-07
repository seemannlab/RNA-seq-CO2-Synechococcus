mydir=featurecounts-no-genes
[ ! -d $mydir ] && mkdir $mydir

#RefSeq riboswitches
gunzip -c /home/projects/rth/co2capture/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/raw-data/PCC7002-genome.gff.gz | awk -F "\t" '$3=="riboswitch"{sub("Dbxref","Name"); print $0}' > $mydir/refseq-no-genes.gff

#RefSeq pseudogenes
gunzip -c /home/projects/rth/co2capture/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/raw-data/PCC7002-genome.gff.gz | awk -F "\t" '$3=="pseudogene"' >> $mydir/refseq-no-genes.gff

#isaR1
awk 'BEGIN{OFS="\t"}NR==3{print $1,"cmsearch","sRNA",$9,$8,".",$10,".","ID="$3}' srna/isar1/isar1.txt >> $mydir/refseq-no-genes.gff

#yfr1
awk 'BEGIN{OFS="\t"}NR==3{print $1,"cmsearch","sRNA",$9,$8,".",$10,".","ID="$3}' srna/yfr1/yfr1.txt >> $mydir/refseq-no-genes.gff

#CRSs
awk -F "\t" 'BEGIN{OFS="\t"}$6=="CRS"{split($1,a,";"); sub(".fna.motif", "", a[1]); print $2,$6,$6,$3,$4,".",$5,".","ID="a[1]}' /home/projects/rth/co2capture/subprojects/RNA-seq_1v15-pct-CO2/RNA-seq-CO2-Synechococcus/data/C_annotation.tsv >> $mydir/refseq-no-genes.gff

#feature counts with both paired-end reads mapping inside the intergenic region
module load subread/2.0.3
for bam in /home/projects/rth/co2capture/subprojects/RNA-seq_1v15-pct-CO2/RNA-browser-Dataset2/analysis/21_sorted/*.bam; do
 bamid=`basename ${bam%_CO2-adaptation.bam*}`
 featureCounts -s 2 -p --fracOverlap 0.5 -B -P -D 700 -O -t riboswitch,pseudogene,sRNA,CRS -g ID -a $mydir/refseq-no-genes.gff -o $mydir/${bamid}-featurecounts.txt $bam 
done

#combine feature counts from all samples
tmp=$(mktemp -d)
# Extract the gene names etc
input=($(ls $mydir/*-featurecounts.txt))
tail +2 ${input[0]} | cut -f 1-6 > $tmp/00_annot
# for each file extract only the counts
for i in "${input[@]}"; do
 bsn=$(basename $i .featureCounts)
 echo $bsn > $tmp/$bsn
 tail +3 $i | cut -f 7 >> $tmp/$bsn
done
# collect columns together
paste $tmp/* > $mydir/refseq-no-genes-featurecounts-readpairs.txt
rm -rf $tmp
