# build CM of sRNA IsaR1 and search PCC 7002
#
#get alignment from Fig S1B in https://doi.org/10.1016/j.cub.2017.04.010 --> isar1.fasta

outdir=srna/isar1

#get sec struct with PETfold
export PETFOLDBIN=~/projects-git/petfold/bin
~/projects-git/petfold/bin/PETfold -f $outdir/isar1.fasta

#get CM
module load infernal/1.1.4
cmbuild -F $outdir/isar1.cm $outdir/isar1.sto

#calibrate CM
cmcalibrate $outdir/isar1.cm

#cmsearch genome
cmsearch --tblout $outdir/isar1.txt --cpu 4 $outdir/isar1.cm /home/projects/rth/co2capture/subprojects/RNA-seq_1v15-pct-CO2/Rfam-candidates/genome.fna.gz

#get sequence of best hit at position 1342543 to 1342601 between genes SYNPCC7002_RS06465 (chlorophyll a/b binding light-harvesting protein) and flavodoxin FldA
subseq.pl 1342543 1342601 /home/projects/rth/co2capture/subprojects/RNA-seq_1v15-pct-CO2/Rfam-candidates/genome.fna.gz | perl -ane '$o=(reverse $F[0] =~ tr/ATGC/TACG/r); print $o."\n"' > $outdir/isar1-pcc7002.seq

#redo the same for RNA cons structure from RNAalifold
module load ViennaRNA2/2.6.3
fa2aln.pl $outdir/isar1.fasta | RNAalifold
