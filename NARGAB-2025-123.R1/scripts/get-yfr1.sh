# build CM of sRNA Yfr1 and search PCC 7002
#
#get alignment from Fig 2 in https://doi.org/10.1186/1471-2164-8-375 --> yfr1.fasta

outdir=srna/yfr1

function fa2sto () {
        local in=$1
        local out=$2

        local tempf=`mktemp`
        fasta2tab $in | cut -f1,2 > $tempf

        local noseq=$(head -n -1 $tempf | wc -l)

        echo "# STOCKHOLM 1.0" > $out
        echo "#=GF AU Realignment" >> $out
        echo -e "#=GF SQ "$noseq"\n" >> $out
        awk '{if($1=="structure"){printf "%-50s %s\n", "#=GC SS_cons", $2}else{printf "%-50s %s\n", $1, $2}}' $tempf >> $out
        echo "//" >> $out
}

#get sec struct with PETfold
module load PETfold/2.2
PETfold -f $outdir/yfr1.fasta -g 0.7 --war | tail -2 | sed 's/-/./g' >> $outdir/yfr1.fasta

#fa2sto
fa2sto $outdir/yfr1.fasta $outdir/yfr1.sto

#get CM
module load infernal/1.1.4
cmbuild -F $outdir/yfr1.cm $outdir/yfr1.sto

#calibrate CM
cmcalibrate $outdir/yfr1.cm

#cmsearch genome
cmsearch --tblout $outdir/yfr1.txt --cpu 4 $outdir/yfr1.cm /home/projects/rth/co2capture/subprojects/RNA-seq_1v15-pct-CO2/Rfam-candidates/genome.fna.gz

#get sequence of best hit at position NC_010475.1:1694670-1694606 between genes SYNPCC7002_RS08085 (trxA, thioredoxin) and SYNPCC7002_RS08090 (glycosyltransferase family 4 protein)
subseq.pl 1694606 1694670 /home/projects/rth/co2capture/subprojects/RNA-seq_1v15-pct-CO2/Rfam-candidates/genome.fna.gz | perl -ane '$o=(reverse $F[0] =~ tr/ATGC/TACG/r); print ">Syn_PCC7002\n".$o."\n"' > $outdir/yfr1-pcc7002.fa

#realign sequences including PCC 7002 to the built covariance model
cat $outdir/yfr1-pcc7002.fa $outdir/yfr1.fasta | head -n -2 | awk '{if(/^>/){print $0}else{gsub("-",""); print $0}}' > $outdir/yfr1.incl-pcc7002.ungapped.fa
cmalign $outdir/yfr1.cm $outdir/yfr1.incl-pcc7002.ungapped.fa > $outdir/yfr1.incl-pcc7002.cm
