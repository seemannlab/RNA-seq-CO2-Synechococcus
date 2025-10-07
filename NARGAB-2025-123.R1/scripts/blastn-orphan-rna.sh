#####################
# annotation of orphan RNAs
####################

outdir=orphan-rna/annotation/ncbi-blastn-nt_prok/
[ ! -d $outdir ] && mkdir -p $outdir

declare -a orphanrna=( \
 #RNAhub output: infernal-cacofold/RNAhub-cand_6.html -> RF00028: Group I catalytic intron
 "gene-SYNPCC7002_RS02500_gene-SYNPCC7002_RS02510" \
 #RNAhub output: infernal-cacofold/RNAhub-cand_31.html -> partly antisense to CDS
 "gene-SYNPCC7002_RS11335_gene-SYNPCC7002_RS11340" \
 #RNAhub output: infernal-cacofold/RNAhub-cand_36.html
 "gene-SYNPCC7002_RS13890_gene-SYNPCC7002_RS13895" \
 "gene-SYNPCC7002_RS02980_gene-SYNPCC7002_RS02985" \
 "gene-SYNPCC7002_RS06725_gene-SYNPCC7002_RS06730" \
 "gene-SYNPCC7002_RS09150_gene-SYNPCC7002_RS09155" \
 "gene-SYNPCC7002_RS10475_gene-SYNPCC7002_RS10480" \
 "gene-SYNPCC7002_RS13490_gene-SYNPCC7002_RS13495" \
 "gene-SYNPCC7002_RS14650_gene-SYNPCC7002_RS14655"
)

##########

module load blast+
module load anaconda3/1

for i in "${orphanrna[@]}"
do
  fasta="orphan-rna/fasta/$i.fa"
  output="$outdir"`basename ${fasta%.*}`".blastn-nt_prok.txt"
  outfa="$outdir"`basename ${fasta%.*}`".blastn-nt_prok.fa"
  echo $output

  ### NCBI blastn to nr_prokaryotes database of orphan RNAs ###

  blastn -db /home/databases/ncbi/blast/db/nt_prok -query $fasta -task blastn -num_threads 4 -outfmt "6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids stitle sseq" -evalue 1 -out $output

  ### get fasta file of all subjects ###

  fasta2tab $fasta | awk '{print ">"$1; print $2}' > $outfa
  perl -F'\t' -lane '$seq = $F[14] =~ s/-//gr; if($F[8]>$F[9]){$seq=(reverse $seq =~ tr/ACGT/TGCA/r)}; print ">".$F[1].":".$F[8]."-".$F[9]." ".$F[13]."\n".$seq' $output >> $outfa

  ### multiple RNA alignment ###

  input=`basename ${outfa%.*}`
  canona=65
  #for aligner in "Tcoffee" "Tcoffee+hmmer3" "Muscle" "Muscle+hmmer3" "nhmmer"; do
  for aligner in "nhmmer"; do
    bash scripts/multiRNAaligner.sh -f $outfa -a $aligner -c $canona -o "$outdir/denovo-structure/$input-$aligner-$canona"
  done
done
