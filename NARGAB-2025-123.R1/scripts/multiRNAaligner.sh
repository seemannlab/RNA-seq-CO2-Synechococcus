#!/bin/bash
#set -euxo pipefail
# Copyright (C) 2025 Stefan E Seemann <seemann@rth.dk>

#▗▖  ▗▖▗▖ ▗▖▗▖ ▗▄▄▄▖▗▄▄▄▖    ▗▄▄▖ ▗▖  ▗▖ ▗▄▖      ▗▄▖ ▗▖   ▗▄▄▄▖ ▗▄▄▖▗▖  ▗▖▗▄▄▄▖▗▄▄▖
#▐▛▚▞▜▌▐▌ ▐▌▐▌   █    █      ▐▌ ▐▌▐▛▚▖▐▌▐▌ ▐▌    ▐▌ ▐▌▐▌     █  ▐▌   ▐▛▚▖▐▌▐▌   ▐▌ ▐▌
#▐▌  ▐▌▐▌ ▐▌▐▌   █    █      ▐▛▀▚▖▐▌ ▝▜▌▐▛▀▜▌    ▐▛▀▜▌▐▌     █  ▐▌▝▜▌▐▌ ▝▜▌▐▛▀▀▘▐▛▀▚▖
#▐▌  ▐▌▝▚▄▞▘▐▙▄▄▖█  ▗▄█▄▖    ▐▌ ▐▌▐▌  ▐▌▐▌ ▐▌    ▐▌ ▐▌▐▙▄▄▖▗▄█▄▖▝▚▄▞▘▐▌  ▐▌▐▙▄▄▖▐▌ ▐▌


source /etc/profile.d/modules.sh
module use /home/local/modulefiles
if [[ -f /home/local/opt/miniconda3/etc/profile.d/conda.sh ]]; then
        source /home/local/opt/miniconda3/etc/profile.d/conda.sh
        conda activate python3
fi

# Function to display the help message
usage() {
  echo "Usage: $(basename "$0") [OPTIONS]"
  echo "Multiple RNA alignment based on"
  echo "(1) sequence alignment (using T-coffee, Muscle, T-coffee+hmmer3, Muscle+hmmer3, and nhmmer)"
  echo "(2) consensus structure prediction (using PETfold)"
  echo "(3) structural realignment (using Infernal)"
  echo "(4) filtering of sequences with high number of non-canonical basepairs."
  echo
  echo "Options:"
  echo "  -h              Display this help message"
  echo "  -f <file>       Input FASTA file (required)"
  echo "  -a <aligner>    Aligner (optional, options: Tcoffee, Tcoffee+hmmer3, Muscle, Muscle+hmmer3, default: Muscle+hmmer3)"
  echo "  -c <perc>       Remove sequences with greater than percentage non-canonical basepairs in the consensus structure (optional, default: 65)"
  echo "  -o <directory>  Directory where output should be written (optional, default: denovo-structure-<file>-<aligner>-<perc>)"
  echo
  echo "Example:"
  echo "  $(basename "$0") -f sequences.fasta -a Muscle+hmmer3 -c 65 -o denovo-structure"
  exit 0
}

# Default values for options
aligner="Muscle+hmmer3"
canona=65
gapperccolremoveseq=95

# Parse command-line options using getopts
while getopts "hf:a:c:o:" opt; do
  case "$opt" in
    h)
      usage
      ;;
    f)
      fasta="$OPTARG"
      ;;
    a)
      aligner="$OPTARG"
      ;;
    c)
      canona="$OPTARG"
      ;;
    o)
      outdir="$OPTARG"
      ;;
    :) # Handle missing argument for options that require it
      echo "Error: Option -$OPTARG requires an argument." >&2
      usage
      exit 1
      ;;
    \?) # Handle invalid options
      echo "Error: Invalid option -$OPTARG" >&2
      usage
      exit 1
      ;;
  esac
done

# Set parameter
[ ! -v fasta ] && usage
input=`basename ${fasta%.*}`
[ ! -v outdir ] && outdir="denovo-structure-$input-$aligner-$canona"
[ ! -d $outdir ] && mkdir -p $outdir
qcout="$outdir/$input.qcout"


### functions

#remove sequences with less than $canona percentage of canonical basepairs in consensus structure
function remove_seq_lt_canon_ss_cons () {
	local input=$1
	local canon=$2

	/home/users/seemann/projects-git/parse-stockholm/parsestockholm.py -f $input.sto -n $canon > "$input-canonicalbp$canon.tsv"
	/home/users/seemann/projects-git/parse-stockholm/parsestockholm.py -f $input.sto -n $canon -o FASTA > "$input-canonicalbp$canon.fasta"
	/home/users/seemann/projects-git/parse-stockholm/parsestockholm.py -f $input.sto -n $canon -o STOCKHOLM > "$input-canonicalbp$canon.sto"
}

#run R-scape two-set statistical test (one test for annotated basepairs, another for all other pairs)
function rscape_twotest () {
	local input=$1
	local id=$2
	local rdir=$3
	[ ! -d $rdir/power ] && mkdir -p $rdir/power
	[ ! -d $rdir/sorted.cov ] && mkdir -p $rdir/sorted.cov
	[ ! -d $rdir/fig ] && mkdir -p $rdir/fig

	#R-scape -s --GTp --C16 --nofigures $input > /dev/null
	R-scape -s --GTp --C16 $input > /dev/null
  mv ${id}*.power $rdir/power
  mv ${id}*.sorted.cov $rdir/sorted.cov
	mv ${id}*.dplot.svg $rdir/fig
	mv ${id}*.R2R.sto.* $rdir/fig
	rm ${id}*
}

#diagnostic of alignment quality
function align_qc () {
	local input=$1
	local rdir=$2
	[ ! -d $rdir ] && mkdir -p $rdir
	local id=`basename $input`

	#remove gaps and consensus structure
	sed 's/-//g' $input.fasta | sed 's/|/-/g' | head -n -2 > $rdir/temp.$id-wogaps.fasta

	#number of sequences
	local noseq=$(grep ">" $rdir/temp.$id-wogaps.fasta | wc -l)

	#consensus structure
	local petfold=$(tail -1 $input.fasta)

	#average pairwise sequence identity
	local apsi=$(seqidentity.pl <(head -n -2 $input.fasta) 100)

	#G+C content
	module load EMBOSS/6.6.0
	local gc=$(infoseq -nohead -only -name -length -pgc $rdir/temp.$id-wogaps.fasta 2> /dev/null | awk 'BEGIN{sum=0}{sum=sum+$3}END{print sum/NR}')

	#R-scape significant basepair and power
	module load rscape/2.0.0.j
	rscape_twotest $input.sto $id $rdir/rscape-v2_0_0_j_twoset

	#parse rscape results
  	if [ -s $rdir/rscape-v2_0_0_j_twoset/sorted.cov/${id}_1.sorted.cov ]; then
   		local covbp=$(awk 'BEGIN{OFS="\t"}$1=="*"{if(bp==""){bp=$2"-"$3}else{bp=bp","$2"-"$3}}END{print bp}' $rdir/rscape-v2_0_0_j_twoset/sorted.cov/${id}_1.sorted.cov)
  	fi
  	if [ -s $rdir/rscape-v2_0_0_j_twoset/power/${id}_1.power ]; then
   		local power=( $(awk 'BEGIN{OFS="\t"}$2=="BPAIRS"{if($3=="observed"){obsbp=$6}else if($3=="expected"){expbp=$6;expbp_conf=$8}else{allbp=$3}}END{print allbp,expbp,expbp_conf,obsbp}' $rdir/rscape-v2_0_0_j_twoset/power/${id}_1.power) )
  	fi

	#average mutual information content across all columns, average relative entropy (Kullback-Leibler divergence) across all columns, and number of compensatory basepair changes and covarying basepairs in the consensus structure
	local mec=( $(/home/users/seemann/projects-git/parse-stockholm/parsestockholm.py -f $input.sto -mec | awk -F "\t" 'BEGIN{OFS="\t"}{if(NR==1){mic=$2};if(NR==2){entr=$2};if(NR==4){print mic,entr,$1,$2}}') )

	#summary of alignment QC
	echo -e $id "\t" $noseq "\t" $petfold "\t" $apsi "\t" $gc "\t" ${mec[0]} "\t" ${mec[1]} "\t" ${mec[2]} "\t" ${mec[3]} "\t" ${power[0]} "\t" ${power[1]} "\t" ${power[2]} "\t" ${power[3]} "\t" $covbp >> $qcout
}

#add PETfold consensus structure to fasta
function add_petfold_ss_cons_to_fasta () {
	local in=$1

  module load PETfold/2.2
	sscons=$(PETfold -f $in -g 0.5 --war | tail -1)
	echo -e ">SS_cons\n"$sscons >> $in
}

#convert fasta to stockholm
function fa2sto () {
	local in=$1
	local out=$2

	local tempf=`mktemp`
	fasta2tab $in | cut -f1,2 > $tempf

	local noseq=$(head -n -1 $tempf | wc -l)

	echo "# STOCKHOLM 1.0" > $out
	echo "#=GF AU Realignment" >> $out
	echo -e "#=GF SQ "$noseq"\n" >> $out
	awk '{if($1=="SS_cons"){printf "%-50s %s\n", "#=GC SS_cons", $2}else{printf "%-50s %s\n", $1, $2}}' $tempf >> $out
	echo "//" >> $out
}

#create CM with cmbuild and realign with cmalign
function realign_cmbuild_cmalign () {
	local input=$1
	local rdir=$2
	[ ! -d $rdir ] && mkdir -p $rdir
	local id=`basename $input`

  #remove gaps and consensus structure
  sed 's/-//g' $input.fasta | head -n -2 > $rdir/temp.$id-wogaps.fasta

  module load infernal/1.1.4
  cmbuild $rdir/$id.cmbuild.cm $input.sto > /dev/null
  cmalign $rdir/$id.cmbuild.cm $rdir/temp.$id-wogaps.fasta > $rdir/$id.cmbuild_cmalign.sto
  /home/users/seemann/projects-git/parse-stockholm/parsestockholm.py -f $rdir/$id.cmbuild_cmalign.sto -u 0 -o FASTA > $rdir/$id.cmbuild_cmalign.fasta
}

#remove sequences with a nucleotide in an alignment column with more than $gapperccolremoveseq percentage of gaps
function remove_seq_in_gapcolumn () {
  local input=$1
  local rdir=$2
  [ ! -d $rdir ] && mkdir -p $rdir
  local gapperccolremoveseq=95
  [[ $3 ]] && gapperccolremoveseq=$3
  local id=`basename $input`

  #clean alignment by removing sequences that introduce gap columns supported by less than 5% of all sequences
  /home/users/seemann/projects-git/parse-stockholm/parsestockholm.py -f $input.sto -q $gapperccolremoveseq -o FASTA > $rdir/$id-cleaned.fasta
  fa2sto $rdir/$id-cleaned.fasta $rdir/$id-cleaned.sto
}

#remove lonely basepairs
function remove_lonely_basepairs () {
  local input=$1
  local rdir=$2
  [ ! -d $rdir ] && mkdir -p $rdir
  local id=`basename $input`

  /home/users/seemann/projects-git/parse-stockholm/parsestockholm.py -f $input.sto -l -o FASTA > $rdir/$id-nolonely-bp.fasta
  fa2sto $rdir/$id-nolonely-bp.fasta $rdir/$id-nolonely-bp.sto
}

#sort sequences in alignment by sequence identity
function sort_by_sequence_identity () {
  local input=$1
  local rdir=$2
  [ ! -d $rdir ] && mkdir -p $rdir
  local id=`basename $input`

  /home/users/seemann/projects-git/parse-stockholm/parsestockholm.py -f $input.sto -w -o FASTA > $rdir/$id-sorted.fasta
  fa2sto $rdir/$id-sorted.fasta $rdir/$id-sorted.sto
}

#update PETfold consensus structure and create CM with cmbuild and realign with cmalign
function realign_petfold_cmbuild_cmalign () {
	local input=$1
	local rdir=$2
	[ ! -d $rdir ] && mkdir -p $rdir
	local id=`basename $input`

	#add PETfold consensus structure to fasta
	head -n -2 $input.fasta > $rdir/$id.petfold.fasta
	add_petfold_ss_cons_to_fasta $rdir/$id.petfold.fasta
	#convert fasta to stockholm
	fa2sto $rdir/$id.petfold.fasta $rdir/$id.petfold.sto

	#remove gaps and consensus structure
	sed 's/-//g' $input.fasta | head -n -2 > $rdir/temp.$id-wogaps.fasta

	module load infernal/1.1.4
	cmbuild $rdir/$id.petfold_cmbuild.cm $rdir/$id.petfold.sto > /dev/null
	cmalign $rdir/$id.petfold_cmbuild.cm $rdir/temp.$id-wogaps.fasta > $rdir/$id.petfold_cmbuild_cmalign.sto
	/home/users/seemann/projects-git/parse-stockholm/parsestockholm.py -f $rdir/$id.petfold_cmbuild_cmalign.sto -u 0 -o FASTA > $rdir/$id.petfold_cmbuild_cmalign.fasta
}

#realign with hmmer3
function realign_hmmer () {
  local msa=$1
  local fasta=$2
  local hmm=${msa%.*}".hmm"
  module load hmmer3/hmmer-3.1b2

  hmmbuild $hmm $msa
  hmmalign $hmm $fasta > "$hmm.sto"
  stk2tab "$hmm.sto" "$hmm.tab"
  tab2fasta "$hmm.tab" > "$hmm.fasta"
}


### main

#align QC header
echo -e "Id\tNo.seq\tPetfold.ss_cons\tAPSI\tGC\tMIC\tRel.entropy\tTotal.comp.bp.changes\tTotal.covarying.bps\tNo.bp\tExp.covarying.bp\tExp.covarying.bp.conf\tObs.covarying.bps\tList.of.covarying.bps" > $qcout

# multiple sequence alignment: Tcoffee, Tcoffee+hmmer3, Muscle, Muscle+hmmer3, nhmmer
case "$aligner" in
  "Tcoffee")
    module load seqan
    seqan_tcoffee -s $fasta -a dna -o "$outdir/$input.tcoffee.fasta" -m global
    ln -s "$input.tcoffee.fasta" "$outdir/$input.msa.fasta"
    ;;
  "Tcoffee+hmmer3")
    module load seqan
    seqan_tcoffee -s $fasta -a dna -o "$outdir/$input.tcoffee.fasta" -m global
    realign_hmmer "$outdir/$input.tcoffee.fasta" $fasta
    ln -s "$input.tcoffee.hmm.fasta" "$outdir/$input.msa.fasta"
    ;;
  "Muscle")
    module load muscle/5.1
    muscle -align $fasta -output "$outdir/$input.muscle.fasta"
    ln -s "$input.muscle.fasta" "$outdir/$input.msa.fasta"
    ;;
  "Muscle+hmmer3")
    module load muscle/5.1
    muscle -align $fasta -output "$outdir/$input.muscle.fasta"
    realign_hmmer "$outdir/$input.muscle.fasta" $fasta
    ln -s "$input.muscle.hmm.fasta" "$outdir/$input.msa.fasta"
    ;;
  "nhmmer")
    module load hmmer3/hmmer-3.1b2
    fasta2tab $fasta | awk 'NR==1{print ">"$1; print $2}' > "$outdir/$input.ref.fa"
    fasta2tab $fasta | awk 'NR>1{print ">"$1; print $2}' > "$outdir/$input.hits.fa"
    nhmmer --rna --incE 1e-10 --tblout "$outdir/$input.nhmmer-1.tbl" --noali -A "$outdir/$input.nhmmer-1.sto" "$outdir/$input.ref.fa" "$outdir/$input.hits.fa" > /dev/null
    nhmmer --rna --incE 1e-10 --tblout "$outdir/$input.nhmmer-2.tbl" --noali -A "$outdir/$input.nhmmer-2.sto" "$outdir/$input.nhmmer-1.sto" "$outdir/$input.hits.fa" > /dev/null
    nhmmer --rna --incE 1e-5 --tblout "$outdir/$input.nhmmer-3.tbl" --noali -A "$outdir/$input.nhmmer-3.sto" "$outdir/$input.nhmmer-2.sto" "$outdir/$input.hits.fa" > /dev/null
    /home/users/seemann/projects-git/parse-stockholm/parsestockholm.py -f "$outdir/$input.nhmmer-3.sto" -o FASTA > "$outdir/$input.nhmmer-3.fasta"
    ln -s "$input.nhmmer-3.fasta" "$outdir/$input.msa.fasta"
    ;;
esac

# check if at least 3 sequences exist - if not then exit
if [ `grep "^>" "$outdir/$input.msa.fasta" | wc -l` -le 3 ]; then
  echo "Blastn has found less then 2 hits.";
  exit 1
fi

# add PETfold consensus structure to fasta
cp -p "$outdir/$input.msa.fasta" "$outdir/$input.msa.petfold.fasta"

add_petfold_ss_cons_to_fasta "$outdir/$input.msa.petfold.fasta"
#convert fasta to stockholm
fa2sto "$outdir/$input.msa.petfold.fasta" "$outdir/$input.msa.petfold.sto"

# 1st diagnostic of alignment quality
align_qc "$outdir/$input.msa.petfold" $outdir

# create CM with cmbuild and realign with cmalign
realign_cmbuild_cmalign "$outdir/$input.msa.petfold" $outdir

# 2nd diagnostic of alignment quality
align_qc "$outdir/$input.msa.petfold.cmbuild_cmalign" $outdir

# remove sequences with less than $canona percentage of canonical basepairs in consensus structure
remove_seq_lt_canon_ss_cons "$outdir/$input.msa.petfold.cmbuild_cmalign" $canona

# 3rd diagnostic of alignment quality
align_qc "$outdir/$input.msa.petfold.cmbuild_cmalign-canonicalbp$canona" $outdir

# remove sequences with a nucleotide in an alignment column with more than $gapperccolremoveseq percentage of gaps
remove_seq_in_gapcolumn "$outdir/$input.msa.petfold.cmbuild_cmalign-canonicalbp$canona" $outdir $gapperccolremoveseq

# remove lonely basepairs
remove_lonely_basepairs "$outdir/$input.msa.petfold.cmbuild_cmalign-canonicalbp$canona-cleaned" $outdir

# sort alignment by sequence identity using the (central) sequence with shortest distance to all sequences as references
sort_by_sequence_identity "$outdir/$input.msa.petfold.cmbuild_cmalign-canonicalbp$canona-cleaned-nolonely-bp" $outdir

# 4th diagnostic of alignment quality
align_qc "$outdir/$input.msa.petfold.cmbuild_cmalign-canonicalbp$canona-cleaned-nolonely-bp-sorted" $outdir

# update PETfold consensus structure and create CM with cmbuild and realign with cmalign
realign_petfold_cmbuild_cmalign "$outdir/$input.msa.petfold.cmbuild_cmalign-canonicalbp$canona-cleaned-nolonely-bp-sorted" $outdir

# 5th diagnostic of alignment quality
align_qc "$outdir/$input.msa.petfold.cmbuild_cmalign-canonicalbp$canona-cleaned-nolonely-bp-sorted.petfold_cmbuild_cmalign" $outdir

#clean
rm -rf $outdir/temp.*
