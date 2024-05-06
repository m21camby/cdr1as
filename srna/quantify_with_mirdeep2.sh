#!/bin/bash

# Small RNA-seq analyses for Cledi quantification of known miRNAs with
# miRDeep2 bash script
#
# Copyright: Marcel Schilling with Seung
#


#######################
# general information #
#######################

# file:    quantify_with_mirdeep2.sh
# created: 2018-11-22
# author:  Marcel Schilling
# purpose: quantify known miRNAs in sRNA-seq data using miRDeep2


######################################
# change log (reverse chronological) #
######################################

# 2018-11-22: initial version


#########
# usage #
#########

# Define usage.
USAGE="$0 <min-read-length> <species-id> <threads> <genome-bowtie-indices> \\
  <mirbase-hairpin-fa> <mirbase-mature-fa> <reads-fastq> <results-tsv> \\
  <mapper-log> <mapper-bowtie-log> <mapper-mapper-log> <quantifier-log>"

# Check for proper usage.
[ $# -eq 12 ] || { printf '%s\n' "usage: $USAGE" >&2; exit 1; } 


##############
# parameters #
##############

# Get minimum length to consider small RNA-seq read for mapping.
MIN_READ_LENGTH=$1

# Get miRBase species ID to pass to miRDeep2.
SPECIES_ID="$2"

# Get number of threads/processes to run at the same time.
THREADS=$3


#########
# paths #
#########

# Get (relative) path of bowtie index (without file extension) for
# genome reference.
GENOME_BOWTIE_INDICES="$4"

# Get (relative) path of (compressed) FASTA file with miRNA hairpin
# sequences extracted from miRBase.
MIRBASE_HAIRPIN_FA="$5"

# Get (relative) path of (compressed) FASTA file with mature miRNA
# sequences extracted from miRBase.
MIRBASE_MATURE_FA="$6"

# Get (relative) path of (compressed) FASTQ file with preprocessed sRNA-seq
# reads.
READS_FASTQ="$7"

# Get (relative) path of compressed TSV file to write results to.
RESULTS_TSV="$8"

# Get (relative) path of TXT file to write miRDeep2 mapper log to.
MAPPER_LOG="$9"

# Get (relative) path of TXT file to write miRDeep2 mapper bowtie log to.
MAPPER_BOWTIE_LOG="${10}"

# Get (relative) path of TXT file to write miRDeep2 mapper mapper log to.
MAPPER_MAPPER_LOG="${11}"

# Get (relative) path of TXT file to write miRDeep2 quantification log to.
QUANTIFIER_LOG="${12}"


#################
# bash settings #
#################

# Exit with error as soon as any command exists with error.
set -e


####################
# absolutify paths #
####################

# Get absolute path of first bowtie index file based on bowtie index path.
GENOME_BOWTIE_INDEX_FILE_1=$(realpath "$GENOME_BOWTIE_INDICES.1.ebwt")

# Derive absolute path of bowtie index from first index file.
GENOME_BOWTIE_INDICES=${GENOME_BOWTIE_INDEX_FILE_1%.1.ebwt}

# Get absolute path of FASTA file with miRNA hairpins extracted from miRBase.
MIRBASE_HAIRPIN_FA=$(realpath "$MIRBASE_HAIRPIN_FA")

# Get absolute path of FASTA file with mature miRNAs extracted from miRBase.
MIRBASE_MATURE_FA=$(realpath "$MIRBASE_MATURE_FA")

# Get absolute path of FASTQ file with preprocessed sRNA-seq reads.
READS_FASTQ=$(realpath "$READS_FASTQ")

# Create (empty) TSV file for miRDeep2 results and get absolute path.
touch "$RESULTS_TSV"
RESULTS_TSV=$(realpath "$RESULTS_TSV")

# Create (empty) miRDeep2 mapper log file and get absolute path.
touch "$MAPPER_LOG"
MAPPER_LOG=$(realpath "$MAPPER_LOG")

# Create (empty) miRDeep2 mapper bowtie log file and get absolute path.
touch "$MAPPER_BOWTIE_LOG"
MAPPER_BOWTIE_LOG=$(realpath "$MAPPER_BOWTIE_LOG")

# Create (empty) miRDeep2 mapper mapper log file and get absolute path.
touch "$MAPPER_MAPPER_LOG"
MAPPER_MAPPER_LOG=$(realpath "$MAPPER_MAPPER_LOG")

# Create (empty) miRDeep2 quantification log file and get absolute path.
touch "$QUANTIFIER_LOG"
QUANTIFIER_LOG=$(realpath "$QUANTIFIER_LOG")


############
# commands #
############

# Define command to use to create a temporary directory.
CREATE_TEMP_DIRECTORY='mktemp --directory'

# Define command to use to uncompress GZIP files (to STDOUT).
ZCAT='zcat'

# Define command to use to uncompress file if necessary compressed (to STDOUT).
DECOMPRESS="$ZCAT --force"

# Define command to use to run awk.
AWK='awk'

# Define awk command to use to convert FASTQ to FASTA format (STDIN to STDOUT).
FASTQ2FASTA_AWK='{print ">" substr($0,2);getline;print;getline;getline}'

# Define awk command to use to adjust FASTA headers for miRDeep2 (STDIN to
# STDOUT).
# The following format is needed for miRDeep2's quantifier module:
# > Please make sure your file is in accordance with the fasta format
# > specifications and does not contain whitespace in IDs or sequences.
ADJUST_FASTA_HEADERS_FOR_MIRDEEP2_AWK='!!(NR%2){NF=1}1'

# Define command to use to change the working directory.
CHANGE_DIRECTORY='cd'

# Define command to use to run miRDeep2's mapper module.
MIRDEEP2_MAPPER='mapper.pl'

# Define parameters to use for miRDeep2's mapper module.
# The parameters where chosen based on the folowing email:
# -------------------- citing the following email -------------------- #
# From: Filippos Klironomos <Filippos.Klironomos@mdc-berlin.de>
# To: Marcel Schilling <marcel.schilling@mdc-berlin.de>
# Date: Tue, 3 Feb 2015 17:08:24 +0100
# Subject: mirdeep2
# -------------------------- citation begin -------------------------- #
#  $1 = {SampleID} (without .fastq part, i.e. {SampleID}.fastq should be
#       the full name of the library)
#  $2 = {Adaptor} (exactly what will be fed to mapper.pl using the -k
#       flag)
#  $3 = {Reference} full pathname of reference genome FASTA file and
#       index-file-prefix
#  $4 = {Hairpins} full pathname of miRBase FASTA file for precursors
#  $5 = {Mature} full pathname of miRBase FASTA file for mature
#  $6 = {Species} [Human,Mouse]
#
#  mapper.pl config -d -e -h -j -k $2 -l 18 -m -p $3 \
#                   -s reads_collapsed.fa -t reads_vs_genome.arf -v \
#                   -o 12  &> mapper_summary.log
#
#  quantifier.pl -p hsa.hairpin.miRBase20.fa -m hsa.mature.miRBase20.fa
#                -r reads_collapsed.fa -P -t hsa -y now
# --------------------------- citation end --------------------------- #
# Note that Filippos used `-e` (see above) for FASTQ input, but here
# FASTA (`-c`) input is used.
MIRDEEP2_MAPPER_OPTIONS="-c -j -l $MIN_READ_LENGTH -m -v -o $THREADS -p '$GENOME_BOWTIE_INDICES'"

# Define command to use to rename files.
RENAME='mv'

# Define command to use to run miRDeep2's quantifier module.
# Note that Filippos didn't use the `-d` and `-j` parameter to not
# create PDFs or an `output.mrd` file but but as those would be deleted
# here anyhow, there is no need to generate them in the first place.
MIRDEEP2_QUANTIFIER="quantifier.pl -P -d -j -t $SPECIES_ID"

# Define command used to run pigz.
PIGZ="pigz --best --processes $THREADS"

# Define command to use to compress data (STDIN to STDOUT).
COMPRESS="$PIGZ"

# Define command to use to remove directories.
REMOVE='rm --force --recursive --'


##########################
# temporary output files #
##########################

# Create temporary output directory.
OUTPUT_DIRECTORY=$($CREATE_TEMP_DIRECTORY)

# Create absolute path of temporary output directory.
OUTPUT_DIRECTORY=$(realpath "$OUTPUT_DIRECTORY")

# Define absolute path of temporary FASTA file with input reads.
RAW_READS_FASTA="$OUTPUT_DIRECTORY/reads.fa"

# Define absolute path of temporary FASTA file with miRNA hairpins.
HAIRPIN_FA="$OUTPUT_DIRECTORY/hairpin.fa"

# Define absolute path of temporary FASTA file with mature miRNAs.
MATURE_FA="$OUTPUT_DIRECTORY/mature.fa"

# Define absolute path of temporary file with processed reads.
PROCESSED_READS_FILE="$OUTPUT_DIRECTORY/reads.processed"

# Define absolute path of temporary ARF file with alignments.
ALIGNMENTS_ARF="$OUTPUT_DIRECTORY/alignments.arf"

# Define filename stem (without leading or final dots (`.`s)) for
# miRDeep2 result file.
RESULTS_STEM="results"


######################
# setup miRDeep2 run #
######################

# Decompress FASTQ file with reads if necessary and convert to miRDeep2
# compatible FASTA file.
$DECOMPRESS \
< "$READS_FASTQ" \
| $AWK "$FASTQ2FASTA_AWK" \
| $AWK "$ADJUST_FASTA_HEADERS_FOR_MIRDEEP2_AWK" \
> "$RAW_READS_FASTA"

# Decompress FASTA file with miRNA hairpins if necessary.
$DECOMPRESS \
< "$MIRBASE_HAIRPIN_FA" \
> "$HAIRPIN_FA"

# Decompress FASTA file with mature miRNAs if necessary.
$DECOMPRESS \
< "$MIRBASE_MATURE_FA" \
> "$MATURE_FA"

# Change working directory to temporary output directory to ensure all
#miRDeep2 output file are create inside this run-specific directory.
$CHANGE_DIRECTORY "$OUTPUT_DIRECTORY"


#################
# mapper module #
#################

# Run miRDeep2 mapper module to generate temporary file with processed
# reads. Write log output to mapper log.
$MIRDEEP2_MAPPER \
  "$RAW_READS_FASTA" \
  $MIRDEEP2_MAPPER_OPTIONS \
  -s "$PROCESSED_READS_FILE" \
  -t "$ALIGNMENTS_ARF" \
2> "$MAPPER_LOG"


# Move bowtie log file to target miRDeep2 mapper bowtie log.
$RENAME \
  "$OUTPUT_DIRECTORY/bowtie.log" \
  "$MAPPER_BOWTIE_LOG"

# Move mapper log file to target miRDeep2 mapper mapper log.
$RENAME \
  "$OUTPUT_DIRECTORY/mapper_logs/"mapper.log_* \
  "$MAPPER_MAPPER_LOG"


#####################
# quantifier module #
#####################

# Run miRDeep2 quantifier using temporary file with processed reads.
# Write log output to quantifier log.
$MIRDEEP2_QUANTIFIER \
  -p "$HAIRPIN_FA" \
  -m "$MATURE_FA" \
  -r "$PROCESSED_READS_FILE" \
  -y "$RESULTS_STEM" \
2> "$QUANTIFIER_LOG"

# Compresss output TSV file and save it to target path.
$COMPRESS \
< "$OUTPUT_DIRECTORY/miRNAs_expressed_all_samples_$RESULTS_STEM.csv" \
> "$RESULTS_TSV"


###########
# cleanup #
###########

# Leave temporary output directory.
$CHANGE_DIRECTORY \
  -

# Remove temporary output directory.
$REMOVE \
  "$OUTPUT_DIRECTORY"
