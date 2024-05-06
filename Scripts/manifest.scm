(use-modules (bimsb packages staging))

(define packages '(
  "fastqc@0.11.5"
  "bwa@0.7.17"
  "bowtie@2.3.4.3"
  "python"
  "python2"
  "python2-pysam"
  "python2-numpy"
  "bedtools"
  "python-pandas@0.25.2"
  "python-numpy@1.17.3"
  "star"
  "subread"
  "samtools"
  "trim-galore"
  "kallisto"
  "rseqc@2.6.1"
  "htseq@0.9.1"
  "bbmap"))


  
(specifications->manifest packages)
