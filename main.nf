#!/usr/bin/env nextflow

/*
 *  Purpose of this nextflow pipeline:
 *    - Take raw fastq files
 *    - Clip adapter sequence
 *    - Trim 4N from 5' and 3'
 *    - Collapse
 *    - Map to genome (for counting sRNA reads)
 *    - Align non-collapsed to miRNA hairpins (for QualityScore cutoffs)
 */

/*
 *  Considerations / ToDo:
 *    - Add preprocessing steps:
 *        - Demultiplexing
 *        - Adapter clipping & Trimming
 *    - Write out json?
 *    - Add R postprocessing!
 */

version = 0.1.1

// Configurable variables -- default values
params.genome = false
// Get genome files depending on --genome matched in impimba.config, if the genome is found and set on CLI.
params.index         = params.genome ? params.genomes[ params.genome ].bowtie ?: false : false
params.hairpin       = params.genome ? params.genomes[ params.genome ].hairpin ?: false : false
params.wholeGenome   = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.genomeAnno    = params.genome ? params.genomes[ params.genome ].genomeAnno ?: false : false
params.conversionIdx = params.genome ? params.genomes[ params.genome ].ucscNames ?: false : false
params.rdna          = params.genome ? params.genomes[ params.genome ].ribosome ?: false : false
params.repeats       = params.genome ? params.genomes[ params.genome ].repeats ?: false : false
params.saveReference = false
params.name = "miRNA-Seq Best practice"
params.outdir = './results'
// if params.mismatches is null, assign 3, otherwise assign mismatches
mismatches = params.mismatches ? ${params.mismatches} ?: 0 : 3

// Check that we have a hairpin and wholeGenome reference Fasta!
if( !params.hairpin ){
    exit 1, "Missing hairpin reference indexes! Is --genome specified?"
}

// Logging
log.info "==========================================="
log.info " Ameres Lab sRNAseq pipeline v${version}"
log.info "==========================================="
log.info "Reads                : ${params.reads}"
if (params.genome) log.info "Genome               : ${params.genome}"
log.info "Genome Annotation    : ${params.genomeAnno}"
log.info "miRBase hairpin      : ${params.hairpin}"
log.info "Index                : ${params.index}"
log.info "Mismatches           : ${mismatches}"
log.info "Current user         : $USER"
log.info "Current path         : $PWD"
log.info "Script dir           : $baseDir"
log.info "Output dir           : ${params.outdir}"
log.info "Config Profile       : ${workflow.profile}"
log.info "==========================================="

/*
 *  SETUP -- Error checking
 */

if( params.hairpin ){
  hairpins = file(params.hairpin)
  if( !hairpins.exists() ) exit 1, "Hairpin file not found: ${params.hairpin}"
}
if( params.wholeGenome ){
  genomeFasta = file(params.wholeGenoeme)
  if( !genomeFasta.exists() ) exit 1, "Genome Fasta file not found: ${params.wholeGenome}"
}

// Do we really need this? We'll be generating the index every time...
if( params.index ){
  bt_index = file("${params.index}.1.ebwt")
  bt_indices = Channel.fromPath( "${params.index}*" ).toList()
}

/* This allows passing wild card tagged files on CLI:
 * `--reads <file*.x>`
 */
Channel
  .fromPath( params.reads )
  .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
  .into { raw_reads }

genomeParts = Channel.from('ribosome', 'mitochondria', 'genome_cleaned')

/*
 * PREPROCESSING - Extract Flybase Annotations:
 *   - Extract anything with Flybase in Col2 and name it
 *     according to the annotaiton in it's parent ID line
 */
process prepareGenome {
  input:
  file genome from genomeFasta
  file rRnaPrecursor from rRNA
  file assReport from assemblyReport

  output:
  file '*.fa' into riboMitoGenomeFa

  script:
  """
  # Change Fasta to tab separated format
  gunzip -c $genome | \\
    awk '/^>/&&NR>1{printf \"%s\", /^>/ ? \$0\" \":\$0}' | \\
    tr ' ' '\\t' > genome_fb-names.txt

  # Get flybase[tab]ucscNames
  awk '(\$0 !~ /^#/){printf \"%s\\t%s\", \$1 \$NF}' $assReport > nameConv.tsv

  # Write out ribosome, mitochondria and
  # exchange other names to ucscNames
  awk -v names=\"nameConv.tsv\" \\
      -v ribo=\"ribosome.txt\" \\
      -v mito=\"mito.fa\" \\
      'BEGIN{
          while((getline L < names) > 0) {
            if (\$NF ~ \"chrM\") nameArray[\"mitochondrion_genome\"] = \$2
            if (\$1 ~ \"rDNA\") nameArray[\"rDNA\"] = \"rDNA\"
            else nameArray[\$1] = \$2
          }
        }
        {
          if (\$1 ~ /mitochondrion/) printf \">%s\\n%s\", nameArray[\$1], \$2 >> ribo
          else if (\$1 ~ /rDNA/) printf \">%s\\n%s\", nameArray[\$1], \$2 >> mito
          else printf \">%s\\n%s\", nameArray[\$1], \$2
        }' > genome_noMito-noRibo_incl-Y.fa

  cat ribosome.txt $rRnaPrecursor > ribosome.fa
  """
}
process prepareAnnos {
  beforeScript 'installdeps.R'

  input:
  file gff from genomeAnno
  file nameConv from assemblyReport
  file rRNA

  output:
  file "fullGenomeAnno.bed" into genomeAnno

  script:
  """
  # Extract Flybase Annos
  prepareref.R -g $gff -n $nameConv -o flybase.bed

  # Add rRNA Annotation

  """
}

process trim_adapter {
  tag "$reads"

  input:
  file reads from raw_reads

  output:
  file "*.adapter_clipped.fq.gz" into adapter_clipped
  file "*.trim_report.txt" into trim_results

  script:
  prefix = reads.toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
  """
  cutadapt \\
    -m 26 \\
    -M 38 \\
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \\
    -o ${prefix}.adapter_clipped.fq.gz \\
    $reads > ${prefix}.trim_report.txt
  """
}

/*
 *  STEP 2: Trim 4-N from 5' and 3'-ends
 */

process trim_4N {
  tag "$reads"

  input:
  file reads from adapter_clipped

  output:
  file "*.trimmed.fq.gz" into trimmed_reads

  script:
  prefix = reads.toString() - ".adapter_clipped.fq.gz"
  """
  /home/imba/reichholf/bin/seqtk trimfq \\
    -b 4 \\
    -e 4 \\
    $reads > ${prefix}.trimmed.fq

  gzip ${prefix}.trimmed.fq
  """
}

/*
 * STEP 3: Align
 */
process bowtie_hairpins {
  tag "$reads"
  publishDir "${params.outdir}/bowtie/ext_hairpins", mode: "copy", pattern: '*.TCtagged_hairpin.bam'

  input:
  file reads from trimmed_reads
  file index from hairpin_index

  output:
  file '*.TCtagged_hairpin.bam' into hairpin_aligned

  script:
  input_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
  prefix = reads.toString() - ~/(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
  """
  bowtie \\
    -a --best --strata \\
    -v $mismatches \\
    -S \\
    -p ${task.cpus} \\
    $input_base \\
    -q <(zcat $reads) | \\
    /groups/ameres/Reichholf/bioinf/scripts/arrayTagTCreads.awk | \\
    samtools view -bS - > ${prefix}.TCtagged_hairpin.bam
  """
}

/*
 *  STEP 4: Postprocess (Sorting/counting)
 */

def wrap_hairpin = { file ->
  if (file.contains("hairpin")) return "ext_hairpins/$file"
}

process post_alignment {
  publishDir "${params.outdir}/bowtie", mode: "copy", saveAs: wrap_hairpin

  input:
  file input from hairpin_aligned

  output:
  file "${input.baseName}.count" into hairpin_counts
  file "${input.baseName}.sorted.bam" into hairpin_sorted_bam
  file "${input.baseName}.sorted.bam.bai" into hairpin_sorted_bai

  script:
  """
  samtools sort ${input.baseName}.bam -o ${input.baseName}.sorted.bam
  samtools index ${input.baseName}.sorted.bam
  samtools idxstats ${input.baseName}.sorted.bam > ${input.baseName}.count
  """
}
