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
 *    - Add preprocessing steps?
 *        - Demultiplexing
 *        - Adapter clipping & Trimming
 *    - Write out json?
 *    - Get hairpin (FASTA) from genome + 20nt downstream
 *    - Build Hairpin index
 *    - align to hairpin
 *    - Add R postprocessing!
 *      - Get max positions
 *      - Get mutation rates
 *      - Get T>C read count
 *      - Get steady state read count
 *      - if available (json/align and count myself?) normalise read count to ppm
 *        if not: norm to miRNA(!) mapped reads
 */

version = "0.2.2"

// Configurable variables -- default values
params.genome = false
params.mismatches = 3
// Get genome files depending on --genome matched in impimba.config, if the genome is found and set on CLI.
// params.index         = params.genome ? params.genomes[ params.genome ].bowtie ?: false : false
params.genomeFasta   = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.genomeAnno    = params.genome ? params.genomes[ params.genome ].genomeAnno ?: false : false
params.mirArmAnno    = params.genome ? params.genomes[ params.genome ].mirArmAnno ?: false : false
params.conversionIdx = params.genome ? params.genomes[ params.genome ].ucscNames ?: false : false
params.rdna          = params.genome ? params.genomes[ params.genome ].ribosome ?: false : false
params.repeats       = params.genome ? params.genomes[ params.genome ].repeats ?: false : false
params.saveReference = false
params.name          = "miRNA-Seq Best practice"
params.outdir        = './results'
// if params.mismatches is null, assign 3, otherwise assign mismatches
mismatches           = params.mismatches ?: 3

// Check that we have a hairpin and wholeGenome reference Fasta!
if( !params.genomeAnno ){
    exit 1, "Missing hairpin reference indexes! Is --genome specified?"
}

// Logging
log.info "==========================================="
log.info " Ameres Lab sRNAseq pipeline v${version}"
log.info "==========================================="
log.info "Reads                : ${params.reads}"
if (params.genome) log.info "Genome               : ${params.genome}"
log.info "Genome Annotation    : ${params.genomeAnno}"
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

genomeFastaFile = params.genomeFasta ? file(params.genomeFasta) : null
genomeAnno      = params.genomeAnno ? file(params.genomeAnno) : null
mirArmAnno      = params.mirArmAnno ? file(params.mirArmAnno) : null

if( !params.genomeFasta || !genomeFastaFile.exists() ) {
   exit 1, "Genome Fasta file not found: ${params.genomeFasta}. Please download from flybase."
}

if ( !params.genomeAnno || !genomeAnno.exists() ) {
  exit 1, "Genome annotation file ${params.genomeAnno} not found. Please download from flybase."
}

if ( !params.mirArmAnno || !mirArmAnno.exists() ) {
  exit 1, "miRNA Arm annotation file ${params.mirArmAnno} not found."
}

/* This allows passing wild card tagged files on CLI:
 * `--reads <file*.x>`
 */
Channel
  .fromPath( params.reads )
  .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
  .set { rawReads }

// Required for genomic alignment for hierarchical counting
// process prepareGenome {
//   input:
//   file genome from genomeFasta
//   file rRnaPrecursor from rRNA
//   file assReport from assemblyReport

//   output:
//   file '*.fa' into riboMitoGenomeFa

//   script:
//   """
//   # Change Fasta to tab separated format
//   gunzip -c $genome | \\
//     awk '/^>/&&NR>1{printf \"%s\", /^>/ ? \$0\" \":\$0}' | \\
//     tr ' ' '\\t' > genome_fb-names.txt

//   # Get flybase[tab]ucscNames
//   awk '(\$0 !~ /^#/){printf \"%s\\t%s\", \$1 \$NF}' $assReport > nameConv.tsv

//   # Write out ribosome, mitochondria and
//   # exchange other names to ucscNames
//   awk -v names=\"nameConv.tsv\" \\
//       -v ribo=\"ribosome.txt\" \\
//       -v mito=\"mito.fa\" \\
//       'BEGIN{
//           while((getline L < names) > 0) {
//             if (\$NF ~ \"chrM\") nameArray[\"mitochondrion_genome\"] = \$2
//             if (\$1 ~ \"rDNA\") nameArray[\"rDNA\"] = \"rDNA\"
//             else nameArray[\$1] = \$2
//           }
//         }
//         {
//           if (\$1 ~ /mitochondrion/) printf \">%s\\n%s\", nameArray[\$1], \$2 >> ribo
//           else if (\$1 ~ /rDNA/) printf \">%s\\n%s\", nameArray[\$1], \$2 >> mito
//           else printf \">%s\\n%s\", nameArray[\$1], \$2
//         }' > genome_noMito-noRibo_incl-Y.fa

//   cat ribosome.txt $rRnaPrecursor > ribosome.fa
//   """
// }

// Prepare annotations for hierarchical counting
// process prepareAnnos {
//   beforeScript 'installdeps.R'

//   input:
//   file gff from genomeAnno
//   file nameConv from assemblyReport
//   file rRNA

//   output:
//   file "fullGenomeAnno.bed" into genomeAnno

//   script:
//   """
//   # Extract Flybase Annos
//   prepareref.R -g $gff -n $nameConv -o flybase.bed

//   # Add rRNA Annotation

//   """
// }

process extractHairpins {
  tag "genomePrep"

  input:
  file genomeFastaFile
  file genomeAnno

  output:
  file "hairpin.fa" into hairpinFasta

  script:
  """
  grep pre_miR $genomeAnno | \
    sed -e 's/\"//g' | \
    awk -v FS="\t" '{OFS=FS} {
        split(\$9, info, "; ")
        for (i in info){
          if (info[i] ~ /^transcript/) {
            split(info[i], ts, " ")
            tsid = ts[2]
          }
        }
        strand = \$7
        if (strand == "+") {
          start = \$4
          end = \$5 + 20
        } else {
          start = \$4 - 20
          end = \$5
        }
        print \$1, start, end, tsid, 0, strand
      }' | tr ' ' '\t' > hairpin_plus20nt.bed
    zcat -f ${genomeFastaFile} > genome.fa
    samtools faidx genome.fa
    bedtools getfasta -s -fi genome.fa -bed hairpin_plus20nt.bed > hairpin.fa
  """
}

process makeIndex {
  tag "Hairpins"

  input:
  file hairpinFasta

  output:
  file "hairpin_idx*" into hairpinIndex

  script:
  """
  samtools faidx hairpin.fa
  bowtie-build ${hairpinFasta} hairpin_idx
  """
}

// Clip 3' Adapter using cutadapt
process trim_adapter {
  tag "$reads"

  input:
  file reads from rawReads

  output:
  file "*.adapter_clipped.fq.gz" into adapterClipped
  file "*.trim_report.txt" into trimResults

  script:
  prefix = reads.toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
  """
  cutadapt \
    -m 26 \
    -M 38 \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG \
    -o ${prefix}.adapter_clipped.fq.gz \
    $reads > ${prefix}.trim_report.txt
  """
}

/*
 *  STEP 2: Trim 4-N from 5' and 3'-ends
 */

process trim_4N {
  tag "$acReads"

  input:
  file acReads from adapterClipped

  output:
  file "*.trimmed.fq.gz" into trimmedReads

  script:
  prefix = acReads.toString() - ".adapter_clipped.fq.gz"
  """
  seqtk trimfq \\
    -b 4 \\
    -e 4 \\
    $acReads > trimmedReads.fq

  gzip -c trimmedReads.fq > ${prefix}.trimmed.fq.gz
  """
}

/*
 * STEP 3: Align
 */
process bowtie_hairpins {
  tag "$trimmedReads"

  publishDir "${params.outdir}/bowtie/ext_hairpins", mode: "copy", pattern: '*.TCtagged_hairpin.bam'

  input:
  file trimmedReads
  file index from hairpinIndex

  output:
  file '*.TCtagged_hairpin.bam' into hairpinAligned

  script:
  index_base = index.toString().tokenize(' ')[0].tokenize('.')[0]
  prefix = trimmedReads.toString() - ~/(_trimmed)?(\.fq)?(\.fastq)?(\.gz)?$/
  """
  echo "${trimmedReads}" >> postBowtie.log
  bowtie \
    -a --best --strata \
    -v $mismatches \
    -S \
    -p ${task.cpus} \
    $index_base \
    -q <(zcat $trimmedReads) | \
    arrayTagTCreads.awk | \
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
  tag "$hairpinAligned"

  publishDir "${params.outdir}/bowtie", mode: "copy", saveAs: wrap_hairpin

  input:
  file hairpinAligned

  output:
  file "${hairpinAligned.baseName}.count" into hairpinCounts
  file "${hairpinAligned.baseName}.sorted.bam" into hairpinSorted
  file "${hairpinAligned.baseName}.sorted.bam.bai" into hairpinSortedIndex

  script:
  """
  echo "${hairpinAligned.baseName}" >> postAlignment.log
  samtools sort ${hairpinAligned.baseName}.bam -o ${hairpinAligned.baseName}.sorted.bam
  samtools index ${hairpinAligned.baseName}.sorted.bam
  samtools idxstats ${hairpinAligned.baseName}.sorted.bam > ${hairpinAligned.baseName}.count
  """
}

process writeJson {
  /*
   * INFO:
   *   `file sortedBams from hairpinSorted.collect()` takes care of staging
   *   This way `post_alignment` can still run in parallel, for each sample.
   *   Then `writeJson` will run after that is all done.
   *   Currently we are NOT allowing normalisation by sRNA-mapping reads(!) and will
   *   simply report raw counts instead.
   *   To 'future proof' the R analysis, we have included the option for normalisation
   *   but are simply setting that count to 1000000
   *   Also, we're not handling time at the moment, but this would be an easy fix
   *   once we move to a sample sheet.
   */

  // !!!!
  // Mode needs to change to "copy" if we're going to use the json file later on
  // !!!!
  publishDir "${params.outdir}/counting", mode: "move", pattern: '*.json'

  input:
  file sortedBams from hairpinSorted.toSortedList()
  file counts from hairpinCounts.toSortedList()

  output:
  file "samples.json" into readCountConfig

  script:
  """
  # awk -v name=${sortedBams} \
  #   '{map += \$3; unmapped += \$4} END {printf "%s\t%d\t%d", name, map, unmapped}' \
  #   > countsum.txt
  echo $counts
  """

  """
  #!/usr/bin/env python
  import json
  import os

  bamfiles = [os.path.basename(b) for b in "$sortedBams".split(' ')]
  jsDict = {"base": "${params.outdir}/bowtie",
            "mir.anno": "$mirArmAnno"}

  # with open('countsum.txt', 'r') as counts:
  #   countStats = counts.readline()

  fullName, mapped, unmapped = countStats.split('\t')

  samples = []

  print('Debug-STR :: {}'.format(bamfiles))
  print('Debug-Counts :: {} - {} - {}'.format(fullName, mapped, unmapped))

  for bam in bamfiles:
    nameElements = bam.split('_')
    idx = int(nameElements[0])
    sRNAreads = 1000000
    samples.append({"id": idx,
                    "time": str(idx) + "h",
                    "align": bam,
                    "posFile": str(idx) + "_pos.tsv",
                    "sRNAreads": sRNAreads})

  jsDict["samples"] = samples

  with open('samples.json', 'w') as fp:
    json.dump(jsDict, fp)
  """
}
