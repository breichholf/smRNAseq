#!/usr/bin/env nextflow

/*
 *  Purpose of this nextflow pipeline:
 *    - Take raw fastq files
 *    - Clip adapter sequence
 *    - Trim 4N from 5' and 3'
 *    - Collapse *** IN PROGRESS ***
 *    - Map to genome (for counting sRNA reads) *** IN PROGRESS ***
 *    - Align non-collapsed to miRNA hairpins (for QualityScore cutoffs)
 */

/*
 *  Considerations / ToDo:
 *    - Add preprocessing steps?
 *        - Demultiplexing
 */

version = "0.7.0"

/*
 * Helper functions
 */

String getOutDir(output_type) {
    new File(params.output.get('base_dir', ''),
             params.output.dirs.get(output_type, output_type)).getCanonicalPath()
}

// Configurable variables -- default values
params.genome = false
params.mismatches = 3
// Get genome files depending on --genome matched in impimba.config, if the genome is found and set on CLI.
params.genomeFasta   = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.genomeAnno    = params.genome ? params.genomes[ params.genome ].genomeAnno ?: false : false
params.mirArmAnno    = params.genome ? params.genomes[ params.genome ].mirArmAnno ?: false : false
params.conversionIdx = params.genome ? params.genomes[ params.genome ].ucscNames ?: false : false
params.rdna          = params.genome ? params.genomes[ params.genome ].ribosome ?: false : false
params.repeats       = params.genome ? params.genomes[ params.genome ].repeats ?: false : false
params.saveReference = false
params.name          = "miRNA-Seq Best practice"
// if params.mismatches is null, assign 3, otherwise assign mismatches
mismatches           = params.mismatches ?: 3

// Check that we have a hairpin and wholeGenome reference Fasta!
if( !params.genomeAnno ){
    exit 1, "Missing hairpin reference indexes! Is --genome specified?"
}

// R library locations
params.rlocation = "$HOME/R/nxtflow_libs/"
nxtflow_libs=file(params.rlocation)
nxtflow_libs.mkdirs()

// Adapter initialising
params.adapter = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG"
adapter = params.adapter

// Set default min/max lengths to trim down to
params.min = 26
params.max = 38
minlen = params.min
maxlen = params.max

// --notrim to prevent trimming
params.notrim = false
notrim = params.notrim


// Logging
log.info "==========================================="
log.info " Ameres Lab sRNAseq pipeline v${version}"
log.info "==========================================="
log.info "Reads                : ${params.reads}"
if (params.genome) log.info "Genome               : ${params.genome}"
log.info "Genome Annotation    : ${params.genomeAnno}"
log.info "Genome Fasta         : ${params.genomeFasta}"
log.info "Mismatches           : ${mismatches}"
log.info "sRNA annotated reads : ${params.annoreads}"
log.info "Adapter Sequence     : ${params.adapter}"
log.info "Min / Max / Trim?    : $minlen / $maxlen / $notrim"
log.info "Current user         : $USER"
log.info "Current path         : $PWD"
log.info "Script dir           : $baseDir"
log.info "R location           : ${params.rlocation}"
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

readcounts      = file(params.annoreads)

// We need to put a check in here, to set read counts to 10000000 or so.
if (!readcounts.exists()) exit 1, "Read counts not provided."

/* This allows passing wild card tagged files on CLI:
 * `--reads <file*.x>`
 */
Channel
  .fromPath( params.reads )
  .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
  .set { rawReads }


/*
 *  STEP 1a - Preprocess reference: Extract hairpin loci from genome
 *  Required: GTF file with hairpins. Formatting and naming scheme widely divergent
 *            in different species, currently only works with fly.
 */
process extractHairpins {
  tag "genomePrep"
  publishDir path: getOutDir('ref'), mode: "copy", pattern: 'hairpin*.{fa,bed}'

  input:
  file genomeFastaFile
  file genomeAnno

  output:
  file "hairpin.fa" into hairpinFasta
  file "hairpin_plus20nt.bed" into hairpinBED

  script:
  """
  grep pre_miR $genomeAnno | \
    sed -e 's/\"//g' | \
    awk -v FS="\t" '{OFS=FS} {
        split(\$9, info, "; ")
        for (i in info){
          if (info[i] ~ /^gene_id/) {
            split(info[i], gid, " ")
            geneid = gid[2]
          }
        }
        strand = \$7
        if (strand == "+") {
          start = \$4 - 1    # GTF start coordinate is 1 based, BED is 0-based!
          end = \$5 + 20
        } else {
          start = \$4 - 20
          end = \$5
        }
        print \$1, start, end, geneid, 0, strand
      }' | tr ' ' '\t' > hairpin_plus20nt.bed
    zcat -f ${genomeFastaFile} > genome.fa
    samtools faidx genome.fa
    bedtools getfasta -s -name -fi genome.fa -bed hairpin_plus20nt.bed | \
      sed -e 's/(+)//g;s/(-)//g' > hairpin.fa
  """
}

/*
 *  STEP 1b - Preprocess reference: Make index
 */
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

/*
 *  STEP 2a - Preprocess reads: Clip 3' adapter
 *  Note: Default adapter sequence is identical to Illumina adapter, custom sequence can be supplied.
 *  Default sequence: `AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG`
 */
process trim_adapter {
  tag "$reads"
  publishDir path: getOutDir('reports'), mode: "copy", pattern: '*.trim_report.txt'

  input:
  file reads from rawReads

  output:
  file "*.adapter_clipped.fq.gz" into adapterClipped
  file "*.trim_report.txt" into trimResults

  script:
  prefix = reads.toString() - ~/(\.fq)?(\.fastq)?(\.gz)?$/
  """
  # zcat $reads | paste -d '' - - - - > column.fq
  # Parallelism doesn't work just yet, joined adapter clipped file is not always grouped properly.
  # Default adapter: AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG
  # parallel --pipepart --line-buffer --round-robin -j ${task.cpus} -a column.fq
  # tr '' '' | cutadapt -m 26 -M 38\
  #  -a ${adapter} \
  #  - 2>> ${prefix}.trim_report.txt | gzip > ${prefix}.adapter_clipped.fq.gz
  cutadapt \
    -m ${minlen} \
    -M ${maxlen} \
    -a ${adapter} \
    -o ${prefix}.adapter_clipped.fq.gz \
    $reads > ${prefix}.trim_report.txt
  """
}

/*
 *  STEP 2b - Preprocess reads: Trim random N from 5' and 3' end
 *  Note: if parameter --norand is passed, trimming will be skipped entirel.
 */
process trim_4N {
  /*
   * This relies on seqtk being installed! http://github.com/lh3/seqtk
   * build requirements are only zlib
   */
  tag "$acReads"

  input:
  file acReads from adapterClipped

  output:
  file "*.trimmed.fq.gz" into trimmedReads

  script:
  prefix = acReads.toString() - ".adapter_clipped.fq.gz"
  if (!notrim) {
    """
    seqtk trimfq -b 4 -e 4 $acReads > trimmedReads.fq
    gzip -c trimmedReads.fq > ${prefix}.trimmed.fq.gz
    """
  } else {
    """
    zcat $acReads | gzip -c > ${prefix}.trimmed.fq.gz
    """
  }
}

/*
 * STEP 3: Align
 * Note: We use bowtie v1.2.2, but likely any aligner could be used
 */
process bowtie_hairpins {
  tag "$trimmedReads"

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
  publishDir path: getOutDir('sortedAlignment'), mode: "copy", saveAs: wrap_hairpin

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

/*
 *  STEP 5a - Postprocessing -- Prepare JSON configuration
 */

process writeJson {
  /*
   * INFO:
   *   `file sortedBams from hairpinSorted.collect()` takes care of staging
   *   This way `post_alignment` can still run in parallel, for each sample.
   *   Then `writeJson` will run after that is all done.
   *   Currently we only map, and take sRNAcounts from an externally provided sample file
   */

  publishDir path: getOutDir('json'), mode: "copy", pattern: '*.json'

  input:
  file sortedBams from hairpinSorted.toSortedList()
  file counts from hairpinCounts.toSortedList()
  file readcounts

  output:
  file "samples.json" into readCountConfig

  script:
  absOutDir = getOutDir('sortedAlignment')

  """
  #!/usr/bin/env python
  import json
  import os
  import pandas as pd

  bamfiles = [os.path.basename(b) for b in "$sortedBams".split(' ')]
  jsDict = {"base": "${absOutDir}/ext_hairpins",
            "mir.anno": "$mirArmAnno"}

  readCounts = pd.read_table("$readcounts")

  # with open('countsum.txt', 'r') as counts:
  #   countStats = counts.readline()

  # fullName, mapped, unmapped = countStats.split('\t')

  samples = []

  # print('Debug-STR :: {}'.format(bamfiles))
  # print('Debug-Counts :: {} - {} - {}'.format(fullName, mapped, unmapped))

  for bam in bamfiles:
    nameElements = bam.split('_')
    idx = int(nameElements[0])
    # sRNAreads = 1000000
    sRNAreads = readCounts[readCounts.idx == idx]['sRNAreads'].values[0]
    time = readCounts[readCounts.idx == idx]['time'].values[0]
    miRNAreads = readCounts[readCounts.idx == idx]['miRNAreads'].values[0]
    samples.append({"id": idx,
                    "time": time,
                    "align": bam,
                    "sRNAreads": sRNAreads,
                    "miRNAreads": miRNAreads})

  jsDict["samples"] = samples

  with open('samples.json', 'w') as fp:
    json.dump(jsDict, fp)
  """
}

/*
 *  STEP 5b - Postprocessing -- Getting alignment statistics for miR hairpin mappers
 *                              1) Count all reads and T>C reads with T>C BQ>27
 *                              2) Get average reads across all samples
 *                              3) Summarise reads per length isoform, and total Reads
 *                              4) Save unfiltered information to `allCounts.tsv`
 *                              5) Keep wide (untidy) tbl with abundance of length isoforms
 *                                 as well as total lengths
 *                              6) Only keep (top 5 most abundant) isoforms that are represented
 *                                 across all datasets
 *                                 `filteredPositions.tsv` -> steady state
 *                                 `filteredLenDis.tsv`    -> length isoforms
 *                              7) Add metadata (seed, UCount, arm type, arm name)
 *  Currently NOT done (anymore):
 *                              -) Any type of filtering
 *                              -) Processing data as 'ready to plot'
 */

process alignmentStats {
  publishDir path: getOutDir('stats'), mode: "copy", pattern: "*.tsv"

  input:
  file readCountConfig
  file hairpinFasta

  output:
  file 'allCounts.tsv' into allCounts
  file 'filteredPositions.tsv' into topPositions
  file 'filteredLenDis.tsv' into lendisTSV

  script:
  """
  export OMP_NUM_THREADS=${task.cpus}
  getAllCounts.R $baseDir ${params.rlocation} $readCountConfig $hairpinFasta
  """
}

/*
 *  STEP 5c - Postprocessing -- Process BAM files in parallel using BiocParallel
 *                              -) Get mutations for all positions in a miR, for all miRNAs
 *                                 provided in `filteredPositions.tsv`
 */
process mutationStats {
  publishDir path: getOutDir('stats'), mode: "copy", pattern: "*.tsv"

  input:
  file readCountConfig
  file topPositions
  file hairpinFasta

  output:
  file 'miRs.wAllMuts.tsv' into mutStats
  file 'mirMutCodes.tsv' into mirMuts
  file 'mirMutsWide.tsv' into wideMuts

  script:
  """
  export OMP_NUM_THREADS=${task.cpus}
  export TMPDIR=/scratch-ii2/users/reichholf/tmp
  export TMP=/scratch-ii2/users/reichholf/tmp
  export TEMP=/scratch-ii2/users/reichholf/tmp
  getMutStats.R $baseDir ${params.rlocation} ${task.cpus} $readCountConfig $topPositions $hairpinFasta
  """
}
