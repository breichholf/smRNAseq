setupRLibs <- function(){
  if(!require(devtools)) install.packages('devtools')
  if(!require(pacman)) install.packages('pacman')

  if(!require(knitr) | !require(tidyverse) | !require(stringr) | !require(optparse)) {
    pacman::p_install(knitr, tidyverse, stringr, optparse)
  }
  if(!require(cowplot)) pacman::p_install_gh(c('wilkelab/cowplot'))
}

prepareFlybaseGFF <- function(gff, conversion, output) {
  require(tidyverse)
  require(stringr)

  # Read in file
  fbGFF <- read_tsv(opt$gff,
                    col_names = c('chr', 'source', 'type', 'start', 'end', 'score', 'strand', 'dot', 'tags'),
                    col_types = list('c', 'c', 'c', 'n', 'n', 'c', 'c', 'c', 'c'))
  fbRows <- fbGFF %>% separate_rows(tags, sep = ";")
  fbIDs <-
    fbRows %>% filter(grepl("ID=FBtr", tags)) %>%
    mutate(id = str_sub(tags, 4, length(tags))) %>%
    select(type, id) %>% rename(parent.type = type)
  
  # miRNAs have "Parent=FBgn", so they're filtered out here
  fbGenicParents <-
    fbRows %>% filter(grepl("Parent=FBtr", tags)) %>%
    mutate(id = str_sub(tags, 8, length(tags))) %>%
    rename(source.type = type) %>% left_join(fbIDs, by = "id")
  
  fbGenic <-
    fbGenicParents %>% filter(!grepl("mito", chr) & grepl("exon|intron|UTR", source.type)) %>%
    mutate(print.start = ifelse(parent.type == "tRNA", ifelse(start < 20, 0, start - 20),
                                start),
           print.end = ifelse(parent.type == "tRNA",
                              end + 20,
                              end),
           print.name = ifelse(grepl("mRNA|pseudogene|ncRNA", parent.type),
                               paste(parent.type, source.type, sep = "_"),
                               parent.type),
           score = 0) %>%
    select(chr, print.start, print.end, print.name, score, strand)
  
  fbMirs <-
    fbGFF %>% filter(!grepl("mito", chr) & grepl("miRNA", type)) %>%
    mutate(score = 0) %>%
    rename(print.start = start, print.end = end, print.name = type) %>%
    select(chr, print.start, print.end, print.name, score, strand)
  
  # Concatenate genic and miRs
  concatBed <- rbind(fbGenic, fbMirs)
  
  # Convert Chr Names
  # Read in file first - first lines contain "#" which are comments.
  assemblyReport <- read_tsv(opt$nameconversion, comment = "#",
                             col_names = c('fb.name', 'sequence.role', 'assigned.molecule',
                                           'molecule.loaction.type', 'genbank.accn', 'relationship',
                                           'RefSeq.accn', 'Assembly.Unit', 'sequence.length', 'ucsc.name'),
                             col_types = list('c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c', 'c'))
  
  # Select appropriate columns and remove leading "chr" from UCSC Names
  chrConversion <-
    assemblyReport %>%
    mutate(new.ucsc = str_replace(ucsc.name, "^chr", "")) %>%
    select(fb.name, new.ucsc) %>%
    mutate(fb.name = ifelse(new.ucsc == "M", "mitochondrion_genome", fb.name))
  
  # Convert Flybase Names to UCSC-compatible names
  outputBed <-
    concatBed %>% left_join(chrConversion, by = c('chr' = 'fb.name')) %>%
    select(new.ucsc, print.start, print.end, print.name, score, strand)
  
  # Write output file, will overwrite old files
  outputBed %>% write_tsv(opt$output, col_names = FALSE)
}
