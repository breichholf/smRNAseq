#!/usr/bin/env Rscript

source('functions.R')

suppressPackageStartupMessages(library("optparse"))

# Provide command line options -g/--gff, -n/--names and -o/--output
options <- list(
  make_option(c('-g', '--gff'), dest = 'gff', required = TRUE, help = 'dmel-all-VERSION.gff - download from ftp://ftp.flybase.org/genomes/dmel/current/gff/'),
  make_option(c('-n', '--names'), dest = 'nameconversion', required = TRUE,
              help = 'UCSC<->FlyBase name conversion file: ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/invertebrate/Drosophila_melanogaster/latest_assembly_versions/GCA_000001215.4_Release_6_plus_ISO1_MT/GCA_000001215.4_Release_6_plus_ISO1_MT_assembly_report.txt'),
  make_option(c('-o', '--output'), dest = 'output', required = TRUE, help = '6-column BED output file')
)

opt <- parse_args(OptionParser(option_list = options))

if (opt$gff & opt$nameconversion & opt$output) {
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(stringr))
  
  # Todo: Check if GFF exists
  #       Check if nameconversion exists -> download otherwise
  
  prepareFlybaseGFF(opt$gff, opt$nameconversion, opt$output)
}