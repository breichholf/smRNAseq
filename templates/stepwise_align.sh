#!/usr/bin/bash

# Argv[1]: virus index (with yes/no built in)
# Argv[2]: ribo index
# Argv[3]: tRNA index
# Argv[4]: snRNA index
# Argv[5]: snoRNA index
# Argv[6]: reads
# Argv[7]: prefix
# Argv[8]: cpus
# Argv[9]: mismatches

vAnno
vIDX
rIDX
tIDX
snIDX
snoIDX
reads
prefix
mismatches
cpus

if [[ -f TCtag_virus.sam ]]; then
    samples="virus ribo tRNA snRNA snoRNA"
else
    samples="ribo tRNA snRNA snoRNA"
fi

samtools faidx $vAnno
bowtie-build $vAnno $vIDX

alignCmd="bowtie -a --best --strata -v $mismatches -S -p $cpus"

$alignCmd $vIDX -q $reads --un vClean.fq | arrayTagTCreads.awk > TCtag_virus.sam

$alignCmd $riboIDX -q vClean.fq --un riboClean.fq | arrayTagTCreads.awk > TCtag_ribo.sam

$alignCmd $tIDX -q riboClean.fq --un tClean.fq | arrayTagTCreads.awk > TCtag_tRNA.sam

$alignCmd $snIDX -q tClean.fq --un snClean.fq | arrayTagTCreads.awk > TCtag_snRNA.sam

$alignCmd $snoIDX -q snClean.fq --un unalClean.fq | arrayTagTCreads.awk > TCtag_snoRNA.sam

for sam in \$samples; do
    samtools view -H TCtag_\${sam}.sam | awk '{
        if(\$1 ~ /^@SQ/) { print \$0 > \${sam}_sq.txt }
        if(\$1 ~ /^@PG/) { print \$0" sample: ${prefix}"  > \${sam}_pg.txt }
      }'
done

cat virus_sq.txt ribo_sq.txt tRNA_sq.txt snRNA_sq.txt snoRNA_sq.txt > ncRNA_sq.txt
cat virus_pg.txt ribo_pg.txt tRNA_pg.txt snRNA_pg.txt snoRNA_pg.txt > ncRNA_pg.txt

for sample in \$samples; do

done

