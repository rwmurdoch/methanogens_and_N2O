#!/bin/bash

qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path reads \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path reads.qza

qiime demux summarize \
--i-data reads.qza \
--o-visualization demux.qzv \
--verbose

## trim the primers (341F/785R)
## stdout is directed to a logfile

qiime cutadapt trim-paired \
--i-demultiplexed-sequences reads.qza \
--p-cores 8 \
--p-front-f CCTACGGGNGGCWGCAG \
--p-front-r GACTACHVGGGTATCTAATCC \
--o-trimmed-sequences data-cutadapt \
--verbose > cutadapt_log.txt

###############################################
# training the feature classifier
###############################################

## working with full-length silva138 database (seqs and taxonomy)
## downloaded from: https://docs.qiime2.org/2021.2/data-resources/?highlight=silva
## Database is first trimmed using the same primers used for the study (341F/785R)

qiime feature-classifier extract-reads \
  --i-sequences silva-138-99-seqs.qza \
  --p-f-primer CCTACGGGNGGCWGCAG \
  --p-r-primer GACTACHVGGGTATCTAATCC \
  --o-reads ref-seqs.qza \
  --verbose

qiime feature-classifier fit-classifier-naive-bayes \
--i-reference-reads ref-seqs.qza \
--i-reference-taxonomy silva-138-99-tax.qza \
--o-classifier classifier.qza


#####################################
## dada2 and classification
####################################33333

rm -R data2-stats
rm -R dada2-stats

qiime dada2 denoise-paired \
--i-demultiplexed-seqs data-cutadapt.qza \
--o-table dada2-table \
--o-representative-sequences dada2_seqs \
--o-denoising-stats dada2-stats \
--p-n-threads 16 \
--p-trunc-len-f 0 \
--p-trunc-len-r 0
#--p-min-fold-parent-over-abundance 4 # I used this in a previous script in which I erroneausly fed seqs without primer removal; check results this time around

qiime feature-table summarize \
--i-table dada2-table.qza \
--o-visualization dada2-seq-stats \
--verbose

qiime feature-classifier classify-sklearn \
  --i-classifier classifier.qza \
  --i-reads dada2_seqs.qza \
  --p-n-jobs 16 \
  --o-classification dada2-taxonomy

qiime metadata tabulate \
  --m-input-file dada2-taxonomy.qza \
  --o-visualization dada2-taxonomy

qiime taxa barplot \
  --i-table dada2-table.qza \
  --i-taxonomy dada2-taxonomy.qza \
  --m-metadata-file metadata.csv \
  --o-visualization dada2-taxa-bar-plots

qiime feature-table summarize \
  --i-table dada2-table.qza \
  --o-visualization dada2-table.qzv \
  --m-sample-metadata-file metadata.csv

qiime feature-table tabulate-seqs \
  --i-data dada2_seqs.qza \
  --o-visualization dada2_seqs.qzv

#######################################
# preparation for alpha div metrics
#######################################

qiime alignment mafft \
  --i-sequences dada2_seqs.qza \
  --o-alignment dada2_seqs-aligned.qza \
  --p-n-threads 16

qiime alignment mask \
 --i-alignment dada2_seqs-aligned.qza \
--o-masked-alignment dada2_seqs-aligned-masked.qza

qiime phylogeny fasttree \
 --i-alignment dada2_seqs-aligned-masked.qza \
 --o-tree unrooted-tree-dada2.qza

qiime phylogeny midpoint-root \
--i-tree unrooted-tree-dada2.qza \
--o-rooted-tree rooted-tree-dada2.qza

################################################
# generate core diversity metrics #
# pay close attention to the sampling depth command #
###################################################

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree-dada2.qza \
  --i-table dada2-table.qza \
  --p-sampling-depth 30000 \
  --m-metadata-file metadata.csv \
  --output-dir core-metrics-results-dada2

## visualization of diversity metrics ##

qiime diversity alpha-rarefaction \
  --i-table dada2-table.qza \
  --i-phylogeny rooted-tree-dada2.qza \
  --p-min-depth 1 \
  --p-max-depth 30000 \
  --m-metadata-file metadata.csv \
  --o-visualization open-alpha-rarefaction-dada2.qzv

###############################################
# exporting data
# you can easily export seqs and taxonomy directly out of
# their corresponding .qza files using the qiime tools export feature.
# to get the feature table, you export the table in biom format and then convert
# into a simple text file
#################################

qiime tools export --input-path dada2-table.qza --output-path exports
qiime tools export --input-path dada2_seqs.qza --output-path exports
qiime tools export --input-path dada2-taxonomy.qza --output-path exports

biom convert -i exports/feature-table.biom \
-o exports/otu_table.txt --to-tsv

## combine everything into a single table
## this uses the R script create.OTU.table.R
## script is available here: https://github.com/rwmurdoch/project.scripts/blob/master/create.OTU.table.R
## download/copy the script into your project directory; it will read and write in the "export" directory

Rscript create.OTU.table.R

## combine various outputs into a results package

tar -cf results.tar.gz \
exports \
core-metrics-results-dada2 \
open-alpha-rarefaction-dada2.qzv \
dada2-taxa-bar-plots.qzv \
metadata.csv \
dada2-seq-stats.qzv \
demux.qzv
