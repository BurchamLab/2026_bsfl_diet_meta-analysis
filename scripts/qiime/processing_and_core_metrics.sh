#This code was written by Reese Saho for handling of the Deblur denoising through the generation of core metrics for the BSFL meta-analysis. 2025/07/20


# # Import the samples from the current study



# for study in liu2020b querejeta2022 li2022 silvaraju2024 tanga2021 wu2020 kluber2022 klammsteiner2025 edwards2024 zhang2020 galassi2021 wu2021 chen2023 ; do

# # current_study klammsteiner2020 wynants2018 vandeweyer2023 klammsteiner2021 gorrens2022 deng2025 reyer2025 li2025  dragone2025 defilippis2023 pei2022 li2021 jiang2019 greenwood2021 shelomi2020 xu2025 cifuentes2020 piersanti2024 bruno2019 marasco2022 yu2023 tegtmeier2021a tegtmeier2021b liu2021 liu2020

#     if [ $study = "current_study" ] || [ $study = "klammsteiner2020" ] || [ $study = "wynants2018" ] || [ $study = "vandeweyer2023" ] || [ $study = "klammsteiner2021" ] || [ $study = "gorrens2022" ] || [ $study = "reyer2025" ] || [ $study = "dragone2025" ] || [ $study = "jiang2019" ] || [ $study = "liu2020" ] || [ $study = "klammsteiner2025" ] || [ $study = "liu2020b" ] ; then

#     # Look for the 515 primer sequence, trim everything before it
#         qiime cutadapt trim-single \
#         --i-demultiplexed-sequences $demux/${study}-single-demux.qza \
#         --p-cores 0 \
#         --p-front GTGYCAGCMGCCGCGGTAA \
#         --o-trimmed-sequences $demux/${study}-single-515_trim-demux.qza

#     else

#         qiime cutadapt trim-single \
#         --i-demultiplexed-sequences $demux/${study}-single-demux.qza \
#         --p-cores 0 \
#         --p-front GTGYCAGCMGCCGCGGTAA \
#         --p-discard-untrimmed \
#         --o-trimmed-sequences $demux/${study}-single-515_trim-demux.qza

#     fi

#     qiime demux summarize \
#     --i-data $demux/${study}-single-515_trim-demux.qza \
#     --o-visualization $demux/${study}-single-515_trim-demux.qzv

#     ### Denoise with deblur, trimmed to 250 to generate 16S V4

#     echo "Starting Deblur denoising on samples from $study."

#     qiime deblur denoise-16S \
#         --i-demultiplexed-seqs $demux/${study}-single-515_trim-demux.qza \
#         --p-trim-length 250 \
#         --o-representative-sequences $deblur/${study}-repseqs.qza \
#         --o-table $deblur/${study}-feature_table.qza \
#         --o-stats $deblur/${study}-deblur_stats.qza \
#         --p-sample-stats \
#         --p-jobs-to-start 0

#     ## Visualize deblur stats

#     qiime deblur visualize-stats \
#         --i-deblur-stats $deblur/${study}-deblur_stats.qza \
#         --o-visualization $deblur/${study}-deblur_stats.qzv

#     ## Summarize deblur table

#     qiime feature-table summarize \
#         --i-table $deblur/${study}-feature_table.qza \
#         --o-visualization $deblur/${study}-feature_table.qzv

#     ## Visualize deblur rep seqs

#     qiime feature-table tabulate-seqs \
#         --i-data $deblur/${study}-repseqs.qza \
#         --o-visualization $deblur/${study}-repseqs.qzv

# done

# # Merge the feature tables and repseqs together from all studies

# qiime feature-table merge \
#     --i-tables $deblur/*-feature_table.qza \
#     --o-merged-table $deblur/merged-feature-table.qza

# qiime feature-table merge-seqs \
#     --i-data $deblur/*-repseqs.qza \
#     --o-merged-data $deblur/merged-repseqs.qza

# # Download GreenGenes 2 Classifier

# wget -O $phylogenetics/gg_classifier.qza https://data.qiime2.org/classifiers/sklearn-1.4.2/greengenes2/2024.09.backbone.v4.nb.sklearn-1.4.2.qza

# #Build taxonomy using SILVA-138 classifier

# qiime feature-classifier classify-sklearn \
#     --i-classifier $phylogenetics/gg_classifier.qza \
#     --i-reads $deblur/merged-repseqs.qza \
#     --o-classification $taxonomy/taxonomy.qza

# # Export the taxonomy for later use

# qiime tools export \
#     --input-path $taxonomy/taxonomy.qza \
#     --output-path $taxonomy/taxa/

# #Filtering out mitochondria & chloroplasts from table

# qiime taxa filter-table \
#     --i-table $deblur/merged-feature-table.qza \
#     --i-taxonomy $taxonomy/taxonomy.qza \
#     --p-exclude mitochondria,chloroplast \
#     --o-filtered-table $deblur/bsf-table-no-mc.qza

# # Filter out blanks

# qiime feature-table filter-samples \
#     --i-table $deblur/bsf-table-no-mc.qza \
#     --m-metadata-file $metadata/bsf_metadata.txt \
#     --p-where "[control]='yes'" \
#     --p-exclude-ids \
#     --o-filtered-table $deblur/bsf-table-no-mc_controls.qza

# # Filter out probiotic manure samples

# qiime feature-table filter-samples \
#     --i-table $deblur/bsf-table-no-mc_controls.qza \
#     --m-metadata-file $metadata/bsf_metadata.txt \
#     --p-where "[diet]='cow_manure_probiotic'" \
#     --p-exclude-ids \
#     --o-filtered-table $deblur/bsf-table-no-mc_controls_manpro.qza

# # Filter out probiotic manure samples

# qiime feature-table filter-samples \
#     --i-table $deblur/bsf-table-no-mc_controls_manpro.qza \
#     --m-metadata-file $metadata/bsf_metadata.txt \
#     --p-where "[diet]='initial_larvae'" \
#     --p-exclude-ids \
#     --o-filtered-table $deblur/bsf-table-no-mc_controls_manpro_il.qza

# qiime feature-table summarize \
#     --i-table $deblur/bsf-table-no-mc_controls_manpro_il.qza \
#     --o-visualization $deblur/bsf-table-no-mc_controls_manpro_il.qzv

# # #Filtering out features that appear less than 3 samples from the table

#     qiime feature-table filter-features \
#         --i-table $deblur/bsf-table-no-mc_controls_manpro_il.qza \
#         --p-min-samples 3 \
#         --o-filtered-table $deblur/bsf-merged-3p-filtered-table.qza

#     #Visualizing the filtered table

#     qiime feature-table summarize \
#         --i-table $deblur/bsf-merged-3p-filtered-table.qza \
#         --o-visualization $deblur/bsf-merged-3p-filtered-table.qzv


# #Filtering out features that appear less than 3 samples from the rep seqs

# qiime feature-table filter-seqs \
#     --i-data $deblur/merged-repseqs.qza \
#     --i-table$deblur/bsf-merged-3p-filtered-table.qza \
#     --o-filtered-data $deblur/bsf-merged-filtered-seqs-no-mc.qza

# #Visualizing the filtered rep seqs
# qiime feature-table tabulate-seqs \
#     --i-data $deblur/bsf-merged-filtered-seqs-no-mc.qza \
#     --o-visualization $deblur/bsf-merged-filtered-seqs-no-mc.qzv

# #Generating tree

# wget -O $phylogenetics/sepp-refs-silva-128.qza "https://data.qiime2.org/2021.8/common/sepp-refs-silva-128.qza"

# qiime fragment-insertion sepp \
#     --i-representative-sequences $deblur/bsf-merged-filtered-seqs-no-mc.qza \
#     --i-reference-database $phylogenetics/sepp-refs-silva-128.qza \
#     --p-threads 10 \
#     --o-tree $phylogenetics/bsf-tree.qza \
#     --o-placements $phylogenetics/bsf-tree-placements.qza

# # Generation of Alpha rarefaction

# qiime diversity alpha-rarefaction \
#     --i-table $deblur/bsf-merged-p_0.01-filtered-table.qza \
#     --i-phylogeny $phylogenetics/bsf-tree.qza \
#     --p-max-depth 20000 \
#     --p-steps 20 \
#     --m-metadata-file $metadata/bsf_metadata.txt \
#     --o-visualization $rarefaction/"bsf-alpha-rarefaction-20k.qzv"

# ##Creating the Core Metrics Files

# qiime feature-table filter-samples \
#     --i-table $deblur/bsf-merged-p_0.01-filtered-table.qza \
#     --m-metadata-file $metadata/bsf_metadata.txt \
#     --o-filtered-table $deblur/bsf-metadata-feature_table.qza

# qiime feature-table filter-samples \
#     --i-table $deblur/bsf-metadata-feature_table.qza \
#     --m-metadata-file $metadata/bsf_metadata.txt \
#     --p-where "[diet]='initial_larvae'" \
#     --p-exclude-ids \
#     --o-filtered-table $deblur/bsf-feature_table-final.qza


# qiime diversity core-metrics-phylogenetic \
#     --i-table $deblur/bsf-feature_table-final.qza \
#     --i-phylogeny $phylogenetics/bsf-tree.qza \
#     --p-sampling-depth 10000 \
#     --p-n-jobs-or-threads 4 \
#     --m-metadata-file $metadata/bsf_metadata.txt \
#     --output-dir $base_dir/core-metrics-phylogeny-10k/

# #Generating Generalized UniFrac distances on the unfiltered table

# qiime diversity beta-phylogenetic \
#     --i-table $deblur/bsf-feature_table-final.qza \
#     --i-phylogeny $phylogenetics/bsf-tree.qza \
#     --p-alpha 0.5 \
#     --p-metric generalized_unifrac \
#     --p-threads 4 \
#     --o-distance-matrix $base_dir/core-metrics-phylogeny-10k/generalized_unifrac.qza

# #Generating PCOA for the Unfiltered Generalized UniFrac

# qiime diversity pcoa \
#     --i-distance-matrix $base_dir/core-metrics-phylogeny-10k/generalized_unifrac.qza \
#     --o-pcoa $base_dir/core-metrics-phylogeny-10k/generalized_unifrac_pcoa.qza

# #Visualizing the Unfiltered Generalized PCOA

# qiime emperor plot \
#     --i-pcoa $base_dir/core-metrics-phylogeny-10k/generalized_unifrac_pcoa.qza \
#     --m-metadata-file $metadata/bsf_metadata.txt \
#     --o-visualization $base_dir/core-metrics-phylogeny-10k/generalized_unifrac_emperor.qzv
