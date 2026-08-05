# Made by Reese Saho for running ADONIS and ANCOM-BC on the BSF meta-analysis samples. 2025/08/25

# Run adonis testing to determine variance attributed to study_id and diet_condensed


# qiime feature-table summarize \
#   --i-table $deblur/bsf-merged-p_0.01-filtered-table.qza \
#   --o-visualization $deblur/bsf-merged-p_0.01-filtered-table.qzv \
#   --m-sample-metadata-file $metadata/bsf_metadata.txt

# qiime metadata filter-samples \
#   --m-metadata-file $metadata/bsf_metadata.txt \
#   --p-where "diet_condensed IS NOT NULL" \
#   --o-filtered-metadata $metadata/md_no_na.qza


# for matrix in $base_dir/core-metrics-phylogeny-10k/unweighted_unifrac_distance_matrix.qza $base_dir/core-metrics-phylogeny-10k/weighted_unifrac_distance_matrix.qza $base_dir/core-metrics-phylogeny-10k/unweighted_unifrac_distance_matrix.qza $base_dir/core-metrics-phylogeny-10k/weighted_unifrac_distance_matrix.qza; do

# #  $base_dir/core-metrics-phylogeny-10k/generalized_unifrac.qza

#     if [ $matrix = "$base_dir/core-metrics-phylogeny-10k/weighted_unifrac_distance_matrix.qza" ] ; then 

#         name="weighted"
        
#         echo "Running ADONIS testing on the $name distance matrix."

#         qiime diversity adonis \
#             --i-distance-matrix ./core-metrics-phylogeny-10k/unweighted_unifrac_distance_matrix.qza \
#             --p-formula "diet_condensed" \
#             --p-n-jobs 4 \
#             --m-metadata-file ./metadata/md_no_na.qza \
#             --o-visualization ./core-metrics-phylogeny-10k/${name}-dc-adonis.qzv

#         qiime diversity adonis \
#             --i-distance-matrix $matrix \
#             --p-formula "diet_condensed" \
#             --p-n-jobs 4 \
#             --m-metadata-file $metadata/md_no_na.qza \
#             --o-visualization $base_dir/core-metrics-phylogeny-10k/${name}-dc-adonis.qzv


#     elif [ $matrix = "$base_dir/core-metrics-phylogeny-10k/unweighted_unifrac_distance_matrix.qza" ] ; then

#         name="unweighted"

#         echo "Running ADONIS testing on the $name distance matrix."

#         qiime diversity adonis \
#             --i-distance-matrix $matrix \
#             --p-formula "diet_condensed" \
#             --p-n-jobs 4 \
#             --m-metadata-file $metadata/md_no_na.qza \
#             --o-visualization $base_dir/core-metrics-phylogeny-10k/${name}-dc-adonis.qzv

#     elif [ $matrix = "$correction/core-metrics-current_study/unweighted_unifrac_distance_matrix.qza" ] ; then
    
#         name="unweighted"

#         echo "Running ADONIS testing on the $name distance matrix."

#         qiime diversity adonis \
#             --i-distance-matrix $matrix \
#             --p-formula "diet_condensed" \
#             --p-n-jobs 4 \
#             --m-metadata-file $metadata/md_no_na.qza \
#             --o-visualization $correction/core-metrics-current_study/${name}-dc-adonis.qzv

#     else

#         name="weighted"

#         echo "Running ADONIS testing on the $name distance matrix."

#         qiime diversity adonis \
#             --i-distance-matrix $matrix \
#             --p-formula "diet_condensed" \
#             --p-n-jobs 4 \
#             --m-metadata-file $metadata/md_no_na.qza \
#             --o-visualization $correction/core-metrics-current_study/${name}-dc-adonis.qzv
#     fi

# echo "Running ADONIS testing on the $name distance matrix."


#     qiime diversity adonis \
#         --i-distance-matrix $matrix \
#         --p-formula "diet_condensed" \
#         --p-n-jobs 4 \
#         --m-metadata-file $metadata/bsf_metadata.txt \
#         --o-visualization $base_dir/core-metrics-phylogeny-10k/${name}-dc-adonis.qzv

    # qiime tools export \
    #       --input-path $matrix \
    #       --output-path $base_dir/core-metrics-phylogeny-10k/unweighted_unifrac_distance_matrix

# echo "ADONIS testing complete."

# done

# Convert ConQuR tables to BIOM format and import them into QIIME2


# for table in ${correction}/correction/taxa_table_penalized_corrected_current_study.tsv; do

#     if [ $table = "${correction}/correction/taxa_table_penalized_corrected_bruno2019.tsv" ] ; then
    
#         name="bruno2019"
        
#     else
    
#         name="current_study"
        
#     fi

#     # qiime feature-table summarize \
#     #   --i-table $correction/core-metrics-${name}/rarefied_table.qza \
#     #   --o-visualization $correction/core-metrics-${name}/rarefied_table.qzv \
#     #   --m-sample-metadata-file $metadata/bsf_metadata.txt


#     # biom convert \
#     #     -i $table \
#     #     -o ${name}-cor-feature_table.biom \
#     #     --table-type="OTU table" \
#     #     --to-hdf5
    
#     # qiime tools import \
#     #     --input-path ${name}-cor-feature_table.biom \
#     #     --type 'FeatureTable[Frequency]' \
#     #     --input-format BIOMV210Format \
#     #     --output-path ${name}-cor-feature_table.qza
    
#     # qiime feature-table summarize \
#     #     --i-table ${name}-cor-feature_table.qza \
#     #     --o-visualization ${name}-cor-feature_table.qzv \
#     #     --m-sample-metadata-file $metadata/bsf_metadata.txt
    
#     # # Remove features that are not found in any samples
    
#     # qiime feature-table filter-features \
#     #     --i-table ${name}-cor-feature_table.qza \
#     #     --p-min-samples 1 \
#     #     --o-filtered-table ${correction}/correction/${name}-cor-feature_table-filt.qza
    
    
#     # qiime feature-table summarize \
#     #     --i-table ${correction}/correction/${name}-cor-feature_table-filt.qza \
#     #     --m-sample-metadata-file $metadata/bsf_metadata.txt \
#     #     --o-visualization ${correction}/correction/${name}-cor-feature_table-filt.qzv
    
#     # # Filter the phylogenetic tree to match the feature table
    
#     # qiime phylogeny filter-tree \
#     #     --i-tree $phylogenetics/bsf-tree.qza \
#     #     --i-table ${correction}/correction/${name}-cor-feature_table-filt.qza \
#     #     --o-filtered-tree $phylogenetics/bsf-tree-${name}-cor.qza
        
#     # # Run core metrics
    
#     # qiime diversity core-metrics-phylogenetic \
#     #     --i-table ${correction}/correction/${name}-cor-feature_table-filt.qza \
#     #     --i-phylogeny $phylogenetics/bsf-tree-${name}-cor.qza \
#     #     --p-sampling-depth 10000 \
#     #     --p-n-jobs-or-threads 0 \
#     #     --m-metadata-file $metadata/bsf_metadata.txt \
#     #     --output-dir $correction/core-metrics-${name}/

#     # Run a core features analysis
    
#     # qiime feature-table core-features \
#     #     --i-table ${correction}/correction/${name}-cor-feature_table-filt.qza \
#     #     --p-min-fraction 0.5 \
#     #     --p-max-fraction 1.0 \
#     #     --o-visualization $correction/core-metrics-${name}/${name}-core_features.qzv
        
        
#     # for matrix in bray_curtis jaccard unweighted_unifrac weighted_unifrac; do

#     #     qiime diversity adonis \
#     #         --i-distance-matrix $correction/core-metrics-${name}/${matrix}_distance_matrix.qza \
#     #         --p-formula "study_id" \
#     #         --p-n-jobs 3 \
#     #         --m-metadata-file $metadata/bsf_metadata.txt \
#     #         --o-visualization $correction/core-metrics-${name}/${name}-${matrix}-dc+si-adonis.qzv
            
            
#     #     # qiime diversity beta-group-significance \
#     #     #     --i-distance-matrix $correction/core-metrics-${name}/${matrix}_distance_matrix.qza \
#     #     #     --m-metadata-file $metadata/bsf_metadata.txt \
#     #     #     --m-metadata-column "diet_condensed" \
#     #     #     --p-method permdisp \
#     #     #     --o-visualization $correction/core-metrics-${name}/${name}-${matrix}-permdisp.qzv


#     # done

#     # qiime tools export \
#     #   --input-path $correction/core-metrics-${name}/unweighted_unifrac_distance_matrix.qza \
#     #   --output-path $correction/core-metrics-${name}/exported_distance

# done        

# Run the Mantel test comparing the uncorrected vs. corrected UniFracs

qiime diversity mantel \
    --i-dm1 ${base_dir}/core-metrics-phylogeny-10k/unweighted_unifrac_distance_matrix.qza \
    --i-dm2 ${correction}/core-metrics-current_study/unweighted_unifrac_distance_matrix.qza \
    --p-intersect-ids \
    --p-label1 "Uncorrected Unweighted UniFrac" \
    --p-label2 "Corrected Unweighted UniFrac" \
    --o-visualization ${correction}/core-metrics-current_study/unweighted_mantel.qzv
    
qiime diversity mantel \
    --i-dm1 ${base_dir}/core-metrics-phylogeny-10k/weighted_unifrac_distance_matrix.qza \
    --i-dm2 ${correction}/core-metrics-current_study/weighted_unifrac_distance_matrix.qza \
    --p-intersect-ids \
    --p-label1 "Uncorrected Weighted UniFrac" \
    --p-label2 "Corrected Weighted UniFrac" \
    --o-visualization ${correction}/core-metrics-current_study/weighted_mantel.qzv
    
    
    
    
    