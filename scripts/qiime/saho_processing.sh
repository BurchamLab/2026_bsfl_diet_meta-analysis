
base_dir="/lustre/isaac24/proj/UTK0406/ReeseSaho/projects/bsf_metaa"
raw_seqs="$base_dir/raw_sequencing_files"
demux="$base_dir/demux"
deblur="$base_dir/deblur"
scripts="$base_dir/scripts"
logs="$base_dir/logs"
phylogenetics="$base_dir/phylogenetics"
taxonomy="$base_dir/taxonomy"
metadata="$base_dir/metadata"
ancom="$base_dir/ancom"
figures="$base_dir/figures"

mkdir -p "$phylogenetics"
mkdir -p "$taxonomy"
mkdir -p "$metadata"
mkdir -p "$ancom"
mkdir -p "$figures"

# # Filtering the no mitochondria/chloroplast feature table to leave only current_study samples

# qiime feature-table filter-samples \
#     --i-table $deblur/bsf-table-no-mc.qza \
#     --m-metadata-file $metadata/bsf_metadata.txt \
#     --p-where "[study_id]='current_study'" \
#     --o-filtered-table $deblur/current_study-table-no-mc.qza

# qiime feature-table summarize \
#     --i-table $deblur/current_study-table-no-mc.qza \
#     --o-visualization $deblur/current_study-table-no-mc.qzv


# #Removing the blanks and controls from the feature table

# qiime feature-table filter-samples \
#     --i-table $deblur/current_study-table-no-mc.qza \
#     --m-metadata-file $metadata/bsf_metadata.txt \
#     --p-where "[diet]='na'" \
#     --p-exclude-ids \
#     --o-filtered-table $deblur/current_study-table-no-mc-controls.qza

# qiime feature-table summarize \
#     --i-table $deblur/current_study-table-no-mc-controls.qza \
#     --o-visualization $deblur/current_study-table-no-mc-controls.qzv


# #Filtering out features that appear less than 1% of all samples from the table

# echo "Filtering out features with >1% prevalence and chloroplasts from the feature table and rep seqs."


# qiime feature-table filter-features \
#     --i-table $deblur/current_study-table-no-mc-controls.qza \
#     --p-min-samples $i \
#     --o-filtered-table $deblur/current_study-${i}p-filtered-table.qza

# #Visualizing the filtered table

# qiime feature-table summarize \
#     --i-table $deblur/current_study-${i}p-filtered-table.qza \
#     --o-visualization $deblur/current_study-${i}p-filtered-table.qzv


# #Filtering the rep seqs to match the feature table - with pro

# qiime feature-table filter-seqs \
#     --i-data $deblur/merged-repseqs.qza \
#     --i-table $deblur/current_study-p-filtered-table.qza \
#     --o-filtered-data $deblur/current_study-p-filtered-seqs.qza

# #Visualizing the filtered rep seqs

# qiime feature-table tabulate-seqs \
#     --i-data $deblur/current_study-p-filtered-seqs.qza \
#     --o-visualization $deblur/current_study-p-filtered-seqs.qzv

# # Removing the manure samples with probiotics

# echo "Removing manure + probiotic samples."

# qiime feature-table filter-samples \
#     --i-table $deblur/current_study-p-filtered-table.qza \
#     --m-metadata-file $metadata/bsf_metadata.txt \
#     --p-where "[diet]='cow_manure_probiotic'" \
#     --p-exclude-ids \
#     --o-filtered-table $deblur/current_study-table-no-probiotic.qza

# qiime feature-table summarize \
#     --i-table $deblur/current_study-table-no-probiotic.qza \
#     --o-visualization $deblur/current_study-table-no-probiotic.qzv

# qiime feature-table filter-samples \
#     --i-table $deblur/current_study-p-filtered-table.qza \
#     --p-min-samples 3
#     --o-filtered-table $deblur/current_study-table-no-probiotic_3p.qza

# qiime feature-table summarize \
#     --i-table $deblur/current_study-table-no-probiotic_3p.qza \
#     --o-visualization $deblur/current_study-table-no-probiotic_3p.qzv


 
# # Filtering the rep seqs to remove the probiotic ASVs

# qiime feature-table filter-seqs \
#     --i-data $deblur/current_study-table-no-probiotic_3p.qza \
#     --i-table $deblur/current_study-table-no-probiotic.qza \
#     --o-filtered-data $deblur/current_study-table-no-probiotic-seqs.qza

#     #Visualizing the filtered rep seqs
# qiime feature-table tabulate-seqs \
#     --i-data $deblur/current_study-table-no-probiotic-seqs.qza \
#     --o-visualization $deblur/current_study-table-no-probiotic-seqs.qzv

# Run through the rest of the data set with either the feature table with probiotic manure samples or without them

# for table in "$deblur/current_study-p-filtered-table.qza" "$deblur/current_study-table-no-probiotic.qza" ; do

#     # Set up an if statement to match the correct rep seqs file to the correct feature table

#     if [ "$table" = "$deblur/current_study-p-filtered-table.qza" ] ; then
#         repseqs="$deblur/current_study-p-filtered-seqs.qza"
#         name="pro"
#         refdiet='diet::cow_manure'

#     else
#         repseqs="$deblur/current_study-table-no-probiotic-seqs.qza"
#         name="no_pro"
#         refdiet='diet::gainesville'

#     fi

#     #Build taxonomy using GreenGenes 2 classifier

#     echo "Building the taxonomy file for current_study."

#     qiime feature-classifier classify-sklearn \
#         --i-classifier $phylogenetics/gg_classifier.qza \
#         --i-reads $repseqs \
#         --o-classification $taxonomy/current_study-$name-taxonomy.qza


#     # Export the taxonomy for later use
#     qiime tools export \
#         --input-path $taxonomy/current_study-$name-taxonomy.qza \
#         --output-path $taxonomy/cs-$name-taxa/

#     sed 's/^Feature ID/Feature_ID/' $taxonomy/cs-$name-taxa/taxonomy.tsv > $taxonomy/cs-$name-taxa/taxonomy-r.tsv

#     #Generating tree

#     echo "Starting generation of the IQ phylogenetic tree."

#     qiime phylogeny align-to-tree-mafft-iqtree \
#         --i-sequences $repseqs \
#         --p-n-threads 6 \
#         --p-seed 76 \
#         --o-alignment $phylogenetics/current_study-aligned-$name-seqs.qza \
#         --o-masked-alignment $phylogenetics/current_study-mask_aligned-$name-seqs.qza \
#         --o-tree $phylogenetics/current_study-$name-tree.qza \
#         --o-rooted-tree $phylogenetics/current_study-$name-rooted_tree.qza \

#     echo "Starting generation of the fragment-insertion phylogenetic tree."

#     qiime fragment-insertion sepp \
#         --i-representative-sequences $repseqs \
#         --i-reference-database $phylogenetics/sepp-refs-silva-128.qza \
#         --p-threads 6 \
#         --o-tree $phylogenetics/current_study-$name-frag_in-tree.qza \
#         --o-placements $phylogenetics/current_study-$name-frag_in-tree-placements.qza

#     # Filter the tree to match the feature table

#     echo "Filtering the phylogenetic trees to match the feature table."

#     qiime phylogeny filter-tree \
#         --i-tree $phylogenetics/current_study-$name-icamp_tree.qza \
#         --i-table $table \
#         --o-filtered-tree $phylogenetics/current_study-$name-icamp_tree-filtered.qza

#     qiime tools export \
#         --input-path $phylogenetics/current_study-$name-icamp_tree-filtered.qza \
#         --output-path $phylogenetics/current_study-$name-icamp_tree-filtered/

#     qiime phylogeny filter-tree \
#         --i-tree $phylogenetics/current_study-$name-frag_in-tree.qza \
#         --i-table $table \
#         --o-filtered-tree $phylogenetics/current_study-$name-frag_in-filtered-tree.qza

#     qiime tools export \
#         --input-path $phylogenetics/current_study-$name-frag_in-filtered-tree.qza \
#         --output-path $phylogenetics/current_study-$name-frag_in-filtered-tree/

    # ##Creating the Core Metrics Files
    # #Core metrics

    # echo "Beginning core metrics analysis on current_study."

    # qiime diversity core-metrics-phylogenetic \
    #     --i-table $table \
    #     --i-phylogeny $phylogenetics/current_study-$name-frag_in-filtered-tree.qza \
    #     --p-sampling-depth 100000 \
    #     --p-n-jobs-or-threads 0 \
    #     --m-metadata-file $metadata/bsf_saho_metadata.txt \
    #     --output-dir $base_dir/core-metrics-phylogeny-current_study-168k-$name-260128/

    # #Generating Generalized UniFrac distances on the unfiltered table
    # qiime diversity beta-phylogenetic \
    #     --i-table $table \
    #     --i-phylogeny $phylogenetics/current_study-$name-frag_in-filtered-tree.qza \
    #     --p-alpha 0.5 \
    #     --p-metric generalized_unifrac \
    #     --p-threads 0 \
    #     --o-distance-matrix $base_dir/core-metrics-phylogeny-current_study-168k-$name-260128/generalized_unifrac.qza

    # #Generating PCOA for the Unfiltered Generalized UniFrac
    # qiime diversity pcoa \
    #     --i-distance-matrix $base_dir/core-metrics-phylogeny-current_study-168k-$name-260128/generalized_unifrac.qza \
    #     --o-pcoa $base_dir/core-metrics-phylogeny-current_study-168k-$name-260128/generalized_unifrac_pcoa.qza

    # #Visualizing the Unfiltered Generalized PCOA
    # qiime emperor plot \
    #     --i-pcoa $base_dir/core-metrics-phylogeny-current_study-168k-$name-260128/generalized_unifrac_pcoa.qza \
    #     --m-metadata-file $metadata/bsf_saho_metadata.txt \
    #     --o-visualization $base_dir/core-metrics-phylogeny-current_study-168k-$name-260128/generalized_unifrac_emperor.qzv

    # qiime tools export \
    #     --input-path $table \
    #     --output-path $deblur/current_study-$name-table/
    # biom convert \
    #     -i $deblur/current_study-$name-table/feature-table.biom \
    #     -o $deblur/current_study-$name-table/feature-table.tsv \
    #     --to-tsv

    # sed '1d; s/^#OTU ID/Feature_ID/' $deblur/current_study-$name-table/feature-table.tsv > $deblur/current_study-$name-table/feature-table-r.tsv

    # Run ANCOM-BC2 for differential abundance on current study samples

    # qiime composition ancombc2 \
    #     --i-table $table \
    #     --m-metadata-file $metadata/bsf_metadata.txt \
    #     --p-reference-levels $refdiet \
    #     --p-fixed-effects-formula diet \
    #     --p-group diet \
    #     --o-ancombc2-output $ancom/saho_bsfl-$name-ancom.qza

    # qiime composition ancombc2-visualizer \
    #     --i-data $ancom/saho_bsfl-$name-ancom.qza \
    #     --i-taxonomy taxonomy/current_study-pro-taxonomy.qza \
    #     --o-visualization $ancom/saho_bsfl-$name-ancom.qzv

    # qiime tools export \
    #     --input-path $ancom/saho_bsfl-$name-ancom.qza \
    #     --output-path $ancom/saho_bsfl-$name-ancom/


    # qiime composition ancombc2 \
    #     --i-table $table \
    #     --m-metadata-file $metadata/bsf_metadata.txt \
    #     --p-reference-levels $refdiet \
    #     --p-fixed-effects-formula diet \
    #     --p-group diet \
    #     --p-structural-zeros "True" \
    #     --o-ancombc2-output $ancom/saho_bsfl-$name-ancom-str_zero.qza

    # qiime composition ancombc2-visualizer \
    #     --i-data $ancom/saho_bsfl-$name-ancom-str_zero.qza \
    #     --i-taxonomy taxonomy/current_study-pro-taxonomy.qza \
    #     --o-visualization $ancom/saho_bsfl-$name-ancom-str_zero.qzv

    # qiime tools export \
    #     --input-path $ancom/saho_bsfl-$name-ancom-str_zero.qza \
    #     --output-path $ancom/saho_bsfl-$name-ancom-str_zero/


# sed '2d; s/^#SampleID/SampleID/' "$metadata/bsf_saho_metadata.txt" > "$metadata/bsf_saho_metadata-r.txt"

# sed 's/^#SampleID/SampleID/; s/^#q2types/#q2:types/' "$metadata/bsf_saho_metadata.txt" > "$metadata/bsf_saho_metadata-q2r.txt"

# Run beta-group significance

    # for matrix in unweighted_unifrac weighted_unifrac; do
    
    #     qiime diversity beta-group-significance \
    #         --i-distance-matrix $base_dir/core-metrics-phylogeny-current_study-168k-$name-260128/${matrix}_distance_matrix.qza \
    #         --m-metadata-file $metadata/bsf_saho_metadata.txt \
    #         --m-metadata-column "diet" \
    #         --p-method permdisp \
    #         --p-pairwise TRUE \
    #         --o-visualization $base_dir/core-metrics-phylogeny-current_study-168k-$name-260128/${name}-${matrix}-permdisp.qzv
    
    #     qiime diversity beta-group-significance \
    #         --i-distance-matrix $base_dir/core-metrics-phylogeny-current_study-168k-$name-260128/${matrix}_distance_matrix.qza \
    #         --m-metadata-file $metadata/bsf_saho_metadata.txt \
    #         --m-metadata-column "diet" \
    #         --p-method permanova \
    #         --p-pairwise TRUE \
    #         --o-visualization $base_dir/core-metrics-phylogeny-current_study-168k-$name-260128/${name}-${matrix}-permanova.qzv

    # done

# Run alpha-group significance
    
    # for metric in evenness faith_pd shannon; do
    
    #     qiime diversity alpha-group-significance \
    #         --i-alpha-diversity $base_dir/core-metrics-phylogeny-current_study-168k-$name-260128/${metric}_vector.qza \
    #         --m-metadata-file $metadata/bsf_saho_metadata.txt \
    #         --o-visualization $base_dir/core-metrics-phylogeny-current_study-168k-$name-260128/${metric}_vector.qzv

    # done
    
# Run a core features analysis

    # qiime feature-table core-features \
    #     --i-table $table \
    #     --p-min-fraction 0.5 \
    #     --p-max-fraction 1.0 \
    #     --o-visualization $deblur/current_study-${name}core_features.qzv

# done

# Run a core features analysis without initial larvae samples

    qiime feature-table core-features \
        --i-table $deblur/current_study_no_initial.qza \
        --p-min-fraction 0.5 \
        --p-max-fraction 1.0 \
        --o-visualization $deblur/current_study_no_initial-core_features.qzv


# Rscript $scripts/qiime/current_study-figs.R

# # Filter the feature table and rep seqs to only include Gainesville and IL

# qiime feature-table filter-samples \
#     --i-table $deblur/current_study-table-no-probiotic.qza \
#     --m-metadata-file $metadata/bsf_saho_metadata.txt \
#     --p-where "[diet] = 'gainesville_initial'" \
#     --p-exclude-ids \
#     --o-filtered-table $deblur/current_study_no_initial.qza
    
# qiime feature-table filter-samples \
#     --i-table $deblur/current_study-table-no-probiotic.qza \
#     --m-metadata-file $metadata/bsf_saho_metadata.txt \
#     --p-where "[diet] = 'gainesville_initial' OR [diet] = 'gainesville'" \
#     --o-filtered-table $deblur/current_study_gv_initial.qza

# for table in $deblur/current_study_no_initial.qza $deblur/current_study_gv_initial.qza; do

#     if [ "$table" = "$deblur/current_study_no_initial.qza" ] ; then

#         name="no_initial"

#     else
#         name="gv_initial"

#     fi

#     qiime feature-table filter-seqs \
#     --i-data $deblur/current_study-table-no-probiotic-seqs.qza \
#     --i-table $table \
#     --o-filtered-data $deblur/${name}_repseqs.qza 

# done
