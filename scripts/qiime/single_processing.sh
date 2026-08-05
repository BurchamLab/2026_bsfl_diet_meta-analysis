# This code was writted by Reese Saho for the purpose of processing single-end samples to be used in the black soldier fly meta-analysis

# Import the single end sequencing data for each study


max_jobs=6
count=0

for study in li2021 jiang2019 greenwood2021 shelomi2020 xu2025 ; do


    echo "Importing single-end FASTQ files for $study."

    # Convert FASTQ files to a QIIME2 artifact

    qiime tools import \
        --type 'SampleData[SequencesWithQuality]' \
        --input-path $raw_seqs/${study}_single_manifest.txt \
        --input-format SingleEndFastqManifestPhred33V2 \
        --output-path $demux/${study}_demux.qza

    qiime demux summarize \
        --i-data $demux/${study}_demux.qza \
        --o-visualization $demux/${study}_demux.qzv

    # Quality Filtering

    qiime quality-filter q-score \
        --i-demux $demux/${study}_demux.qza \
        --p-min-quality 25 \
        --o-filtered-sequences $demux/${study}-single-demux.qza \
        --o-filter-stats $demux/${study}_filtered_stats.qza 

    # Summarize filtered seqs
    qiime demux summarize \
        --i-data $demux/${study}-single-demux.qza  \
        --o-visualization $demux/${study}-single-demux.qzv

    # Summarize filtered stats
    qiime metadata tabulate \
        --m-input-file $demux/${study}_filtered_stats.qza \
        --o-visualization $demux/${study}_filtered_stats.qzv 

done

