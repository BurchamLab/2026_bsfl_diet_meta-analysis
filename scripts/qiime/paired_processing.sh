# This code was writted by Reese Saho for the purpose of processing paired-end samples to be used in the black soldier fly meta-analysis

# Import the paired end sequencing data for each study


max_jobs=23
count=0

for study in liu2020b querejeta2022 li2022 silvaraju2024 tanga2021 wu2020 kluber2022 klammsteiner2025 edwards2024 zhang2020 galassi2021 wu2021 chen2023; do
# liu2020 dragone2025 liu2021 deng2024 cifuentes2020 piersanti2024 bruno2019 yu2023 tegtmeier2021a tegtmeier2021b defilippis2023 deng2025 klammsteiner2020 klammsteiner2021 pei2022 vandeweyer2023 wynants2018 current_study reyer2025 gorrens2022 li2025 marasco2022 

    echo "Importing paired-end FASTQ files for $study."

    if [ $study = "current_study" ] ; then
        :
    else
        # Convert FASTQ files to a QIIME2 artifact

        qiime tools import \
            --type 'SampleData[PairedEndSequencesWithQuality]' \
            --input-path $raw_seqs/${study}_paired_manifest.txt \
            --input-format PairedEndFastqManifestPhred33V2 \
            --output-path $demux/${study}_demux.qza

        qiime demux summarize \
            --i-data $demux/${study}_demux.qza \
            --o-visualization $demux/${study}_demux.qzv

    fi
        # Trim 3' ends of each sample in order to maximize merging efficiency

    echo "Trimming 3' ends of samples from $study."

    if [ $study = "klammsteiner2020" ] || [ $study = "klammsteiner2021" ] ; then #V4 only

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -10 \
            --p-reverse-cut -50 
        
        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv


    elif [ $study = "current_study" ] || [ $study = "wynants2018" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -50 \
            --p-reverse-cut -50 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "yu2023" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -10 \
            --p-reverse-cut -10 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "reyer2025" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -2 \
            --p-reverse-cut -20 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv


    elif [ $study = "marasco2022" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -20 \
            --p-reverse-cut -80 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "pei2022" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -10 \
            --p-reverse-cut -60 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "liu2020" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -10 \
            --p-reverse-cut -50 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "piersanti2024" ] ; then # Lot of low quality seqs

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -20 \
            --p-reverse-cut -85 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "gorrens2022" ] ; then # Just doesn't merge very well

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -90 \
            --p-reverse-cut -90 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv
    elif [ $study = "deng2024" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -5 \
            --p-reverse-cut -5 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "tegtmeier2021a" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -15 \
            --p-reverse-cut -20 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv
    
    elif [ $study = "vandeweyer2023" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -60 \
            --p-reverse-cut -60 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "deng2025" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -1 \
            --p-reverse-cut -1 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "deng2024" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -1 \
            --p-reverse-cut -15 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "dragone2025" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -1 \
            --p-reverse-cut -1 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "li2025" ] ; then # Just doesn't merge very well

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -2 \
            --p-reverse-cut -2 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv
    
    elif [ $study = "bruno2019" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -1 \
            --p-reverse-cut -10 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "liu2021" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -1 \
            --p-reverse-cut -1 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv
    
    elif [ $study = "tegtmeier2021b" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -15 \
            --p-reverse-cut -10 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "cifuentes2020" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -20 \
            --p-reverse-cut -35 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "zhang2020" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -20 \
            --p-reverse-cut -100 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "li2022" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -1 \
            --p-reverse-cut -1 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv


    elif [ $study = "chen2023" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -10 \
            --p-reverse-cut -80 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "kluber2022" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -45 \
            --p-reverse-cut -35 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "querejeta2022" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -10 \
            --p-reverse-cut -30 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "wu2020" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -40 \
            --p-reverse-cut -40 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "wu2021" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -5 \
            --p-reverse-cut -5 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "liu2020b" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -1 \
            --p-reverse-cut -1 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "galassi2021" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -20 \
            --p-reverse-cut -40 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "tanga2021" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -40 \
            --p-reverse-cut -80 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    elif [ $study = "silvaraju2024" ] ; then

        qiime cutadapt trim-paired \
            --i-demultiplexed-sequences $demux/${study}_demux.qza \
            --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
            --p-forward-cut -40 \
            --p-reverse-cut -100 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv

    else # Studies that don't cover the V3-V4 region should not be trimmed or they are unable to successfully merge

        qiime cutadapt trim-paired \
        --i-demultiplexed-sequences $demux/${study}_demux.qza \
        --o-trimmed-sequences $demux/${study}_trimmed_seqs.qza \
        --p-forward-cut -10 \
        --p-reverse-cut -10 

        qiime demux summarize \
            --i-data $demux/${study}_trimmed_seqs.qza \
            --o-visualization $demux/${study}_trimmed_seqs.qzv


    fi

# Quality Filtering

echo "Quality filtering samples from $study."

    qiime quality-filter q-score \
        --i-demux $demux/${study}_trimmed_seqs.qza \
        --p-min-quality 25 \
        --o-filtered-sequences $demux/${study}_filtered_seqs.qza \
        --o-filter-stats $demux/${study}_filtered_stats.qza 

    # Summarize filtered seqs
    qiime demux summarize \
        --i-data $demux/${study}_filtered_seqs.qza \
        --o-visualization $demux/${study}_filtered_seqs.qzv

    # Summarize filtered stats
    qiime metadata tabulate \
        --m-input-file $demux/${study}_filtered_stats.qza \
        --o-visualization $demux/${study}_filtered_stats.qzv 

    # Merge the paired end reads together

    echo "Merging samples from $study."

        qiime vsearch merge-pairs \
            --i-demultiplexed-seqs $demux/${study}_filtered_seqs.qza \
            --p-minovlen 20 \
            --o-merged-sequences $demux/${study}_merged-demux.qza \
            --o-unmerged-sequences $demux/${study}_unmerged-demux.qza

        qiime demux summarize \
            --i-data $demux/${study}_merged-demux.qza \
            --o-visualization $demux/${study}_merged-demux.qzv


    # Export merged FASTQ files so that they can be reimported as single-end for Deblur

    echo "Exporting samples from $study."

    qiime tools export \
        --input-path $demux/${study}_merged-demux.qza \
        --output-path $demux/${study}_merged/

    # Re-import merged FASTQ files as single-end sequences

    echo "Importing samples from $study as single-end sequences."

    qiime tools import \
        --type 'SampleData[SequencesWithQuality]' \
        --input-path $raw_seqs/${study}_single_manifest.txt \
        --input-format SingleEndFastqManifestPhred33V2 \
        --output-path $demux/${study}-single-demux.qza

    echo "Finished paired-end processing on $study samples."
        echo "-----------------------------------------"

done

