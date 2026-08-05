# This code was written by Reese Saho for the black soldier fly meta-analysis study. 2025-07-13

# Initialize conda

eval "$(/lustre/isaac24/scratch/rsaho/Miniconda/bin/conda shell.bash hook)"


# Making necessary directories

base_dir="/lustre/isaac24/proj/UTK0406/ReeseSaho/projects/bsf_metaa"
raw_seqs="$base_dir/raw_sequencing_files"
demux="$base_dir/demux"
deblur="$base_dir/deblur"
scripts="$base_dir/scripts"
logs="$base_dir/logs"
phylogenetics="$base_dir/phylogenetics"
taxonomy="$base_dir/taxonomy"
rarefaction="$base_dir/rarefaction"
metadata="$base_dir/metadata"
networks="$base_dir/networks"
correction="$base_dir/batch_correction"
figures="$base_dir/figures"
tools="/lustre/isaac24/scratch/rsaho/tools"

export scripts raw_seqs demux deblur logs phylogenetics taxonomy rarefaction metadata networks correction figures tools

# Create directories

mkdir -p "$demux"
mkdir -p "$deblur"
mkdir -p "$phylogenetics"
mkdir -p "$taxonomy"
mkdir -p "$rarefaction"
mkdir -p "$logs"
mkdir -p "$metadata"
mkdir -p "$networks"
mkdir -p "$figures"

mkdir -p "$raw_seqs/paired"
mkdir -p "$raw_seqs/single-end"
mkdir -p "$raw_seqs/merged_files"
mkdir -p "$raw_seqs/unmerged_files"

mkdir -p "$scripts/qiime"
mkdir -p "$scripts/batch_correction"
mkdir -p "$scripts/networks"
mkdir -p "$scripts/icamp"

cd $base_dir

# I. SAMPLE PROCESSING

# Activate the SRA Toolkit Conda environment

conda activate sra-tools

# Transfer accession files to respective directories and import SRAs. Note that current study samples were already imported

cd $raw_seqs/paired

paired_accession_file=$raw_seqs/paired_master_accession_list-new.csv

while IFS= read -r line; do
    # Skip empty lines or comment lines
    [[ -z "$line" || "$line" =~ ^# ]] && continue

    accession=$(echo "$line" | tr -d '\r')

    echo "Fetching: $accession"
    prefetch "$accession"

    echo "Converting to FASTQ: $accession"
    fasterq-dump --split-files "$accession"
    gzip "${accession}"*.fastq

done < "$paired_accession_file"


cd $raw_seqs/single-end

single_accession_file=$raw_seqs/single_master_accession_list.csv

while IFS= read -r line; do
    # Skip empty lines or comment lines
    [[ -z "$line" || "$line" =~ ^# ]] && continue

    accession=$(echo "$line" | tr -d '\r')

    echo "Fetching: $accession"
    prefetch "$accession"

    echo "Converting to FASTQ: $accession"
    fasterq-dump --split-files "$accession"
    gzip "${accession}"*.fastq

done < "$single_accession_file"

# Activate the QIIME2 Conda environment

conda deactivate

conda activate /lustre/isaac24/proj/UTK0406/shared_conda_envs/qiime2-amplicon-2025.4

# Run the paired-end processing script

source $scripts/qiime/paired_processing.sh

# Run the single-end processing script

source $scripts/qiime/single_processing.sh

# Run the denoising - core metrics script

source $scripts/qiime/processing_and_core_metrics.sh

# # Run the batch correction script

conda deactivate

conda activate conqur-r

sed '2d; s/^#SampleID/SampleID/' "$metadata/bsf_metadata.txt" > "$metadata/bsf_metadata-r.txt"

Rscript $scripts/batch_correction/ConQuR-Default.R \
    --tree $phylogenetics/bsf-tree.qza \
    --table $deblur/bsf-feature_table-final.qza \
    --metadata $metadata/bsf_metadata-r.txt \
    --output_dir $base_dir/batch_correction/

Rscript $scripts/batch_correction/ConQuR-Default-PERMANOVA-R2.R \
    --output_dir $base_dir/batch_correction/

# Run core metrics and adonis on post corrected data

conda deactivate

conda activate /lustre/isaac24/proj/UTK0406/shared_conda_envs/qiime2-amplicon-2025.4

source $scripts/qiime/bsf_post_correction.sh

# Run the Saho processing script

source $scripts/qiime/saho_processing.sh

# Make the networks

conda deactivate

conda activate sparcc-env

bash $networks/meta_a_networks/run_scripts.sh \


    
    
