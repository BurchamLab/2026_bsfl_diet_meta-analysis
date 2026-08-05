
variable=""
Variable=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --variable)  variable="$2"; shift 2 ;;
    --Variable)  Variable="$2"; shift 2 ;;
    *) echo "Unknown option: $1"; exit 1 ;;
  esac
done

# Safety check
if [[ -z "$variable" || -z "$Variable" ]]; then
    echo "Error: --variable and --Variable must be provided. --Variable is an uppercase version of --variable used for naming in figures."
    exit 1
fi

for file in $scripts/networks/meta_a_networks/base-pre_sparcc.R $scripts/networks/meta_a_networks/base-network.R $scripts/networks/meta_a_networks/base-figures.R; do

    new_file_name="${file/base/$variable}"

    sed \
    -e "s/variable/$variable/g" \
    -e "s/Variable/$Variable/g" \
    "$file" > "$new_file_name"

done

echo "Starting analysis on $Variable samples."

Rscript $scripts/networks/meta_a_networks/$variable-pre_sparcc.R

i=1

for table in $networks/$variable/$variable-table_*.txt; do
    output_dir=$networks/$variable/subsample/${i}
    mkdir -p "$output_dir"

    python /lustre/isaac24/scratch/rsaho/tools/SparCC/$variable-General_Execution.py \
        --data-input "$table" \
        --outpath "$output_dir" \
        --outfile-pvals "${output_dir}/pvals_one_sided.csv" \
        --save-corr-file "${output_dir}/cor_sparcc.csv" 

done

Rscript $scripts/networks/meta_a_networks/$variable-network.R
Rscript $scripts/networks/meta_a_networks/$variable-figures.R

echo "Analysis has been completed."
