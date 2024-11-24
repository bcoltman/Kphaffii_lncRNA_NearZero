#!/usr/bin/env bash
#### Download GS115 GO files from GprofileR, convert to CBS7435
#### inputs are: 1) output directory 2) cbs7435_gs115 mapping

while getopts o:i: option
do
    case "${option}" in
    o) GO_DIR=${OPTARG};;
    i) KP_MAP=${OPTARG};;
    esac
done

mkdir -p "$GO_DIR"  # Create the output directory if it doesn't exist

# Download gprofiler file for GS115
wget --quiet -nc https://biit.cs.ut.ee/gprofiler//static/gprofiler_full_kpastoris.name.gmt -P "$GO_DIR"

# Define files
gene_list_file="$GO_DIR/gprofiler_full_kpastoris.name.gmt"
output_file="$GO_DIR/gprofiler_converted.gmt"

homolog_file="$1"
input_file="$2"

# Use awk to perform the replacement
awk -F'\t' 'NR==FNR{homolog[$2]=$1; next} {output=$1 "\t" $2; for (i=3; i<=NF; i++) { if ($i in homolog) { output = output "\t" homolog[$i]; } } print output }' OFS='\t' "$KP_MAP" "$gene_list_file" > "$output_file"

echo "Replacement complete. Output saved in $output_file"

