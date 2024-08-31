#!/bin/bash

# Two cbImport objects - one for each assay - CITE
INPUT_cbImport_directory_RNA="$1"
INPUT_cbImport_directory_ADT="$2"

# Function to unzip, convert files, and apply specific modifications based on directory
convert_files() {
  local directory="$1"
  local assay_type="$2"
  
  # Specify the files to process
  local files=("counts_exprMatrix.tsv.gz" "data_exprMatrix.tsv.gz")
  
  for file in "${files[@]}"; do
    # Construct the full path to the file
    local full_path="${directory}/${file}"
    
    # Check if the file exists
    if [[ -f "$full_path" ]]; then
      # Unzip the file while keeping the original (.gz)
      gunzip -k "$full_path"
      
      # Remove the .gz extension from the filename
      local tsv_file="${full_path%.gz}"
      
      # Replace TSV extension with CSV for the output filename
      local csv_file="${tsv_file%.tsv}.csv"
      
      # Convert TSV to CSV by replacing tabs with commas
      sed 's/\t/,/g' "$tsv_file" > "$csv_file"
      
      # If processing ADT directory files, apply specific modifications
      if [[ "$assay_type" == "ADT" ]]; then
        # Prepend "_ADT" to every value in the first field, except the header
        awk 'BEGIN {FS=OFS=","} {if (NR==1) {print} else {$1="ADT_"$1; print}}' "$csv_file" > "${csv_file%.csv}_ADT_modified.csv"
        echo "Converted and modified: ${csv_file%.csv}_ADT_modified.csv"
      else
        echo "Converted: $csv_file"
      fi
    else
      echo "File not found: $full_path"
    fi
  done
}

# Process the directories with type indication
convert_files "$INPUT_cbImport_directory_RNA" "RNA"
convert_files "$INPUT_cbImport_directory_ADT" "ADT"


# After processing the directories, merge, convert, and gzip the final file
merge_and_compress() {
  # Merge the files, assuming the names follow a specific pattern and are located in the specified directories
  (head -n 1 "$INPUT_cbImport_directory_RNA/data_exprMatrix.csv" && \
   tail -n +2 "$INPUT_cbImport_directory_RNA/data_exprMatrix.csv" && \
   tail -n +2 "$INPUT_cbImport_directory_ADT/data_exprMatrix_ADT_modified.csv") > combined_RNA_ADT_norm_EXPdata.csv

  # Convert the merged CSV back into a TSV
  sed 's/,/\t/g' combined_RNA_ADT_norm_EXPdata.csv > combined_RNA_ADT_norm_EXPdata.tsv

  # GZIP the TSV
  gzip combined_RNA_ADT_norm_EXPdata.tsv

  echo "Merged RNA and labeled ADT assays, converted to TSV, and gzipped successfully."
}

# Call the new function after processing the individual files
merge_and_compress


# ARGUMENT CALL: bash Cell_Browser_cbImport_Processing.sh <path/to/cbImport/RNA-assay/directory> <path/to/cbImport/ADT-assay/directory>

