#!/bin/bash
set -e

# Input directories
INPUT_DIRS=(
  "/burg/cgl/users/ls3456_v2/20250426_crispra_pcr_hash/1-output-fastq-correct"
  "/burg/cgl/users/ls3456_v2/20240930_CRISPRA_CRISPRI/1-output-fastq_a"
  "/burg/cgl/users/ls3456_v2/20251023_CRISPRA_cDNA/1-output-fastq"
  "/burg/cgl/users/ls3456_v2/20240805_CRISPRA_CRISPRI/1-output-fastq_a"
)

# Output directory
OUTDIR="/burg/cgl/users/ls3456_v2/GBM_Tcells_CRISPRA/hash/combined_fastqs"

# Clean up previous output files FIRST before anything else
rm -rf "$OUTDIR"
mkdir -p "$OUTDIR"

# Collect all matching FASTQ files
ALL_FILES=()
for DIR in "${INPUT_DIRS[@]}"; do
  while IFS= read -r FILE; do
    [[ -n "$FILE" ]] && ALL_FILES+=("$FILE")
  done < <(find "$DIR" -type f \( -name '03C*R1*.fastq.gz' -o -name '04C*R1*.fastq.gz' -o -name '03C*R2*.fastq.gz' -o -name '04C*R2*.fastq.gz' \) 2>/dev/null)
done

echo "Found ${#ALL_FILES[@]} files total"

# Initialize grouping
declare -A FILE_GROUPS

# Group files by well ID (A01, A02, ..., H12) and read direction
for FILE in "${ALL_FILES[@]}"; do
  FILENAME=$(basename "$FILE")
  
  if [[ "$FILENAME" =~ ^(03C|04C)([A-H][0-9]{2})_S[0-9]+_R([12])_001\.fastq\.gz$ ]]; then
    WELL="${BASH_REMATCH[2]}"
    READ="R${BASH_REMATCH[3]}"
    KEY="${WELL}_${READ}"
    FILE_GROUPS["$KEY"]+="$FILE "
  else
    echo "  No match: $FILENAME"
  fi
done

echo "Found ${#FILE_GROUPS[@]} groups"

# Combine files by group
for KEY in "${!FILE_GROUPS[@]}"; do
  OUTFILE="${OUTDIR}/03C${KEY}.fastq.gz"
  echo "Writing: $OUTFILE"
  echo "${FILE_GROUPS[$KEY]}" | tr ' ' '\n' | grep -v '^$' | xargs cat > "$OUTFILE"
done

echo "✅ All 03C/04C FASTQ files combined by well ID and read pair."
