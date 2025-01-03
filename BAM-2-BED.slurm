#!/bin/bash
Recommended nodes 5
Recommended cpus/task 20
# Define input directory containing BAM files
INPUT_DIR="all5seqbams"
OUTPUT_DIR="all5seqBeds"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through each BAM file in the input directory
for BAM_FILE in "$INPUT_DIR"/*.bam; do
    if [ ! -f "$BAM_FILE" ]; then
        echo "No BAM files found in $INPUT_DIR."
        exit 1
    fi
    
    BASE_NAME=$(basename "$BAM_FILE" .bam)            # Base name for output files

    # Define output files based on the input BAM file name, directing them to the output directory
    BEDGRAPH_POS="${OUTPUT_DIR}/${BASE_NAME}_positive.bedgraph"        # Output BEDGraph for positive strand
    BEDGRAPH_NEG="${OUTPUT_DIR}/${BASE_NAME}_negative.bedgraph"        # Output BEDGraph for negative strand
    BEDGRAPH_STRANDED="${OUTPUT_DIR}/${BASE_NAME}_stranded.bedgraph"    # Combined stranded BEDGraph
    BED_FILE="${OUTPUT_DIR}/${BASE_NAME}.bed"                            # Output BED file

    # Step 1: Calculate total mapped reads for normalization (RPM)
    echo "Processing $BAM_FILE..."
    TOTAL_READS=$(samtools view -c -F 4 "$BAM_FILE")
    SCALING_FACTOR=$(echo "scale=6; 1000000 / $TOTAL_READS" | bc)

    echo "Total reads for $BASE_NAME: $TOTAL_READS"
    echo "Scaling factor: $SCALING_FACTOR"

    # Step 2: Generate strand-specific normalized coverage in BEDGraph format
    echo "Generating strand-specific normalized BEDGraph for $BASE_NAME..."
    bedtools genomecov -ibam "$BAM_FILE" -bg -scale "$SCALING_FACTOR" -strand + > "$BEDGRAPH_POS"
    bedtools genomecov -ibam "$BAM_FILE" -bg -scale "$SCALING_FACTOR" -strand - > "$BEDGRAPH_NEG"

    # Annotate strands in the BEDGraph files and combine them
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t+"}' "$BEDGRAPH_POS" > "${OUTPUT_DIR}/${BASE_NAME}_positive_annotated.bedgraph"
    awk '{print $1"\t"$2"\t"$3"\t"$4"\t-"}' "$BEDGRAPH_NEG" > "${OUTPUT_DIR}/${BASE_NAME}_negative_annotated.bedgraph"

    cat "${OUTPUT_DIR}/${BASE_NAME}_positive_annotated.bedgraph" "${OUTPUT_DIR}/${BASE_NAME}_negative_annotated.bedgraph" > "$BEDGRAPH_STRANDED"

    # Convert BAM to BED using bedtools
    echo "Converting BAM to BED for $BASE_NAME..."
    bedtools bamtobed -i "$BAM_FILE" > "$BED_FILE"

    # Check if the conversion was successful
    if [ -s "$BED_FILE" ]; then
        echo "Conversion successful for $BASE_NAME. Output file: $BED_FILE"
    else
        echo "Conversion failed for $BASE_NAME. Output file is empty."
    fi

done

echo "Processing complete for all BAM files."

