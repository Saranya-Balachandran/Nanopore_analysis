#!/bin/bash


CSV_FILE="epi2me_data.csv"


if [ ! -f "$CSV_FILE" ]; then
    echo "Error: CSV file '$CSV_FILE' not found."
    exit 1
fi

read -r HEADER < "$CSV_FILE"
echo "Skipping header: $HEADER"
echo "-------------------------------------"

tail -n +2 "$CSV_FILE" | while IFS=',' read -r COL1 COL2 COL3
do
    # Check if a line was successfully read
    if [ -z "$COL1" ]; then
        continue
    fi

    # Trim leading/trailing whitespace and remove double quotes if necessary.
    
    P1=$(echo "$COL1" | sed 's/^"//;s/"$//')
    P2=$(echo "$COL2" | sed 's/^"//;s/"$//')
    P3=$(echo "$COL3" | sed 's/^"//;s/"$//;s/\r$//')

    echo "Passing Row to Target: ($P1) ($P2) ($P3)"

    sbatch epi2me.sh "$P1" "$P2" "$P3"

done

echo "-------------------------------------"
echo "CSV processing complete."
