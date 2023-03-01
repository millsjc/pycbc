#!/bin/bash

# Parse command line arguments
injfile=$1
fixed_snr=$2
inplace=$3

if [ "$inplace" = "--inplace" ]; then
    # Use Python script to modify distance columns in-place
    python fix_snr.py --injfile $injfile --fixed-snr $fixed_snr

else
    # Create a new file with the same contents
    newfile="${injfile%.*}_${fixed_snr}.${injfile##*.}"
    cp "$injfile" "$newfile"
    
    # Use Python script to modify distance columns in the new file
    python fix_snr.py --injfile $newfile --fixed-snr $fixed_snr
    
    # Print the name of the new file
    echo "Created new file: $newfile"
fi

