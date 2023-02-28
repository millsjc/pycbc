#!/bin/bash

# Parse command line arguments
injfile=$1
fixed_snr=$2
inplace=$3

if [ "$inplace" = "--inplace" ]; then
    # Use Python script to modify distance columns in-place
    python << END
import h5py

with h5py.File("$injfile", "r+") as f:
    fixed_snr = $fixed_snr
    for key in f.keys():
        if "distance" in key or "snr" in key:
            f[key][()] = fixed_snr / f["network_optimal_snr"][()] * f[key][()]
END

else
    # Create a new file with the same contents
    newfile="${injfile%.*}_${fixed_snr}.${injfile##*.}"
    cp "$injfile" "$newfile"
    
    # Use Python script to modify distance columns in the new file
    python << END
import h5py

with h5py.File("$newfile", "r+") as f:
    fixed_snr = $fixed_snr
    for key in f.keys():
        if "distance" in key or "snr" in key:
            f[key][()] = fixed_snr / f["network_optimal_snr"][()] * f[key][()]
END
    
    # Print the name of the new file
    echo "Created new file: $newfile"
fi

