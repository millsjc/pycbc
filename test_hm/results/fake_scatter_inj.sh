#!/bin/bash

# Check if FOUNDMISSED_FILE argument is passed
if [ -z "$1" ]
then
    echo "Error: FOUNDMISSED_FILE argument is missing."
    exit 1
fi

# Set OUTPUTFILE variable with changed extension to .png
FOUNDMISSED_FILE=$1
OUTPUTFILE="${FOUNDMISSED_FILE%.*}.png"

# Run pycbc_hm_scatter_inj_results with the updated arguments
../../bin/hm/pycbc_hm_scatter_inj_results --inj-file "${FOUNDMISSED_FILE}" --bank-file "../bank/mini-bank.hdf" --output-file "${OUTPUTFILE}"
