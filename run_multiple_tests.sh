#!/bin/bash

> res.out

# Check if the correct number of arguments is provided
if [ "$#" -ne 7 ]; then
    echo "Usage: $0 <ws> <os> <Cm> <Cx> <Co> <Ci> <Cd>"
    exit 1
fi

# Assign command line arguments to variables
ws=$1
os=$2
Cm=$3
Cx=$4
Co=$5
Ci=$6
Cd=$7

# Read each line of all_tests.txt
while IFS=, read -r file size real_size error; do
    # Execute the command with each line's parameters
    ./quickedaffine_align "test_datasets/$file" res.out "$size" "$ws" "$os" "$Cm" "$Cx" "$Co" "$Ci" "$Cd" "$real_size" "$error"
done < test_datasets/all_tests.txt
