#!/bin/bash

> res.out

# Check if the correct number of arguments is provided
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <Cm> <Cx> <Co> <Ci> <Cd>"
    exit 1
fi

# Assign command line arguments to variables
Cm=$1
Cx=$2
Co=$3
Ci=$4
Cd=$5

# Read each line of all_tests.txt
while IFS=, read -r file size real_size error; do
    
    # Skip lines starting with "//"
    [[ "$file" =~ ^// ]] && continue
    
    # Execute the command with each line's parameters
    ./quickedaffine_align "test_datasets/$file" res.out "$size" "32" "8" "$Cm" "$Cx" "$Co" "$Ci" "$Cd" "$real_size" "$error"
    ./quickedaffine_align "test_datasets/$file" res.out "$size" "64" "16" "$Cm" "$Cx" "$Co" "$Ci" "$Cd" "$real_size" "$error"
    ./quickedaffine_align "test_datasets/$file" res.out "$size" "64" "32" "$Cm" "$Cx" "$Co" "$Ci" "$Cd" "$real_size" "$error"
    ./quickedaffine_align "test_datasets/$file" res.out "$size" "128" "32" "$Cm" "$Cx" "$Co" "$Ci" "$Cd" "$real_size" "$error"
done < test_datasets/all_tests.txt
