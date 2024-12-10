#!/bin/bash

> res.out

# echo "                               |           SWG           |         Windowed        |          Banded         " > res.out
# echo " WS|OS| Pens  | Len  | Sim/Real| Elap  | Mem| Cells|Score| Elap  | Mem| Cells|Score| Elap  | Mem| Cells|Score" >> res.out

# Read each line of all_tests.txt
while IFS=, read -r file type; do
    
    # Skip lines starting with "//"
    [[ "$file" =~ ^// ]] && continue
    
    # Execute the command with each line's parameters
    ./quickaffine -i "test_datasets/$file" -o res.out -a SWG+Windowed+Banded -p Bowtie2 -ws 32 -os 8 -t "$type"
    ./quickaffine -i "test_datasets/$file" -o res.out -a SWG+Windowed+Banded -p Bowtie2 -ws 64 -os 16 -t "$type"
    ./quickaffine -i "test_datasets/$file" -o res.out -a SWG+Windowed+Banded -p Bowtie2 -ws 128 -os 32 -t "$type"
    
done < test_datasets/all_tests.txt
