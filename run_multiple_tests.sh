#!/bin/bash

> res.out

# echo "                               |           SWG           |         Windowed        |          Banded         |       Parasail Scan     |       Parasail Diag     " > res.out
# echo " WS|OS| Pens  | Len  | Sim/Real| Elap  | Mem| Cells|Score| Elap  | Mem| Cells|Score| Elap  | Mem| Cells|Score| Elap  | Mem| Cells|Score| Elap  | Mem| Cells|Score" >> res.out

# Read each line of all_tests.txt
while IFS=, read -r file type; do
    
    # Skip lines starting with "//"
    [[ "$file" =~ ^// ]] && continue
    
    # Execute the command with each line's parameters
    echo "Running test_datasets/v3/$file"
    # ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed+Banded+Parasail -p Bowtie2 -ws 32 -os 8 -t "$type"
    # ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed+Banded+Parasail -p Bowtie2 -ws 32 -os 16 -t "$type"
    # ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed+Banded+Parasail -p Bowtie2 -ws 64 -os 16 -t "$type"
    # ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed+Banded+Parasail -p Bowtie2 -ws 64 -os 32 -t "$type"
    # ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed+Banded+Parasail -p Bowtie2 -ws 128 -os 32 -t "$type"
    # ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed+Banded+Parasail -p Bowtie2 -ws 128 -os 64 -t "$type"

    ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed+Banded -p Bowtie2 -ws 32 -os 8 -t "$type"
    ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed+Banded -p Bowtie2 -ws 32 -os 16 -t "$type"
    ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed+Banded -p Bowtie2 -ws 64 -os 16 -t "$type"
    ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed+Banded -p Bowtie2 -ws 64 -os 32 -t "$type"
    ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed+Banded -p Bowtie2 -ws 128 -os 32 -t "$type"
    ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed+Banded -p Bowtie2 -ws 128 -os 64 -t "$type"

    # Run quickaffine for scores.out results
    # ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed -p Bowtie2 -ws 32 -os 8 -t "$type"
    # ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed -p Bowtie2 -ws 32 -os 16 -t "$type"
    # ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed -p Bowtie2 -ws 64 -os 16 -t "$type"
    # ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed -p Bowtie2 -ws 64 -os 32 -t "$type"
    # ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed -p Bowtie2 -ws 128 -os 32 -t "$type"
    # ./quickaffine -i "test_datasets/v3/$file" -o res.out -so scores.out -a SWG+Windowed -p Bowtie2 -ws 128 -os 64 -t "$type"
    
done < test_datasets/v3/all_tests.txt
