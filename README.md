# QuickedAffine

## Overview

**QuickedAffine** is a sequence alignment algorithm designed to efficiently compute optimal global alignments with gap-affine penalties using a combination of algebraic transformations, windowed heuristics, and exact banded gap-affine algorithms. This project is aimed at accelerating the sequence alignment process, making it suitable for long genomic sequences by reducing the time and memory complexity without sacrificing biological accuracy.

By combining innovative techniques such as windowed heuristics to estimate an upper-bound alignment score and an exact gap-affine alignment algorithm, QuickedAffine achieves a balance between speed and accuracy, outperforming traditional algorithms such as Smith-Waterman-Gotoh in large datasets.

## Key Features

- **Algebraic Transformation**: Converts alignment penalties into non-negative values for correct windowed gap-affine computation.
- **Windowed Heuristic**: A heuristic-based algorithm that calculates an upper-bound on the optimal score within small windows of the dynamic programming matrix, ensuring faster computation of large sequences.
- **Exact Banded Gap-Affine Algorithm**: Refines the solution within a computed band to guarantee that the final alignment score and the CIGAR string are biologically accurate.
- **Backtracking**: Efficient traceback mechanism to compute the final optimal alignment and CIGAR string.

## How It Works

The **QuickedAffine** algorithm follows three main phases:

1. **Algebraic Transformation**:
   - This step transforms all penalties to non-negative values, which ensures that the gap penalties and match/mismatch penalties can be processed consistently. It also eliminates issues with local optima that arise from having mixed positive and negative scores in the dynamic programming matrix.

2. **Windowed Heuristic**:
   - The windowed heuristic method restricts the alignment computation to small windows of the dynamic programming matrix, significantly reducing the memory requirements. This method estimates an upper-bound score for the alignment and progressively shifts the window to cover the entire sequence.

3. **Banded Gap-Affine Alignment**:
   - Once the upper-bound score is computed, the exact banded gap-affine algorithm refines the solution by working within a restricted band to calculate the optimal score and CIGAR string. This final phase ensures that the alignment produced is both optimal and biologically meaningful.

## Algorithms Used

- **Windowed Heuristic**: Used to reduce the computational overhead by focusing on smaller windows of the DP matrix.
- **Banded Gap-Affine Algorithm**: Ensures that the final computed alignment is accurate with biologically relevant gap penalties.
- **Exact Backtrace**: Once the alignment score is computed, the backtrace reconstructs the optimal alignment and outputs the CIGAR string.

## Installation

To use **QuickedAffine**, clone the repository and ensure you have the required dependencies installed. 

```bash
git clone https://github.com/your-username/QuickedAffine.git
cd QuickedAffine
```

## Usage

You can run the QuickedAffine algorithm on any compatible sequence dataset by executing the following command:

```bash
make
./quickaffine -i test_datasets/real/Nanopore_100s.seq -o res.out -a SWG+Windowed+Banded
```






