# EbolaSeq

A command-line tool that downloads and filters Ebola virus sequences from GenBank for phylogenetic analysis. Automatically handles sequence filtering, metadata formatting, and BEAST input preparation.

## Installation

Install using pip:
```
pip install ebolaseq
```

## Usage

Simply run:
```
ebolaseq
```

Follow the prompts to select:
1. Virus species (Zaire, Sudan, Bundibugyo, Tai Forest, or Reston)
2. Genome completeness (complete, partial, or all)
3. Host species (Human, Chimpanzee, Gorilla, or all)
4. Metadata filters (location, date, or both)
5. BEAST output option (yes/no)

## Output Structure

The tool organizes outputs in the following structure:

1. **FASTA Directory** (`FASTA/`)
   - `filtered_[virus]_[completeness]_[host].fasta`: Main sequence file
   - `outgroup_[virus]_[completeness]_[host].fasta`: Outgroup sequence

2. **BEAST Files** (`BEAST_input/`) - *Optional*
   - BEAST-formatted sequence files
   - Location data for phylogeographic analysis

3. **Summary File** (`summary_[virus]_[completeness]_[host].txt`)
   - Run parameters
   - Sequence counts
   - Location statistics

## Requirements
- Python ≥3.6
- Biopython (installed automatically)

## Source Code
[https://github.com/DaanJansen94/ebolaseq](https://github.com/DaanJansen94/ebolaseq)

## License
MIT License

## Issues
Report bugs on the [GitHub issues page](https://github.com/yourusername/ebolaseq/issues)
