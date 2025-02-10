# EbolaSeq
A command-line tool that downloads and filters Ebola virus sequences from GenBank for phylogenetic analysis. Automatically handles sequence filtering, metadata formatting, and BEAST input preparation.

## Installation

Before installing from source, ensure you have:
- Python 3.6 or higher
- Biopython package

```bash
git clone https://github.com/DaanJansen94/ebolaseq
cd ebolaseq
pip install .
```

## Re-installation

To reinstall the package with latest changes:

```bash
cd ebolaseq
pip uninstall ebolaseq
pip install .
```

## Usage

### Basic Command
```bash
ebolaseq -output_loc <output_directory>
```

### Interactive Options
1. Virus Species Selection:
   - 1: Zaire ebolavirus
   - 2: Sudan ebolavirus
   - 3: Bundibugyo ebolavirus
   - 4: Tai Forest ebolavirus
   - 5: Reston ebolavirus

2. Genome Completeness:
   - 1: Complete genomes only
   - 2: Partial genomes (specify minimum completeness)
   - 3: All genomes

3. Host Selection:
   - 1: Human (Homo sapiens)
   - 2: Chimpanzee (Pan troglodytes)
   - 3: Gorilla (Gorilla sp.)
   - 4: All hosts

4. Metadata Filters:
   - 1: Location data only
   - 2: Collection date only
   - 3: Both location and date
   - 4: All sequences

5. BEAST Output:
   - 1: Standard FASTA only
   - 2: Generate BEAST input

## Output Files

The tool creates:
- `FASTA/filtered_[virus]_[completeness]_[host].fasta`: Main sequence file
- `FASTA/outgroup_[virus]_[completeness]_[host].fasta`: Outgroup sequence
- `BEAST_input/`: BEAST-formatted files (optional)
- `summary_[virus]_[completeness]_[host].txt`: Run parameters and statistics

## License
MIT License
