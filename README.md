# EbolaSeq

EbolaSeq is a command-line tool that simplifies the process of analyzing Ebola virus sequences. It automates the complete workflow from downloading sequences to creating phylogenetic trees. The tool retrieves Ebola virus sequences from NCBI GenBank, processes them according to user specifications, performs multiple sequence alignment and generates phylogenetic trees.

## Installation

### Prerequisites
First, install conda if you haven't already:
```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

Then, ensure you have the required channels:
```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

### Option 1: Using Conda (Recommended)
Install [EbolaSeq via Conda](https://anaconda.org/bioconda/ebolaseq):
```bash
conda create -n ebolaseq -c conda-forge -c bioconda ebolaseq -y
conda activate ebolaseq
```

### Option 2: From Source Code
1. Create and activate a new conda environment:
   ```bash
   conda create -n ebolaseq -c conda-forge -c bioconda python=3.9 mafft trimal iqtree=2.4.0 biopython minimap2 pal2nal
   conda activate ebolaseq
   ```

2. Install ebolaseq:
   ```bash
   git clone https://github.com/DaanJansen94/ebolaseq.git
   cd ebolaseq
   pip install .
   ```

3. Re-installation (when updates are available):
   ```bash
   conda activate ebolaseq  # Make sure you're in the right environment
   cd ebolaseq
   git pull  # Get the latest updates from GitHub
   pip uninstall ebolaseq
   pip install .
   ```
   Note: Any time you modify the code or pull updates from GitHub, you need to reinstall the package using these commands for the changes to take effect.

## Usage

First, make sure your conda environment is activated:
```bash
conda activate ebolaseq
```

EbolaSeq can be run in two modes:
1. Interactive Mode (default) - for interactive use
2. Non-Interactive Mode - for HPC submissions or automated runs

### Interactive Mode

Basic command structure:
```bash
ebolaseq --output-dir my_analysis [optional arguments]
```

Example commands:
```bash
# Basic usage
ebolaseq --output-dir my_analysis

# With phylogenetic analysis
ebolaseq --output-dir my_analysis --phylogeny

# Complete analysis with sequence removal and consensus
ebolaseq --output-dir my_analysis \
         --consensus-file path/to/consensus.fasta \
         --remove remove.txt \
         --phylogeny

# Complete analysis with sequence removal and consensus and phylogenetic analysis 
ebolaseq -o my_analysis -c path/to/consensus.fasta -p
```

The tool will interactively prompt you for choices about:
- Virus species
- Genome completeness
- Host type
- Metadata options

### Non-Interactive Mode (HPC)

For HPC submissions or automated runs, specify all parameters via command line:

```bash
ebolaseq --output-dir my_analysis \
         --virus 1 \
         --genome 2 \
         --completeness 80 \
         --host 1 \
         --metadata 3 \
         --beast 2 \
         --phylogeny \
         --consensus-file path/to/consensus.fasta \
         --remove remove.txt
```

Required parameters for non-interactive mode:
- `--virus`: Virus type
  - 1 = Zaire ebolavirus
  - 2 = Sudan ebolavirus
  - 3 = Bundibugyo ebolavirus
  - 4 = Tai Forest ebolavirus
  - 5 = Reston ebolavirus

- `--genome`: Genome completeness
  - 1 = Complete genomes only
  - 2 = Partial genomes (requires --completeness)
  - 3 = All genomes

- `--completeness`: Required when --genome=2
  - Value between 1-100 (percentage)

- `--host`: Host filter
  - 1 = Human only
  - 2 = Non-human only
  - 3 = All hosts

- `--metadata`: Metadata filter
  - 1 = Location only
  - 2 = Date only
  - 3 = Both location and date
  - 4 = None

- `--beast`: Required when --metadata is 2 or 3
  - 1 = No
  - 2 = Yes

Optional arguments (both modes):
- `--output-dir`: Output directory (required)
- `--consensus-file`: Path to consensus FASTA file
- `--remove`: Path to sequence removal list
- `--phylogeny`: Create phylogenetic tree

### Input File Formats

1. Consensus file (optional):
   - FASTA format
   - Example:
   ```
   >Consensus_sequence_name
   ATGCATGCATGC...
   ```

2. Remove file (optional):
   - Text file with one GenBank accession number per line
   - Example:
   ```
   KM034562.1
   KM034563.1
   MK114118.1
   ```

## Important Notes

It is strongly recommended to always use a `remove.txt` file with the `--remove` option when running EbolaSeq. This file should contain sequence IDs that should be excluded from the analysis, particularly sequences obtained from cell culture passages, laboratory-adapted strains, artificially modified sequences, sequences from experimental infections, and other non-natural viral sequences. These sequences can bias analyses as they may not represent natural viral diversity.

When creating large phylogenetic trees for Zaire ebolavirus, it is recommended to root the tree using sequences from the 1976 Yambuku outbreak, as this represents the first documented outbreak of the virus.

## Dependencies

- Python ≥ 3.6
- BioPython ≥ 1.79
- NumPy ≥ 1.21.0
- MAFFT
- TrimAl
- IQTree2

## Citation

If you use EbolaSeq in your research, please cite:

```
Jansen, D., & Vercauteren, K. (2025). EbolaSeq: A Command-Line Tool for Downloading, Processing, and Analyzing Ebola Virus Sequences for Phylogenetic Analysis (v0.1.1). Zenodo. https://doi.org/10.5281/zenodo.14851686
```

## License

This project is licensed under the GNU General Public License v3.0 (GPL-3.0) - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

If you encounter any problems or have questions, please open an issue on GitHub.
