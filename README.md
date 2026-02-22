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

### Option 2: From source
```bash
conda create -n ebolaseq -c conda-forge -c bioconda python=3.9 mafft trimal iqtree=2.4.0 biopython minimap2 pal2nal
conda activate ebolaseq
git clone https://github.com/DaanJansen94/ebolaseq.git
cd ebolaseq
pip install .
```

## Usage

```bash
conda activate ebolaseq
ebolaseq -o OUTPUT_DIR [options] 
```

EbolaSeq can be run in two modes:
1. Interactive Mode (default) - for interactive use
2. Non-Interactive Mode - for HPC submissions or automated runs

### Options (submission reference)

**Required**

**`-o`, `--output-dir`** — Output directory for results

**`--virus`** — Virus / species

1 = Zaire ebolavirus  
2 = Sudan ebolavirus  
3 = Bundibugyo ebolavirus  
4 = Tai Forest ebolavirus  
5 = Reston ebolavirus  
6 = Pan-Ebola: all 5 species  
Comma-separated = multiple species (e.g. 1,2 for Zaire+Sudan; 1,2,3 for Zaire+Sudan+Bundibugyo)

**`--genome`** — Genome completeness

1 = Complete genomes only  
2 = Partial genomes (requires `--completeness`)  
3 = All genomes

**`--completeness`** — Required when `--genome`=2  

Value between 1–100 (percentage)

**`--host`** — Host filter

1 = Human only  
2 = Non-human only  
3 = All hosts

**`--metadata`** — Metadata filter

1 = Location only  
2 = Date only  
3 = Both location and date  
4 = None

**`--beast`** — Required when `--metadata` is 2 or 3

1 = No  
2 = Yes

**Consensus FASTA per species** — Path to a FASTA file; used in alignment/phylogeny when provided

`--c_z` = Zaire  
`--c_s` = Sudan  
`--c_r` = Reston  
`--c_b` = Bundibugyo  
`--c_t` = Tai Forest

**`--alignment`, `-a`** — Alignment type

1 = Whole-genome alignment  
2 = Protein (CDS) alignment  
3 = No alignment

**`--proteins`, `-pr`** — For alignment 2 only; comma-separated

1 = L  
2 = NP  
3 = VP35  
4 = VP40  
5 = VP30  
6 = VP24  
(Or use names: L, NP, VP35, VP40, VP30, VP24)

**Optional (both modes)**

**`--remove`** — Path to file listing sequence IDs/headers to exclude  
**`--phylogeny`, `-p`** — Create phylogenetic tree from alignment

### Examples

```bash
# Interactive (prompts for all choices)
ebolaseq -o my_analysis

# Non-interactive: Zaire, complete genomes, human, location+date, whole-genome alignment + phylogeny
ebolaseq -o my_analysis --virus 1 --genome 1 --host 1 --metadata 3 --alignment 1 --phylogeny

# Pan-Ebola, protein alignment (L and NP), phylogeny per protein, consensus for Zaire and Sudan
ebolaseq -o my_analysis --virus 6 --genome 1 --host 3 --metadata 4 \
  --c_z consensus_zaire.fasta --c_s consensus_sudan.fasta \
  --alignment 2 -pr L,NP --phylogeny

# Exclude specific sequences
ebolaseq -o my_analysis --virus 1 --genome 1 --host 1 --metadata 4 --remove exclude.txt
```

## Output 

- **FASTA/** — Filtered sequences and `location.txt`.
- **Alignment/** — For whole-genome: `FASTA/`, `MAFFT/`, `Trimmed/`. For protein: `pan/` (or species name) with e.g. `L/`, `NP/` each containing `cds_aligned.fasta`.
- **Phylogeny/** — IQTree2 results (whole-genome: one tree; protein: one folder per protein).
- **summary_*.txt** — Run summary and location counts.

## Notes

- Use `--remove` with a list of IDs to exclude cell-culture, lab-adapted, or other non-natural sequences.
- For large Zaire trees, consider rooting with 1976 Yambuku outbreak sequences.

## Dependencies

- Python ≥ 3.9
- Biopython ≥ 1.81
- MAFFT, TrimAl, IQTree2  
- For protein alignment: minimap2, pal2nal

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