# EbolaSeq

EbolaSeq is a command-line tool that simplifies the process of analyzing Ebola virus sequences. It automates the complete workflow from downloading sequences to creating phylogenetic trees. The tool retrieves Ebola virus sequences from NCBI GenBank, processes them according to user specifications, performs multiple sequence alignment, and generates phylogenetic trees. Users can choose from five different Ebola virus species and have the option to include their own consensus sequences or exclude specific sequences from the analysis.

## Installation

1. First, install conda if you haven't already:
   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```

2. Create and activate a new conda environment:
   ```bash
   conda create -n ebola_env python=3.9
   conda activate ebola_env
   ```

3. Install required tools in the environment:
   ```bash
   conda install -c bioconda mafft trimal iqtree
   ```

4. Install ebolaseq:
   ```bash
   git clone https://github.com/DaanJansen94/ebolaseq.git
   cd ebolaseq
   pip install .
   ```

## Re-installation

When updates are pushed to GitHub, or when you want to use your own modifications to the code, you'll need to reinstall the package:

```bash
conda activate ebola_env  # Make sure you're in the right environment
cd ebolaseq
git pull  # Get the latest updates from GitHub
pip uninstall ebolaseq
pip install .
```

Note: Any time you modify the code or pull updates from GitHub, you need to reinstall the package using these commands for the changes to take effect.

## Usage

First, make sure your conda environment is activated:
```bash
conda activate ebola_env
```

### Basic Command Structure

```bash
ebolaseq --output-dir <output_directory> [optional arguments]
```

### Required Arguments

- `--output-dir`: Directory where all output files will be saved

### Optional Arguments

- `--consensus-file`: Path to a consensus FASTA file to include in the analysis
- `--remove`: Path to a text file containing GenBank accession numbers to exclude (one per line). Check 'remove.txt' in the repository for an example format.
- `--phylogeny`: Flag to perform phylogenetic analysis (MAFFT alignment, TrimAl trimming, and IQTree2 tree construction)

### Example Commands

1. Basic usage (only download and filter sequences):
   ```bash
   ebolaseq --output-dir my_analysis
   ```

2. Full analysis with phylogenetic tree:
   ```bash
   ebolaseq --output-dir my_analysis --phylogeny
   ```

3. Complete analysis with sequence removal and consensus:
   ```bash
   ebolaseq --output-dir my_analysis \
            --consensus-file path/to/consensus.fasta \
            --remove remove.txt \
            --phylogeny
   ```

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

### Interactive Steps

When running EbolaSeq, you will be prompted to make the following choices:

1. Select an Ebola virus species:
   ```
   Available Ebola virus species:
   1. Zaire ebolavirus
   2. Sudan ebolavirus
   3. Bundibugyo ebolavirus
   4. Tai Forest ebolavirus
   5. Reston ebolavirus
   ```

2. Choose genome completeness:
   ```
   Genome completeness options:
   1. Complete genomes only
   2. Partial genomes (specify minimum completeness)
   3. All genomes (both complete and partial)
   ```
   If option 2 is selected, you'll be asked to enter a minimum completeness percentage (1-100)

3. Select host:
   ```
   Host options:
   1. Human (Homo sapiens)
   2. Chimpanzee (Pan troglodytes)
   3. Gorilla (Gorilla sp.)
   4. All hosts
   ```

4. Choose metadata filter:
   ```
   Metadata filter options:
   1. Location data only
   2. Collection date only
   3. Both location and date
   4. All sequences (no metadata filter)
   ```

5. BEAST format option:
   ```
   Do you want to generate BEAST input format?
   1. No
   2. Yes
   ```

The tool will then proceed with downloading and processing the sequences based on your choices.

### Output Files and Visualization

After running ebolaseq, you'll find the following key output files:

1. Sequence files:
   - `<output_dir>/FASTA/filtered_*_complete_1.fasta`: Filtered sequences
   - `<output_dir>/FASTA/MAFFT_output/Ebola_Combined_aligned.fasta`: MAFFT alignment
   - `<output_dir>/FASTA/TrimAl_output/Ebola_trimmed.fasta`: Trimmed alignment

2. Tree files (in <output_dir>/FASTA/IQTree_output/):
   - `Ebola_tree.treefile`: Final phylogenetic tree (Newick format)
   - `Ebola_tree.iqtree`: IQTree log file with analysis details
   - `Ebola_tree.contree`: Consensus tree with bootstrap support values

3. Metadata:
   - `<output_dir>/FASTA/location.txt`: Tab-separated file with sequence locations (for FigTree)
   - `<output_dir>/BEAST_input/`: Directory containing BEAST-formatted files (optional, created when selecting BEAST output format)
     - `location.txt`: Location data in BEAST format
     - Additional BEAST input files for phylogeographic analysis

### Visualizing the Tree

You can visualize the phylogenetic tree using FigTree (available at: http://tree.bio.ed.ac.uk/software/figtree/):

1. Open FigTree and load `Ebola_tree.treefile` from the IQTree_output directory
2. To color branches by location:
   - Click on 'File' → 'Import Annotations'
   - Select either `location.txt` file (from FASTA or BEAST_input directory)
   - In the sidebar, under 'Appearance', select 'Colour by' → 'location'

This will create a colored tree where branches are grouped by geographical location, helping to visualize the spatial distribution of the virus lineages.

Note: BEAST-formatted output files are optionally created when selecting the BEAST format during the interactive prompts. These files can be used for more detailed phylogeographic analyses in BEAST.

## Dependencies

- Python ≥ 3.6
- BioPython ≥ 1.79
- NumPy ≥ 1.21.0
- MAFFT
- TrimAl
- IQTree2

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Support

If you encounter any problems or have questions, please open an issue on GitHub.
