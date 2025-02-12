import argparse
import os

def cli_main():
    """Entry point for command line interface"""
    parser = argparse.ArgumentParser(description='Download and analyze Ebola virus sequences')
    
    # Required argument
    parser.add_argument('--output-dir', type=str, required=True, 
                       help='Output directory for results')
    
    # Optional arguments
    parser.add_argument('--consensus-file', type=str, 
                       help='Path to consensus FASTA file to include')
    parser.add_argument('--remove', type=str,
                       help='Path to text file containing headers/accession IDs to remove')
    parser.add_argument('--phylogeny', action='store_true', 
                       help='Create phylogenetic tree using IQTree2')
    
    # Non-interactive mode arguments
    parser.add_argument('--virus', type=str, choices=['1', '2', '3', '4', '5'],
                       help='Virus choice: 1=Zaire, 2=Sudan, 3=Bundibugyo, 4=Tai Forest, 5=Reston')
    parser.add_argument('--genome', type=str, choices=['1', '2', '3'],
                       help='Genome type: 1=Complete, 2=Partial, 3=All')
    parser.add_argument('--completeness', type=float,
                       help='Minimum completeness percentage (1-100) when using --genome 2')
    parser.add_argument('--host', type=str, choices=['1', '2', '3'],
                       help='Host: 1=Human, 2=Non-human, 3=All')
    parser.add_argument('--metadata', type=str, choices=['1', '2', '3', '4'],
                       help='Metadata filter: 1=Location, 2=Date, 3=Both, 4=None')
    parser.add_argument('--beast', type=str, choices=['1', '2'],
                       help='BEAST format: 1=No, 2=Yes')
    
    args = parser.parse_args()
    
    # Convert paths to absolute paths
    args.output_dir = os.path.abspath(args.output_dir)
    if args.consensus_file:
        args.consensus_file = os.path.abspath(args.consensus_file)
    if args.remove:
        args.remove = os.path.abspath(args.remove)
    
    # Create output directory
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)
    
    # Change to output directory
    original_dir = os.getcwd()
    os.chdir(args.output_dir)
    
    try:
        # Check if all required non-interactive parameters are provided
        non_interactive = all([args.virus, args.genome, args.host, args.metadata])
        if args.genome == '2' and args.completeness is None:
            non_interactive = False
        if args.metadata in ['2', '3'] and args.beast is None:
            non_interactive = False
            
        # Run in appropriate mode
        if non_interactive:
            main(args, non_interactive=True)
        else:
            main(args)
    finally:
        os.chdir(original_dir)

def main(args, non_interactive=False):
    """Main function that can run in both interactive and non-interactive modes"""
    # Make these dictionaries global so they're accessible
    global virus_options, virus_processing_names
    
    # ... (existing dictionary definitions) ...
    
    if non_interactive:
        # Use command-line arguments directly
        choice = args.virus
        virus_display_name = virus_options[choice]
        virus_choice = virus_processing_names[choice]
        virus_type = virus_choice.split('_')[0]
        genome_choice = args.genome
        completeness_threshold = args.completeness if args.genome == '2' else 0
        host_choice = args.host
        metadata_choice = args.metadata
        beast_choice = args.beast if args.metadata in ['2', '3'] else '1'
    else:
        # Interactive mode (existing code)
        print("\nAvailable Ebola virus species:")
        for key, value in virus_options.items():
            print(f"{key}. {value}")
        choice = input("\nSelect the virus type (1-5): ")
        # ... (rest of existing interactive code) ... 