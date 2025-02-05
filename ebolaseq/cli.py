import argparse
from ebolaseq.core import main

def cli_main():
    """Entry point for the command line interface."""
    parser = argparse.ArgumentParser(
        description='Download and process Ebola virus sequences from GenBank.',
        usage='%(prog)s -output_loc <output_directory>'
    )
    
    parser.add_argument('-output_loc', '--output_location', 
                        type=str,
                        required=True,
                        help='Output directory for results (use "." for current directory)')
    
    # If no arguments are provided, print help and exit
    import sys
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
        
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    import os
    os.makedirs(args.output_location, exist_ok=True)
    
    # Pass the output location to main
    main(output_location=args.output_location)

if __name__ == "__main__":
    cli_main()
