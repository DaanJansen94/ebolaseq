import argparse
from ebolaseq.core import main

def cli_main():
    """Entry point for the command line interface."""
    parser = argparse.ArgumentParser(description='Process Ebola virus sequences.')
    parser.add_argument('-output_loc', '--output_location', 
                        type=str,
                        required=True,
                        help='Output directory for results (use "." for current directory)')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    import os
    os.makedirs(args.output_location, exist_ok=True)
    
    # Pass the output location to main
    main(output_location=args.output_location)

if __name__ == "__main__":
    cli_main()
