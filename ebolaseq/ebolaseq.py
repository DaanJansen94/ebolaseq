from Bio import Entrez
from Bio import SeqIO
import sys
import subprocess
import os
import re
import datetime
import time
import argparse
import shutil
import tempfile
from collections import OrderedDict

# Import version from package
try:
    from . import __version__
except ImportError:
    # Fallback for development/testing
    try:
        import ebolaseq
        __version__ = ebolaseq.__version__
    except:
        __version__ = "unknown"

# Copy ALL your original functions here exactly as they were
def get_sequence_length(record):
    return len(record.seq)

def standardize_country(location):
    # Dictionary mapping provinces/regions to their countries
    province_to_country = {
        # Sierra Leone regions
        'western_rural': 'Sierra_Leone',
        'western_urban': 'Sierra_Leone',
        'western_area': 'Sierra_Leone',
        'Bombali': 'Sierra_Leone',
        'Port_Loko': 'Sierra_Leone',
        'Kambia': 'Sierra_Leone',
        'Kailahun': 'Sierra_Leone',
        'Kenema': 'Sierra_Leone',
        'Tonkolili': 'Sierra_Leone',
        # Liberia regions
        'Montserrado': 'Liberia',
        'Margibi': 'Liberia',
        'Sinoe': 'Liberia',
        'Bong': 'Liberia',
        'Grand_Bassa': 'Liberia',
        'Grand_Kru': 'Liberia',
        # Guinea regions
        'Gueckedou': 'Guinea',
        'Mandiana': 'Guinea',
        'Meliandou': 'Guinea',
        'Kissidougou': 'Guinea',
        'Macenta': 'Guinea',
        'Nzerekore': 'Guinea',
        'Conakry': 'Guinea',
        'Forecariah': 'Guinea',
        # DRC regions
        'Makanza': 'Democratic_Republic_of_the_Congo',
        'Kikwit': 'Democratic_Republic_of_the_Congo',
        'Kinshasa': 'Democratic_Republic_of_the_Congo',
        'Kivu': 'Democratic_Republic_of_the_Congo',
        'Equateur': 'Democratic_Republic_of_the_Congo'
    }
    
    # Dictionary for standardizing country names
    country_standardization = {
        'Zaire': 'Democratic_Republic_of_the_Congo',
        'DRC': 'Democratic_Republic_of_the_Congo',
        'Guinea': 'Republic_of_Guinea',
        'Great_Britain': 'United_Kingdom',
        'UK': 'United_Kingdom',
        'England': 'United_Kingdom',
        'Britain': 'United_Kingdom',
        "Cote_d'Ivoire": 'Ivory_Coast',
        "Côte_d'Ivoire": 'Ivory_Coast'
    }
    
    # Clean up location string
    location = location.replace(" ", "_")
    location = location.split(',')[0].split(':')[0].strip()
    
    # Check if it's a province we know
    if location in province_to_country:
        return province_to_country[location]
        
    # Check if it's a country name that needs standardization
    if location in country_standardization:
        return country_standardization[location]
        
    # Handle case where full DRC name is in string
    if 'Democratic_Republic_of_the_Congo' in location:
        return 'Democratic_Republic_of_the_Congo'
    
    return location

def get_location(record):
    for feature in record.features:
        if feature.type == "source":
            if 'country' in feature.qualifiers:
                return standardize_country(feature.qualifiers['country'][0])
            elif 'geo_loc_name' in feature.qualifiers:
                return standardize_country(feature.qualifiers['geo_loc_name'][0])
    return "Unknown"

def get_collection_date(record):
    for feature in record.features:
        if feature.type == "source":
            if 'collection_date' in feature.qualifiers:
                date_str = feature.qualifiers['collection_date'][0]
                
                # Skip if date is obviously invalid
                if date_str in ['', 'unknown', 'Unknown', 'NA', 'missing']:
                    return 9999
                
                try:
                    # Handle various date formats
                    date_str = date_str.replace('_', '-').replace('/', '-')
                    
                    # Try to extract year
                    if '-' in date_str:
                        parts = date_str.split('-')
                        # Check if year is in the first or last position
                        potential_year = max(int(parts[0]), int(parts[-1]))
                        # Validate year is reasonable (between 1900 and current year)
                        if 1900 <= potential_year <= 2024:
                            return potential_year
                    else:
                        # Try direct conversion
                        year = int(date_str)
                        if 1900 <= year <= 2024:
                            return year
                        
                except (ValueError, IndexError):
                    # If we can't parse the date, look for a 4-digit year in the string
                    year_match = re.search(r'(19|20)\d{2}', date_str)
                    if year_match:
                        return int(year_match.group())
                    
    return 9999  # Return large number for unknown dates

# ebolaS species names (for protein pipeline)
PROCESSING_TO_EBOLAS_SPECIES = {
    'Zaire_ebolavirus': 'zaire',
    'Sudan_ebolavirus': 'sudan',
    'Reston_ebolavirus': 'reston',
    'Bundibugyo_ebolavirus': 'bundibugyo',
    'Tai_Forest_ebolavirus': 'tai_forest',
}
# NCBI / GenBank sometimes use virus name or variant in the ID; map to same ebolaS species
SPECIES_ID_ALIASES = {
    # Virus name (no "ebolavirus")
    'Bundibugyo_virus': 'bundibugyo',
    'Zaire_virus': 'zaire',
    'Sudan_virus': 'sudan',
    'Reston_virus': 'reston',
    'Tai_Forest_virus': 'tai_forest',
    # Sudan variants
    'Sudan_ebolavirus_-_Nakisamata': 'sudan',
    'Sudan_ebolavirus_-_Sudan': 'sudan',
    'Sudan_ebolavirus_-_Gulu': 'sudan',
    # Reston variants
    'Reston_ebolavirus_-_Reston': 'reston',
    'Reston_ebolavirus_-_Pennsylvania': 'reston',
    # Zaire variants
    'Zaire_ebolavirus_-_Zaire': 'zaire',
    'Zaire_ebolavirus_-_Mayinga': 'zaire',
    'Zaire_ebolavirus_-_Ecran': 'zaire',
    'Zaire_ebolavirus_-_Kikwit': 'zaire',
    'Zaire_ebolavirus_-_Makona': 'zaire',
    # Bundibugyo (already have Bundibugyo_virus)
    'Bundibugyo_ebolavirus_-_Bundibugyo': 'bundibugyo',
    # Tai Forest / Côte d'Ivoire / Ivory Coast
    'Tai_Forest_ebolavirus_-_Tai_Forest': 'tai_forest',
    'Tai_Forest_ebolavirus_-_Cote_d_Ivoire': 'tai_forest',
    'Cote_d_Ivoire_ebolavirus': 'tai_forest',
    'Ivory_Coast_ebolavirus': 'tai_forest',
}

# CDS pipeline: RefSeq and coordinates per species (protein alignment)
REFSEQ_IDS_CDS = {
    "zaire": "NC_002549.1", "sudan": "NC_006432.1", "reston": "NC_004161.1",
    "bundibugyo": "NC_014373.1", "tai_forest": "NC_014372.1",
}
SPECIES_REGIONS_CDS = {
    "zaire": {"L": (11581, 18219), "NP": (470, 2689), "VP35": (3129, 4151), "VP40": (4479, 5459), "VP30": (8509, 9375), "VP24": (10345, 11100)},
    "sudan": {"L": (11535, 18167), "NP": (458, 2674), "VP35": (3138, 4127), "VP40": (4454, 5434), "VP30": (8441, 9307), "VP24": (10299, 11054)},
    "reston": {"L": (11550, 18188), "NP": (464, 2683), "VP35": (3155, 4144), "VP40": (4485, 5480), "VP30": (8490, 9353), "VP24": (10303, 11058)},
    "bundibugyo": {"L": (11567, 18199), "NP": (458, 2677), "VP35": (3108, 4133), "VP40": (4461, 5441), "VP30": (8496, 9365), "VP24": (10335, 11090)},
    "tai_forest": {"L": (11566, 18198), "NP": (464, 2683), "VP35": (3114, 4139), "VP40": (4467, 5447), "VP30": (8503, 9372), "VP24": (10339, 11094)},
}
PROTEIN_DESCRIPTIONS_CDS = {
    "L": "RNA-dependent RNA polymerase", "NP": "nucleoprotein", "VP35": "polymerase cofactor",
    "VP40": "matrix protein", "VP30": "minor nucleoprotein", "VP24": "membrane-associated protein",
}

# Map GenBank organism strings to internal processing names (with underscores)
# NCBI SOURCE/ORGANISM can be "Bundibugyo virus" or "Bundibugyo ebolavirus"; include both
ORGANISM_TO_PROCESSING = {
    'Zaire ebolavirus': 'Zaire_ebolavirus',
    'Zaire Ebolavirus': 'Zaire_ebolavirus',
    'Zaire virus': 'Zaire_ebolavirus',
    'Sudan ebolavirus': 'Sudan_ebolavirus',
    'Sudan Ebolavirus': 'Sudan_ebolavirus',
    'Sudan virus': 'Sudan_ebolavirus',
    'Bundibugyo ebolavirus': 'Bundibugyo_ebolavirus',
    'Bundibugyo virus': 'Bundibugyo_ebolavirus',
    'Tai Forest ebolavirus': 'Tai_Forest_ebolavirus',
    'Tai Forest virus': 'Tai_Forest_ebolavirus',
    "Cote d'Ivoire ebolavirus": 'Tai_Forest_ebolavirus',
    'Reston ebolavirus': 'Reston_ebolavirus',
    'Reston Ebolavirus': 'Reston_ebolavirus',
    'Reston virus': 'Reston_ebolavirus',
}

# Cache for accession -> ebolaS species (avoids repeated NCBI lookups)
_accession_species_cache = {}

def get_ebolas_species_for_accession(accession):
    """Fetch ORGANISM from NCBI for this accession and return ebolaS species (zaire, sudan, etc.) or None. Cached."""
    accession = (accession or "").strip().split(".")[0]  # MT680258.1 -> MT680258
    if not accession:
        return None
    if accession in _accession_species_cache:
        return _accession_species_cache[accession]
    try:
        Entrez.email = "anonymous@example.com"
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
        text = handle.read()
        handle.close()
        if isinstance(text, bytes):
            text = text.decode("utf-8", errors="replace")
        # Parse ORGANISM line: "  ORGANISM  Bundibugyo virus"
        organism = None
        for line in text.splitlines():
            if line.startswith("  ORGANISM"):
                organism = line.replace("  ORGANISM", "").strip().strip(".")
                break
        if not organism:
            _accession_species_cache[accession] = None
            return None
        proc_name = ORGANISM_TO_PROCESSING.get(organism)
        if proc_name:
            ebolaS = PROCESSING_TO_EBOLAS_SPECIES.get(proc_name)
            _accession_species_cache[accession] = ebolaS
            return ebolaS
        _accession_species_cache[accession] = None
        return None
    except Exception:
        _accession_species_cache[accession] = None
        return None

def get_record_species_name(record):
    """Get the virus species processing name from a GenBank record (e.g. Zaire_ebolavirus)."""
    for feature in record.features:
        if feature.type == "source" and 'organism' in feature.qualifiers:
            org = feature.qualifiers['organism'][0]
            return ORGANISM_TO_PROCESSING.get(org, org.replace(' ', '_'))
    return None

def parse_virus_selection(choice, virus_processing_names):
    """
    Parse user virus selection. choice can be:
    - '1' to '5': single species
    - '6': all 5 species (Pan-Ebola)
    - '1,2' or '1, 2, 3': comma-separated species numbers
    Returns (list of processing names, base name for filenames).
    """
    choice = choice.strip().lower()
    all_species = [virus_processing_names[k] for k in sorted(virus_processing_names)]
    # Option 6: all 5 species (Pan-Ebola)
    if choice == '6':
        base_name = 'pan'
        return all_species, base_name
    # Comma-separated numbers
    if ',' in choice:
        parts = [p.strip() for p in choice.split(',') if p.strip()]
        if not parts:
            return None, None
        valid = []
        for p in parts:
            if p in virus_processing_names:
                valid.append(virus_processing_names[p])
            else:
                return None, None
        # unique, preserve order
        seen = set()
        ordered = []
        for v in valid:
            if v not in seen:
                seen.add(v)
                ordered.append(v)
        # base name: zaire_sudan, zaire_sudan_bundibugyo, etc.
        ref_key = {'Zaire_ebolavirus': 'zaire', 'Sudan_ebolavirus': 'sudan', 'Bundibugyo_ebolavirus': 'bundibugyo',
                   'Tai_Forest_ebolavirus': 'tai_forest', 'Reston_ebolavirus': 'reston'}
        base_name = '_'.join(ref_key[s] for s in ordered)
        return ordered, base_name
    # Single digit 1-5
    if choice in virus_processing_names:
        name = virus_processing_names[choice]
        return [name], name.lower()
    return None, None

def get_outgroup_reference(virus_choice_or_list):
    """Get outgroup reference. virus_choice_or_list: single species name or list of species (for pan/multi)."""
    Entrez.email = "anonymous@example.com"
    virus_list = virus_choice_or_list if isinstance(virus_choice_or_list, list) else [virus_choice_or_list]
    virus_set = set(virus_list)

    # Dictionary of reference sequences for each species (RefSeq accessions)
    refseq_dict = {
        'Zaire_ebolavirus': {'outgroup': 'Sudan_ebolavirus', 'refseq': 'NC_006432'},
        'Sudan_ebolavirus': {'outgroup': 'Zaire_ebolavirus', 'refseq': 'NC_002549'},
        'Bundibugyo_ebolavirus': {'outgroup': 'Zaire_ebolavirus', 'refseq': 'NC_002549'},
        'Tai_Forest_ebolavirus': {'outgroup': 'Zaire_ebolavirus', 'refseq': 'NC_002549'},
        'Reston_ebolavirus': {'outgroup': 'Zaire_ebolavirus', 'refseq': 'NC_002549'}
    }
    # For pan (all 5), use Bombali ebolavirus as outgroup (not in the 5 species)
    bombali_refseq = 'NC_039345'

    if len(virus_set) == 1:
        virus_choice = virus_list[0]
        outgroup_info = refseq_dict[virus_choice]
        refseq_id = outgroup_info['refseq']
        outgroup_species = outgroup_info['outgroup']
    elif len(virus_set) >= 2:
        # Pan or any multi-species: use Bombali ebolavirus as sole outgroup
        refseq_id = bombali_refseq
        outgroup_species = 'Bombali_ebolavirus'
    else:
        # Multi (2–4 species): use first species not in the set
        for sp in refseq_dict:
            if sp not in virus_set:
                outgroup_info = refseq_dict[sp]
                refseq_id = outgroup_info['refseq']
                outgroup_species = outgroup_info['outgroup']
                break
        else:
            refseq_id = refseq_dict['Zaire_ebolavirus']['refseq']
            outgroup_species = 'Sudan_ebolavirus'

    handle = Entrez.efetch(db="nucleotide", id=refseq_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()

    return record, outgroup_species

def get_sequences_to_remove(remove_file):
    """Read list of sequences to remove from file."""
    # Get absolute path
    remove_file = os.path.abspath(remove_file)
    
    if not os.path.exists(remove_file):
        print(f"Warning: Remove file {remove_file} not found")
        return set()
    
    to_remove = set()
    try:
        with open(remove_file) as f:
            for line in f:
                # Strip whitespace and ignore empty lines or comments
                line = line.strip()
                if line and not line.startswith('#'):
                    to_remove.add(line)
        
        if to_remove:
            print(f"Found {len(to_remove)} sequences to remove")
        else:
            print("No sequences to remove found in file")
            
    except Exception as e:
        print(f"Error reading remove file: {str(e)}")
        return set()
    
    return to_remove

def cli_main():
    """Entry point for command line interface"""
    parser = argparse.ArgumentParser(
        description=f'Download and analyze Ebola virus sequences ({__version__})',
        prog='ebolaseq'
    )
    parser.add_argument('--version', action='version', version=f'ebolaseq {__version__}', help="show program's version number")

    req = parser.add_argument_group('Required')
    req.add_argument('-o', '--output-dir', type=str, required=True, metavar='DIR',
                     help='Output directory for all results')

    opt = parser.add_argument_group('Optional — download filters (non-interactive)')
    opt.add_argument('--virus', type=str, metavar='CHOICE',
                     help='Species: 1 Zaire, 2 Sudan, 3 Bundibugyo, 4 Tai Forest, 5 Reston, 6 all five, or e.g. 1,2')
    opt.add_argument('--genome', type=str, choices=['1', '2', '3'], metavar='{1,2,3}',
                     help='1 complete only, 2 partial (use --completeness), 3 all')
    opt.add_argument('--completeness', type=float, metavar='PCT',
                     help='Min completeness percent (1–100) when genome=2')
    opt.add_argument('--host', type=str, choices=['1', '2', '3'], metavar='{1,2,3}',
                     help='1 human, 2 non-human, 3 all hosts')
    opt.add_argument('--metadata', type=str, choices=['1', '2', '3', '4'], metavar='{1,2,3,4}',
                     help='1 location only, 2 date only, 3 both, 4 no filter')
    opt.add_argument('--beast', type=str, choices=['1', '2'], metavar='{1,2}',
                     help='BEAST input: 1 no, 2 yes (only if metadata 2 or 3)')

    opt_cons = parser.add_argument_group('Optional — consensus FASTA per species')
    opt_cons.add_argument('--c_z', type=str, dest='consensus_zaire', metavar='FASTA', help='Zaire ebolavirus')
    opt_cons.add_argument('--c_s', type=str, dest='consensus_sudan', metavar='FASTA', help='Sudan ebolavirus')
    opt_cons.add_argument('--c_r', type=str, dest='consensus_reston', metavar='FASTA', help='Reston ebolavirus')
    opt_cons.add_argument('--c_b', type=str, dest='consensus_bundibugyo', metavar='FASTA', help='Bundibugyo ebolavirus')
    opt_cons.add_argument('--c_t', type=str, dest='consensus_tai', metavar='FASTA', help='Tai Forest ebolavirus')

    opt_align = parser.add_argument_group('Optional — alignment and phylogeny')
    opt_align.add_argument('--alignment', '-a', type=str, choices=['1', '2', '3'], default=None, metavar='{1,2,3}',
                           help='1 whole-genome alignment, 2 protein (CDS) alignment, 3 no alignment')
    opt_align.add_argument('--proteins', '-pr', type=str, default=None, metavar='LIST',
                           help='For alignment=2: comma-separated L,NP,VP35,VP40,VP30,VP24')
    opt_align.add_argument('--phylogeny', '-p', action='store_true', help='Build phylogeny from alignment')
    opt_align.add_argument('-m', '--min-cds-fraction', type=float, default=0.5, metavar='F',
                           help='Min fraction of reference CDS length to keep a sequence (default: 0.5). E.g. 0.2 keeps more partial, 0.8 is stricter.')
    opt_align.add_argument('-t', '--threads', type=int, default=1, metavar='N',
                           help='Threads for minimap2, MAFFT, and IQTree2 (default: 1). E.g. -t 64 on a 64-core node. 0 = use all CPUs.')

    opt_other = parser.add_argument_group('Optional — other')
    opt_other.add_argument('--remove', type=str, metavar='FILE', help='File listing sequence IDs/headers to exclude')

    args = parser.parse_args()

    args.output_dir = os.path.abspath(args.output_dir)
    if args.remove:
        args.remove = os.path.abspath(args.remove)
    for attr in ('consensus_zaire', 'consensus_sudan', 'consensus_reston', 'consensus_bundibugyo', 'consensus_tai'):
        p = getattr(args, attr, None)
        if p:
            setattr(args, attr, os.path.abspath(p))
    # Default BEAST to No when metadata is 2 or 3 but --beast not given (so non-interactive can run without prompting)
    if args.metadata in ('2', '3') and args.beast is None:
        args.beast = '1'

    # Override output directory: remove if exists and create fresh (no leftover from previous run)
    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir)
    os.makedirs(args.output_dir)

    # Change to output directory
    original_dir = os.getcwd()
    os.chdir(args.output_dir)
    
    try:
        # Check if all required non-interactive parameters are provided
        non_interactive = all([args.virus, args.genome, args.host, args.metadata])
        if args.genome == '2' and args.completeness is None:
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
    
    # Available virus choices (store with underscores but display without)
    virus_options = {
        '1': 'Zaire ebolavirus',
        '2': 'Sudan ebolavirus',
        '3': 'Bundibugyo ebolavirus',
        '4': 'Tai Forest ebolavirus',
        '5': 'Reston ebolavirus'
    }
    
    # Internal mapping for processing (with underscores)
    virus_processing_names = {
        '1': 'Zaire_ebolavirus',
        '2': 'Sudan_ebolavirus',
        '3': 'Bundibugyo_ebolavirus',
        '4': 'Tai_Forest_ebolavirus',
        '5': 'Reston_ebolavirus'
    }
    
    # Genome type options
    genome_options = {
        '1': 'Complete genomes only',
        '2': 'Partial genomes (specify minimum completeness)',
        '3': 'All genomes (both complete and partial)'
    }
    
    # Host options
    host_options = {
        '1': 'Human (Homo sapiens)',
        '2': 'Non-human (all animal hosts)',
        '3': 'All hosts'
    }
    
    if non_interactive:
        choice = args.virus.strip().lower()
        virus_choices, base_name = parse_virus_selection(choice, virus_processing_names)
        if not virus_choices:
            print(f"Invalid virus selection: {args.virus}. Use 1-6 or comma-separated (e.g. 1,2).")
            return
        virus_display_name = base_name if len(virus_choices) > 1 else virus_options.get(
            next((k for k, v in virus_processing_names.items() if v == virus_choices[0]), ''), base_name
        )
        genome_choice = args.genome
        completeness_threshold = args.completeness if args.genome == '2' else 0
        host_choice = args.host
        metadata_choice = args.metadata
        beast_choice = args.beast if args.metadata in ['2', '3'] else '1'
        if args.alignment is None:
            args.alignment = '3'
        if args.alignment == '2' and not args.proteins:
            args.proteins = 'L'
    else:
        # Interactive mode
        print("\nAvailable Ebola virus species:")
        for key, value in virus_options.items():
            print(f"{key}. {value}")
        print("6. Pan (all 5 species)")
        print("   Or enter multiple species as comma-separated numbers (e.g. 1,2 or 1,2,3)")
        choice = input("\nSelect the virus type (1-6, or e.g. 1,2): ").strip()
        virus_choices, base_name = parse_virus_selection(choice, virus_processing_names)
        if not virus_choices:
            print("Invalid selection. Use 1-6 (6 = all 5 species) or e.g. 1,2 for multiple.")
            return
        virus_display_name = "Pan (all 5)" if base_name == "pan" else (
            ", ".join(virus_options.get(next((k for k, v in virus_processing_names.items() if v == s), ''), s) for s in virus_choices)
            if len(virus_choices) > 1 else virus_options.get(choice, base_name)
        )
        
        # Show genome completeness options
        print("\nGenome completeness options:")
        for key, value in genome_options.items():
            print(f"{key}. {value}")
        genome_choice = input("\nSelect genome type (1-3): ")
        
        # If partial genomes selected, get threshold
        completeness_threshold = 0
        if genome_choice == '2':
            while True:
                try:
                    completeness_threshold = float(input("\nEnter minimum completeness percentage (1-100): "))
                    if 1 <= completeness_threshold <= 100:
                        break
                    print("Please enter a value between 1 and 100")
                except ValueError:
                    print("Please enter a valid number")
        
        # Show host options
        print("\nHost options:")
        for key, value in host_options.items():
            print(f"{key}. {value}")
        host_choice = input("\nSelect host (1-3): ")
        
        # Show metadata filter options
        print("\nMetadata filter options:")
        metadata_options = {
            '1': 'Location data only',
            '2': 'Collection date only',
            '3': 'Both location and date',
            '4': 'All sequences (no metadata filter)'
        }
        for key, value in metadata_options.items():
            print(f"{key}. {value}")
        metadata_choice = input("\nSelect metadata filter (1-4): ")
        if metadata_choice in ('2', '3'):
            print("\nBEAST format output?")
            print("1. No")
            print("2. Yes")
            beast_input = input("Select (1-2, default 1): ").strip() or '1'
            beast_choice = beast_input if beast_input in ('1', '2') else '1'
        else:
            beast_choice = '1'

        # Alignment type and protein selection (interactive only)
        print("\nDo you want to perform alignment?")
        print("1. Whole genome alignment")
        print("2. Protein (CDS-based, protein-level alignment)")
        print("3. No")
        alignment_choice = input("\nSelect alignment option (1-3): ").strip()
        while alignment_choice not in ('1', '2', '3'):
            print("Please enter 1, 2, or 3")
            alignment_choice = input("Select alignment option (1-3): ").strip()
        args.alignment = alignment_choice
        args.proteins = None
        if alignment_choice == '2':
            print("\nSelect protein(s) to align (comma-separated for multiple, e.g. 1,2 or 1,2,3):")
            print("1. L (RNA-dependent RNA polymerase)")
            print("2. NP (nucleoprotein)")
            print("3. VP35 (polymerase cofactor)")
            print("4. VP40 (matrix protein)")
            print("5. VP30 (minor nucleoprotein)")
            print("6. VP24 (membrane-associated protein)")
            protein_choice = input("\nSelect protein(s) (1-6 or e.g. 1,2): ").strip()
            protein_map = {'1': 'L', '2': 'NP', '3': 'VP35', '4': 'VP40', '5': 'VP30', '6': 'VP24'}
            parts = [p.strip() for p in protein_choice.replace(',', ' ').split() if p.strip()]
            prots = []
            for p in parts:
                if p in protein_map:
                    prots.append(protein_map[p])
                else:
                    print("Invalid protein option; use 1-6 or e.g. 1,2")
                    break
            else:
                if prots:
                    args.proteins = ','.join(prots)
            if not args.proteins:
                args.proteins = 'L'
                print("Defaulting to L.")
        if alignment_choice in ('1', '2'):
            print("\nDo you want to build a phylogeny?")
            print("1. Yes")
            print("2. No")
            phylogeny_choice = input("\nSelect phylogeny option (1-2): ").strip()
            while phylogeny_choice not in ('1', '2'):
                print("Please enter 1 or 2")
                phylogeny_choice = input("Select phylogeny option (1-2): ").strip()
            args.phylogeny = (phylogeny_choice == '1')
        else:
            args.phylogeny = False
    
    # Download sequences
    output_file, query = download_sequences(virus_choices, genome_choice, host_choice, metadata_choice, completeness_threshold)
    
    # Exit if no sequences found
    if output_file is None:
        print("\nExiting due to no sequences found.")
        return
    
    # Read GenBank file
    records = list(SeqIO.parse(output_file, "genbank"))
    if not records:
        print("\nNo valid sequences found after parsing.")
        if os.path.exists(output_file):
            os.remove(output_file)
        return
        
    print(f"\nTotal records found: {len(records)}")
    
    # Get sequences to remove if file provided
    sequences_to_remove = set()
    if args.remove:
        sequences_to_remove = get_sequences_to_remove(args.remove)
    
    # Get summary information
    location_counts, unknown_count, oldest_record, oldest_year, total_sequences = summarize_locations(records, completeness_threshold, virus_choices)
    
    genome_type_str = "complete" if genome_choice == '1' else "all"
    host_str = host_choice if host_choice != '3' else "allhosts"
    base_filename = f"{base_name}_{genome_type_str}_{host_str}"
    fasta_filename = f"filtered_{base_filename}.fasta"
    summary_filename = f"summary_{base_filename}.txt"
    
    # Outgroup only when --phylogeny (for tree rooting)
    outgroup_record = outgroup_species = outgroup_filename = None
    if args.phylogeny:
        outgroup_record, outgroup_species = get_outgroup_reference(virus_choices)
        outgroup_filename = f"outgroup_{outgroup_species.lower()}_{genome_type_str}_{host_str}.fasta"
    
    # Create output directories and files
    os.makedirs("FASTA", exist_ok=True)
    if beast_choice == '2':
        os.makedirs("BEAST_input", exist_ok=True)
    
    if args.phylogeny and outgroup_record is not None:
        with open(outgroup_filename, "w") as outgroup_file:
            if metadata_choice == '1':
                outgroup_record.id = f"{outgroup_record.id}/{outgroup_species}/Uganda"
            elif metadata_choice == '2':
                outgroup_record.id = f"{outgroup_record.id}/04-Aug-2004"
            elif metadata_choice == '3':
                outgroup_record.id = f"{outgroup_record.id}/{outgroup_species}/Uganda/04-Aug-2004"
            else:
                outgroup_record.id = f"{outgroup_record.id}/{outgroup_species}/Uganda"
            outgroup_record.description = ""
            SeqIO.write(outgroup_record, outgroup_file, "fasta")
    
    # Process sequences and write FASTA files
    sequences_written = 0
    sequences_removed = 0
    fasta_path = os.path.join("FASTA", fasta_filename)
    with open(fasta_path, "w") as output:
        for record in records:
            try:
                record_species = (get_record_species_name(record) or virus_choices[0]) if len(virus_choices) > 1 else virus_choices[0]
                if record.id in sequences_to_remove:
                    print(f"Removing sequence: {record.id}")
                    sequences_removed += 1
                    continue
                
                formatted_id = format_record_id(record, record_species, metadata_choice)
                if formatted_id:
                    record.id = formatted_id
                    record.description = ""
                    SeqIO.write(record, output, "fasta")
                    sequences_written += 1
            except Exception as e:
                print(f"Error processing record {record.id}: {str(e)}")
    
    print(f"\nSequences written: {sequences_written}")
    if sequences_removed > 0:
        print(f"Sequences removed: {sequences_removed}")

    run_log = ["Sequences downloaded: %d" % sequences_written]

    # Create location file in FASTA directory
    fasta_location_file = os.path.join("FASTA", "location.txt")
    with open(fasta_location_file, "w") as loc_file:
        for record in records:
            record_species = (get_record_species_name(record) or virus_choices[0]) if len(virus_choices) > 1 else virus_choices[0]
            if record.id != "MF102255.1":
                location = get_location(record)
                if location != "Unknown":
                    formatted_id = format_record_id(record, record_species, metadata_choice)
                    if formatted_id:
                        loc_file.write(f"{formatted_id}\t{location}\n")
        if outgroup_record is not None:
            outgroup_location = "Uganda"
            formatted_outgroup_id = format_record_id(outgroup_record, outgroup_species, metadata_choice)
            loc_file.write(f"{formatted_outgroup_id}\t{outgroup_location}\n")
    
    # Create BEAST files if requested
    if beast_choice == '2':
        beast_fasta = os.path.join("BEAST_input", f"beast_{os.path.basename(fasta_filename)}")
        with open(beast_fasta, "w") as output:
            for record in records:
                record_species = (get_record_species_name(record) or virus_choices[0]) if len(virus_choices) > 1 else virus_choices[0]
                if record.id != "MF102255.1":
                    formatted_id = format_record_id(record, record_species, '2')
                    if formatted_id:
                        date_str = formatted_id.split('/')[-1]
                        decimal_date = convert_to_decimal_date(date_str)
                        if decimal_date:
                            beast_record = record[:]
                            beast_record.id = f"{record.id.split('/')[0]}/{decimal_date}"
                            beast_record.description = ""
                            SeqIO.write(beast_record, output, "fasta")
        if outgroup_record is not None:
            beast_outgroup_path = os.path.join("BEAST_input", f"beast_{os.path.basename(outgroup_filename)}")
            with open(beast_outgroup_path, "w") as output:
                beast_outgroup_rec = outgroup_record[:]
                beast_outgroup_rec.id = f"{outgroup_record.id.split('/')[0]}/2004.591"
                beast_outgroup_rec.description = ""
                SeqIO.write(beast_outgroup_rec, output, "fasta")
        if metadata_choice == '3':
            location_file = os.path.join("BEAST_input", "location.txt")
            with open(location_file, "w") as loc_file:
                for record in records:
                    record_species = (get_record_species_name(record) or virus_choices[0]) if len(virus_choices) > 1 else virus_choices[0]
                    if record.id != "MF102255.1":
                        location = get_location(record)
                        if location != "Unknown":
                            formatted_id = format_record_id(record, record_species, '2')
                            if formatted_id:
                                date_str = formatted_id.split('/')[-1]
                                decimal_date = convert_to_decimal_date(date_str)
                                if decimal_date:
                                    beast_id = f"{record.id.split('/')[0]}/{decimal_date}"
                                    loc_file.write(f"{beast_id}\t{location}\n")
                if outgroup_record is not None:
                    outgroup_date = "2004.591"
                    beast_outgroup_id = f"{outgroup_record.id.split('/')[0]}/{outgroup_date}"
                    loc_file.write(f"{beast_outgroup_id}\tUganda\n")
    
    # Write summary file
    with open(summary_filename, "w") as summary:
        if outgroup_record is not None:
            summary.write("=== Suggested Outgroup (for tree rooting when using --phylogeny) ===\n")
            summary.write(f"Species: {outgroup_species}\n")
            summary.write(f"Sequence ID: {outgroup_record.id}\n")
            summary.write("Type: Reference sequence (RefSeq)\n")
            summary.write(f"File: FASTA/{os.path.basename(outgroup_filename)}\n")
        else:
            summary.write("Outgroup not included (use --phylogeny to add outgroup for tree rooting).\n")
        
        # Add location counts if available
        if location_counts:
            summary.write("\n=== Location Summary ===\n")
            for location, count in sorted(location_counts.items()):
                summary.write(f"{location}: {count} sequences\n")
            
            if unknown_count > 0:
                summary.write(f"Unknown location: {unknown_count} sequences\n")
        
        summary.write(f"\nTotal sequences in final dataset: {total_sequences}\n")
        
        # Only include filtering summary if sequences were actually removed
        removed_count = len(records) - sequences_written
        if removed_count > 0:
            summary.write("\n=== Sequence Filtering Summary ===\n")
            summary.write(f"Initial sequences downloaded: {len(records)}\n")
            summary.write(f"Final sequences in dataset: {sequences_written}\n")
        
        if beast_choice == '2':
            summary.write("\n=== BEAST Input Files ===\n")
            summary.write(f"Main sequences: BEAST_input/beast_{os.path.basename(fasta_filename)}\n")
            if outgroup_record is not None:
                summary.write(f"Outgroup sequence: BEAST_input/beast_{os.path.basename(outgroup_filename)}\n")
            if metadata_choice == '3':
                summary.write(f"Location data: BEAST_input/location.txt\n")
    
    try:
        # Create FASTA directory and move files
        os.makedirs("FASTA", exist_ok=True)
        
        # Move outgroup FASTA to FASTA/ when --phylogeny was used
        if outgroup_filename and os.path.exists(outgroup_filename):
            os.rename(outgroup_filename, os.path.join("FASTA", outgroup_filename))
        
        # Copy per-species consensus FASTA if provided (--c_z, --c_s, etc.)
        consensus_copied = []
        for name, attr in [('zaire', 'consensus_zaire'), ('sudan', 'consensus_sudan'), ('reston', 'consensus_reston'),
                          ('bundibugyo', 'consensus_bundibugyo'), ('tai', 'consensus_tai')]:
            path = getattr(args, attr, None)
            if path:
                dst = copy_consensus_file(path, "FASTA", dest_name="consensus_%s.fasta" % name)
                if dst:
                    consensus_copied.append(name)
                else:
                    print("Failed to copy consensus: %s" % path)
        if consensus_copied:
            run_log.append("Consensus added: %d (%s)" % (len(consensus_copied), ", ".join(consensus_copied)))

        fasta_path = os.path.join("FASTA", fasta_filename)
        outgroup_path = os.path.join("FASTA", outgroup_filename) if outgroup_filename else None
        if os.path.exists(fasta_path) and outgroup_path and os.path.exists(outgroup_path):
            print("Creating location file in FASTA directory...")
            create_fasta_location_file(fasta_path, outgroup_path)
        
        # Run alignment and/or phylogeny based on alignment type
        nthreads = getattr(args, 'threads', 1)
        if nthreads <= 0:
            nthreads = os.cpu_count() or 1
        if args.alignment == '1':
            if args.phylogeny:
                print("\nStarting phylogenetic analysis (whole genome)...")
                create_phylogenetic_tree("FASTA", threads=nthreads)
                print("Phylogenetic analysis completed!")
            else:
                print("\nStarting whole-genome alignment...")
                create_alignment_only("FASTA", threads=nthreads)
                print("Alignment completed!")
        elif args.alignment == '2':
            run_protein_pipeline("FASTA", args.proteins, args.phylogeny, virus_choices, base_name, args, run_log=run_log)
        # Write run log
        if run_log:
            log_path = "ebolaseq_run.log"
            with open(log_path, "w") as f:
                f.write("\n".join(run_log))
            print("Run log written to %s" % log_path)
    except Exception as e:
        print(f"Error during processing: {str(e)}")
        raise  # Add this to see full error traceback
    
    # Clean up
    if os.path.exists("downloaded_genomes.gb"):
        os.remove("downloaded_genomes.gb")
    
    print("\nProcessing complete!")
    print(f"- FASTA files saved in FASTA directory")
    print(f"- Summary information saved to {summary_filename}")
    if beast_choice == '2':
        print(f"- BEAST files saved in BEAST_input directory")
    if args.alignment == '1':
        if args.phylogeny:
            print(f"- Phylogeny output in Phylogeny/")
        else:
            print(f"- Alignment output in Alignment/ (FASTA, MAFFT, Trimmed)")
    elif args.alignment == '2':
        print(f"- Protein alignment output in Alignment/")
        if args.phylogeny:
            print(f"- Phylogeny per protein in Phylogeny/<protein>/")

def download_sequences(virus_choices, genome_choice, host_choice, metadata_choice, completeness_threshold=0):
    """Download sequences from GenBank. virus_choices: list of species (e.g. ['Zaire_ebolavirus'] or multiple for pan/multi)."""
    Entrez.email = "anonymous@example.com"
    if isinstance(virus_choices, str):
        virus_choices = [virus_choices]

    # Reference sequence IDs for each virus type
    reference_sequences = {
        'Zaire_ebolavirus': 'NC_002549.1',
        'Sudan_ebolavirus': 'NC_006432.1',
        'Bundibugyo_ebolavirus': 'NC_014373.1',
        'Tai_Forest_ebolavirus': 'NC_014372.1',
        'Reston_ebolavirus': 'NC_004161.1'
    }

    # Build query with all filters
    query_parts = []

    # Virus species filter: one or more species (OR). NCBI uses organism names with spaces.
    ORGANISM_QUERY = {
        'Zaire_ebolavirus': '("Zaire ebolavirus"[organism] OR "Zaire Ebolavirus"[organism] OR "ZEBOV"[All Fields])',
        'Sudan_ebolavirus': '("Sudan ebolavirus"[organism] OR "Sudan Ebolavirus"[organism])',
        'Bundibugyo_ebolavirus': '("Bundibugyo ebolavirus"[organism] OR "Bundibugyo Ebolavirus"[organism])',
        'Tai_Forest_ebolavirus': '("Tai Forest ebolavirus"[organism] OR "Tai Forest Ebolavirus"[organism])',
        'Reston_ebolavirus': '("Reston ebolavirus"[organism] OR "Reston Ebolavirus"[organism])',
    }
    organism_clauses = [ORGANISM_QUERY.get(v, f'"{v}"[organism]') for v in virus_choices]
    query_parts.append("(" + " OR ".join(organism_clauses) + ")")

    # Add host filter with multiple possible formats
    if host_choice == '1':  # Human
        query_parts.append('("Homo sapiens"[host] OR "human"[host]) NOT ("macaque"[host] OR "monkey"[host] OR "Pan troglodytes"[host] OR "Gorilla"[host] OR "bat"[host] OR "Unknown"[host] OR "cynomolgus macaque"[host] OR "cynomolgus_macaque"[host] OR "cynomolgusmacaque"[host] OR "Macaca fascicularis"[host] OR "Macaca mulatta"[host] OR "rhesus"[host])')
    elif host_choice == '2':  # Non-human
        query_parts.append('("macaque"[host] OR "monkey"[host] OR "Pan troglodytes"[host] OR "Gorilla"[host] OR "bat"[host] OR "cynomolgus macaque"[host] OR "cynomolgus_macaque"[host] OR "cynomolgusmacaque"[host] OR "Macaca fascicularis"[host] OR "Macaca mulatta"[host] OR "rhesus"[host]) NOT ("Homo sapiens"[host] OR "human"[host])')
    # For host_choice == '3' (All hosts), don't add any host filter

    # Add genome completeness filter
    if genome_choice == '1':  # Complete genomes only
        query_parts.append('("complete genome"[Title] OR "complete sequence"[Title])')
    elif genome_choice == '2' and completeness_threshold > 0:  # Partial genomes with threshold
        ref_lengths = [get_reference_length(v) for v in virus_choices]
        min_ref = min(ref_lengths)
        max_ref = max(ref_lengths)
        min_length = int(min_ref * (completeness_threshold / 100))
        query_parts.append(f'"{min_length}"[SLEN]:"{max_ref}"[SLEN]')

    # Combine all query parts
    query = " AND ".join(query_parts)

    output_file = "downloaded_genomes.gb"

    handle = Entrez.esearch(db="nucleotide", term=query, retmax=10000)
    record = Entrez.read(handle)
    handle.close()
    id_list = record["IdList"]

    # Remove reference IDs if they appear in search results (no refs added to download; refs only used as outgroup when --phylogeny)
    ref_ids = [reference_sequences[v] for v in virus_choices]
    for ref_id in ref_ids:
        if ref_id in id_list:
            id_list.remove(ref_id)

    print(f"\nFound {len(id_list)} sequences matching criteria.")

    with open(output_file, "w") as out_handle:
        batch_size = 100
        n_batches = max(1, (len(id_list) - 1) // batch_size + 1)
        for i in range(0, len(id_list), batch_size):
            batch = id_list[i:i + batch_size]
            print(f"Downloading batch {(i//batch_size)+1} of {n_batches}...")
            fetch_handle = Entrez.efetch(db="nucleotide",
                                         id=batch,
                                         rettype="gb",
                                         retmode="text")
            out_handle.write(fetch_handle.read())
            fetch_handle.close()
            time.sleep(1)

    return output_file, query

# Map processing names to reference-length keys (for get_reference_length)
PROCESSING_TO_REF_KEY = {
    'Zaire_ebolavirus': 'Zaire',
    'Sudan_ebolavirus': 'Sudan',
    'Bundibugyo_ebolavirus': 'Bundibugyo',
    'Tai_Forest_ebolavirus': 'Tai Forest',
    'Reston_ebolavirus': 'Reston',
}

def get_reference_length(virus_type):
    """virus_type: either ref key (e.g. 'Zaire') or processing name (e.g. 'Zaire_ebolavirus')."""
    reference_lengths = {
        'Zaire': 18959,
        'Sudan': 18875,
        'Bundibugyo': 18940,
        'Tai Forest': 18935,
        'Reston': 18891
    }
    key = PROCESSING_TO_REF_KEY.get(virus_type, virus_type)
    return reference_lengths.get(key, 18959)

def summarize_locations(records, threshold, virus_choice_or_list):
    """virus_choice_or_list: single species name or list of species (for pan/multi). Counts all records so summary matches FASTA/location.txt."""
    location_counts = {}
    unknown_count = 0
    oldest_record = None
    oldest_year = 9999
    total_sequences = 0

    for record in records:
        # Count every record so summary total and location breakdown match FASTA and location.txt
        location = get_location(record)
        total_sequences += 1
        if location == "Unknown":
            unknown_count += 1
        else:
            location_counts[location] = location_counts.get(location, 0) + 1
        year = get_collection_date(record)
        if year != 9999:
            if year < oldest_year:
                oldest_year = year
                oldest_record = record
            elif year == oldest_year:
                if get_location(record) == 'Democratic_Republic_of_the_Congo':
                    oldest_record = record
                    oldest_year = year
    
    return location_counts, unknown_count, oldest_record, oldest_year, total_sequences

def format_record_id(record, virus_choice, metadata_choice):
    """Format record ID based on metadata choice. Uses accession only as base so repeated calls don't duplicate."""
    base_id = record.id.split('/')[0]
    location = get_location(record)
    date = get_collection_date(record)
    virus_name = virus_choice

    # Manually set location for Zaire ebolavirus reference
    if base_id == "NC_002549.1":
        location = "Democratic_Republic_of_the_Congo"
    
    # Always include virus name and location (if available)
    formatted_id = f"{base_id}/{virus_name}"
    if location != "Unknown":
        formatted_id += f"/{location}"
    else:
        formatted_id += "/Unknown"
        
    # Add date if available (for options 2 and 3)
    if metadata_choice in ['2', '3']:
        for feature in record.features:
            if feature.type == "source":
                if 'collection_date' in feature.qualifiers:
                    date_str = feature.qualifiers['collection_date'][0]
                    if date_str not in ['', 'unknown', 'Unknown', 'NA', 'missing']:
                        try:
                            if date_str.isdigit() and len(date_str) == 4:
                                date_str = f"XX-XX-{date_str}"
                            elif '-' in date_str:
                                parts = date_str.split('-')
                                if len(parts) == 1:
                                    date_str = f"XX-XX-{parts[0]}"
                                elif len(parts) == 2:
                                    date_str = f"XX-{parts[1]}-{parts[0]}"
                            date_str = date_str.replace(" ", "_")
                            formatted_id += f"/{date_str}"
                        except:
                            pass
    
    return formatted_id

def convert_to_decimal_date(date_str):
    """Convert date string to decimal date format"""
    try:
        # Handle XX-XX-YEAR format
        if date_str.startswith('XX-XX-'):
            year = int(date_str.split('-')[-1])
            return float(year)
        
        # Parse the date string
        if '-' in date_str:
            parts = date_str.split('-')
            # Handle different date formats
            if len(parts) == 3:
                # Convert month name to number if needed
                month_map = {
                    'Jan': 1, 'Feb': 2, 'Mar': 3, 'Apr': 4, 'May': 5, 'Jun': 6,
                    'Jul': 7, 'Aug': 8, 'Sep': 9, 'Oct': 10, 'Nov': 11, 'Dec': 12
                }
                
                day = int(parts[0])
                month = parts[1] if parts[1].isdigit() else month_map.get(parts[1], 1)
                year = int(parts[2])
                
                # Create date object
                date = datetime.datetime(year, int(month), day)
                day_of_year = date.timetuple().tm_yday
                
                # Calculate if it's a leap year
                days_in_year = 366 if year % 4 == 0 and (year % 100 != 0 or year % 400 == 0) else 365
                
                # Calculate decimal date
                decimal_date = year + ((day_of_year - 1) / days_in_year)
                return round(decimal_date, 3)
            
        # If we can't parse the date, return just the year if present
        year_match = re.search(r'(19|20)\d{2}', date_str)
        if year_match:
            return float(year_match.group())
            
    except (ValueError, IndexError):
        return None
    
    return None

def create_beast_header(record_id, metadata_choice):
    """Create BEAST format header based on metadata choice"""
    parts = record_id.split('/')
    accession = parts[0]
    
    if metadata_choice == '2':  # Date only
        date_str = parts[1]
        decimal_date = convert_to_decimal_date(date_str)
        if decimal_date:
            return f"{accession}/{decimal_date}"
    elif metadata_choice == '3':  # Both
        location = parts[2]  # Location already has underscores
        date_str = parts[3]
        decimal_date = convert_to_decimal_date(date_str)
        if decimal_date:
            return f"{accession}/{location}/{decimal_date}"
    
    return record_id  # Return original if can't convert

def create_fasta_location_file(fasta_path, outgroup_path):
    """Create location.txt in FASTA directory based on FASTA headers."""
    fasta_dir = os.path.dirname(fasta_path)
    location_file = os.path.join(fasta_dir, "location.txt")
    
    with open(location_file, "w") as loc_file:
        # Write header first
        loc_file.write("taxon\tlocation\n")
        
        # Process main FASTA file
        with open(fasta_path) as f:
            for line in f:
                if line.startswith('>'):
                    header = line[1:].strip()
                    parts = header.split('/')
                    if len(parts) >= 3:
                        location = parts[2]
                        loc_file.write(f"{header}\t{location}\n")
        
        # Process outgroup FASTA file
        with open(outgroup_path) as f:
            for line in f:
                if line.startswith('>'):
                    header = line[1:].strip()
                    parts = header.split('/')
                    if len(parts) >= 3:
                        location = parts[2]
                        loc_file.write(f"{header}\t{location}\n")
                    break  # Only need first header from outgroup

def copy_consensus_file(consensus_file, fasta_dir, dest_name=None):
    """Copy consensus FASTA to FASTA directory. dest_name: optional output filename (e.g. consensus_zaire.fasta)."""
    if not os.path.exists(consensus_file):
        print(f"Warning: Consensus file {consensus_file} not found")
        return None
    filename = dest_name if dest_name else os.path.basename(consensus_file)
    dst = os.path.join(fasta_dir, filename)
    shutil.copy2(consensus_file, dst)
    print(f"Copied consensus: {filename}")
    return dst

def check_dependencies():
    """Check if required external tools are available"""
    dependencies = ['mafft', 'trimal', 'iqtree2']
    missing = []
    
    for tool in dependencies:
        if subprocess.run(['which', tool], capture_output=True).returncode != 0:
            missing.append(tool)
    
    if missing:
        print(f"Error: Missing required tools: {', '.join(missing)}")
        print("Please install these tools and ensure they are in your PATH")
        sys.exit(1)

def check_alignment_dependencies():
    """Check if required external tools are available for alignment only"""
    dependencies = ['mafft', 'trimal']
    missing = []
    
    for tool in dependencies:
        if subprocess.run(['which', tool], capture_output=True).returncode != 0:
            missing.append(tool)
    
    if missing:
        print(f"Error: Missing required tools: {', '.join(missing)}")
        print("Please install these tools and ensure they are in your PATH")
        sys.exit(1)

def check_protein_cds_dependencies():
    """Check minimap2 and pal2nal for protein (CDS) alignment pipeline."""
    for tool in ['minimap2', 'mafft']:
        if subprocess.run(['which', tool], capture_output=True).returncode != 0:
            print("Error: %s not found. Required for protein alignment." % tool)
            sys.exit(1)
    if subprocess.run(['which', 'pal2nal.pl'], capture_output=True).returncode != 0:
        print("Error: pal2nal.pl not found. Install with e.g. conda install -c bioconda pal2nal")
        sys.exit(1)

def _fetch_refseq_cds(accession):
    """Fetch RefSeq FASTA from NCBI."""
    try:
        from urllib.request import urlopen
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=%s&rettype=fasta&retmode=text" % accession
        with urlopen(url, timeout=30) as r:
            return r.read().decode()
    except Exception as e:
        raise RuntimeError("Failed to fetch RefSeq %s: %s" % (accession, e))

def _parse_cigar(cigar_str):
    for m in re.finditer(r"(\d+)([MIDNSHP=X])", cigar_str):
        yield m.group(2), int(m.group(1))

def _extract_region_cds(seq, ref_start_1, ref_end_1, align_ref_start_1, cigar_str):
    out, ref_pos, query_pos = [], align_ref_start_1, 0
    for op, length in _parse_cigar(cigar_str):
        if op in "M=X":
            for _ in range(length):
                if ref_start_1 <= ref_pos <= ref_end_1:
                    out.append(seq[query_pos])
                ref_pos += 1
                query_pos += 1
        elif op == "I":
            if ref_start_1 <= ref_pos - 1 <= ref_end_1:
                for i in range(length):
                    out.append(seq[query_pos + i])
            query_pos += length
        elif op in "D N":
            ref_pos += length
        elif op in "S H P":
            query_pos += length
    return "".join(out)

def _extract_cds_from_sam(sam_path, ref_start, ref_end, query_ids, out_path, min_cds_fraction=0.5):
    ref_len = ref_end - ref_start + 1
    min_len = int(ref_len * min_cds_fraction)
    seqs = OrderedDict()
    alignment_index = 0
    with open(sam_path) as sam:
        for line in sam:
            if line.startswith("@"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 11:
                continue
            qname, flag, rname, pos_s, _, cigar, _, _, _, seq = parts[:10]
            if int(flag) & 0x900 or not seq or seq == "*":
                continue
            pos = int(pos_s)
            if pos > ref_end:
                continue
            extracted = _extract_region_cds(seq, ref_start, ref_end, pos, cigar)
            if len(extracted) < min_len:
                continue
            name = query_ids[alignment_index] if query_ids and alignment_index < len(query_ids) else qname
            alignment_index += 1
            seqs[name] = extracted
    with open(out_path, "w") as f:
        for name, s in seqs.items():
            f.write(">%s\n" % name)
            for i in range(0, len(s), 60):
                f.write(s[i:i+60] + "\n")
    n_input = len(query_ids) if query_ids else 0
    kept = set(seqs.keys())
    dropped_ids = [q for q in (query_ids or []) if q not in kept]
    return (len(seqs), n_input, dropped_ids)

def _translate_cds_to_protein(rec):
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    s = rec.seq
    if len(s) % 3 != 0:
        s = s[: (len(s) // 3) * 3]
    aa = str(s.translate(table=1, to_stop=False, cds=False))
    return SeqRecord(Seq(aa), id=rec.id, description="")

def _translate_cds_fasta(cds_path, protein_path):
    records = list(SeqIO.parse(cds_path, "fasta"))
    if not records:
        open(protein_path, "w").close()
        return
    proteins = [_translate_cds_to_protein(rec) for rec in records]
    with open(protein_path, "w") as f:
        SeqIO.write(proteins, f, "fasta")

def _backtranslate_pal2nal(protein_aln_path, cds_path, out_path):
    """Back-translate protein alignment to codon alignment using PAL2NAL."""
    from Bio import AlignIO, SeqIO
    from Bio.SeqRecord import SeqRecord
    cds_by_id = {}
    for rec in SeqIO.parse(cds_path, "fasta"):
        cds_by_id[rec.id] = rec.seq
        cds_by_id[rec.id.replace("/", "_")] = rec.seq
        cds_by_id[rec.id.replace("_", "/")] = rec.seq
    aln = AlignIO.read(protein_aln_path, "fasta")
    cds_ordered = []
    for rec in aln:
        seq = cds_by_id.get(rec.id) or cds_by_id.get(rec.id.replace("/", "_")) or cds_by_id.get(rec.id.replace("_", "/"))
        if seq is None:
            return False
        cds_ordered.append(SeqRecord(seq, id=rec.id, description=""))
    fd, tmp_cds = tempfile.mkstemp(suffix=".fasta")
    os.close(fd)
    try:
        with open(tmp_cds, "w") as f:
            SeqIO.write(cds_ordered, f, "fasta")
        ret = subprocess.run(
            ["pal2nal.pl", protein_aln_path, tmp_cds, "-output", "fasta"],
            capture_output=True, text=True, timeout=600,
        )
        if ret.returncode != 0:
            return False
        with open(out_path, "w") as f:
            f.write(ret.stdout)
        return True
    finally:
        try:
            os.unlink(tmp_cds)
        except OSError:
            pass

def run_protein_cds_pipeline(outroot_abs, pairs, proteins, base_name=None, min_cds_fraction=0.5, threads=1):
    """
    Run CDS extraction, protein alignment, back-translate (inlined from ebolaS).
    Writes only to outroot_abs/MAFFT/<protein>/ and outroot_abs/Trimmed/<protein>/.
    outroot_abs: absolute path to output root (e.g. Alignment dir).
    pairs: list of (species_name, absolute_path_to_fasta).
    base_name: unused (kept for API compatibility).
    """
    work_dir = os.path.join(outroot_abs, "_work")
    os.makedirs(work_dir, exist_ok=True)
    # extraction_counts[(species, protein)] = (n_passed, n_input)
    extraction_counts = {}
    for species, fpath in pairs:
        species_dir = os.path.join(work_dir, species)
        os.makedirs(species_dir, exist_ok=True)
        ref_fa = os.path.join(species_dir, "ref.fa")
        try:
            ref_fasta = _fetch_refseq_cds(REFSEQ_IDS_CDS[species])
            with open(ref_fa, "w") as f:
                f.write(ref_fasta)
        except RuntimeError as e:
            with open(fpath) as fin, open(ref_fa, "w") as fout:
                fout.write(fin.readline())
                for line in fin:
                    if line.startswith(">"):
                        break
                    fout.write(line)
        query_ids = [rec.id for rec in SeqIO.parse(fpath, "fasta")]
        sam_path = os.path.join(species_dir, "_aln.sam")
        with open(sam_path, "w") as f:
            subprocess.run(
                ["minimap2", "-x", "asm5", "-a", "-t", str(threads), ref_fa, fpath],
                stdout=f, stderr=subprocess.DEVNULL,
            )
        protein_regions = SPECIES_REGIONS_CDS[species]
        for protein in proteins:
            if protein not in protein_regions:
                continue
            ref_start, ref_end = protein_regions[protein]
            protdir = os.path.join(species_dir, protein)
            os.makedirs(protdir, exist_ok=True)
            n_passed, n_input, dropped_ids = _extract_cds_from_sam(sam_path, ref_start, ref_end, query_ids, os.path.join(protdir, "cds.fasta"), min_cds_fraction)
            extraction_counts[(species, protein)] = (n_passed, n_input, dropped_ids)
        try:
            os.unlink(sam_path)
        except OSError:
            pass
    subdirs = sorted([os.path.join(work_dir, s) for s, _ in pairs])
    protein_dropped_ids = {}
    print("  (Per protein: only CDS coverage step drops sequences; each protein can keep a different number.)")
    for protein in proteins:
        protdir = os.path.join(work_dir, protein)
        os.makedirs(protdir, exist_ok=True)
        combined_cds = os.path.join(protdir, "cds.fasta")
        with open(combined_cds, "w") as fout:
            for subdir in subdirs:
                cds_file = os.path.join(subdir, protein, "cds.fasta")
                if not os.path.isfile(cds_file):
                    continue
                prefix = os.path.basename(subdir)
                for rec in SeqIO.parse(cds_file, "fasta"):
                    rec.id = "%s|%s" % (prefix, rec.id)
                    rec.description = ""
                    SeqIO.write(rec, fout, "fasta")
        if os.path.getsize(combined_cds) == 0:
            continue
        _translate_cds_fasta(combined_cds, os.path.join(protdir, "protein.fasta"))
        protein_aln = os.path.join(protdir, "protein_aln.fasta")
        with open(protein_aln, "w") as f:
            r = subprocess.run(
                ["mafft", "--auto", "--thread", str(min(threads, 32)), os.path.join(protdir, "protein.fasta")],
                stdout=f, stderr=subprocess.PIPE, text=True,
            )
            if r.returncode != 0:
                print("MAFFT failed for %s: %s" % (protein, r.stderr))
                continue
        cds_aligned = os.path.join(protdir, "cds_aligned.fasta")
        if not _backtranslate_pal2nal(protein_aln, combined_cds, cds_aligned):
            print("PAL2NAL failed for %s" % protein)
            continue
        n_in_alignment = len(list(SeqIO.parse(cds_aligned, "fasta")))
        total_input = sum(extraction_counts.get((s, protein), (0, 0, []))[1] for s, _ in pairs)
        total_passed = sum(extraction_counts.get((s, protein), (0, 0, []))[0] for s, _ in pairs)
        dropped = total_input - total_passed
        dropped_ids_this_protein = []
        for s, _ in pairs:
            dropped_ids_this_protein.extend(extraction_counts.get((s, protein), (0, 0, []))[2])
        protein_dropped_ids[protein] = dropped_ids_this_protein
        per_species = " ".join("%s %d→%d" % (s, extraction_counts.get((s, protein), (0, 0, []))[1], extraction_counts.get((s, protein), (0, 0, []))[0]) for s, _ in pairs)
        print("  %s: %d assigned to species → %d in alignment (%d dropped: CDS coverage < %.0f%% of reference or no alignment). Per species: %s" % (
            protein, total_input, n_in_alignment, dropped, min_cds_fraction * 100, per_species))
        print("  Wrote %s" % cds_aligned)
        # Mirror into Alignment/MAFFT/<protein>/ (originals + aligned) and Alignment/Trimmed/<protein>/
        mafft_prot_dir = os.path.join(outroot_abs, "MAFFT", protein)
        trim_prot_dir = os.path.join(outroot_abs, "Trimmed", protein)
        os.makedirs(mafft_prot_dir, exist_ok=True)
        os.makedirs(trim_prot_dir, exist_ok=True)
        protein_fasta = os.path.join(protdir, "protein.fasta")
        for name, src in [
            ("cds.fasta", combined_cds),
            ("protein.fasta", protein_fasta),
            ("protein_aln.fasta", protein_aln),
            ("cds_aligned.fasta", cds_aligned),
        ]:
            shutil.copy2(src, os.path.join(mafft_prot_dir, name))
        trimmed_fasta = os.path.join(trim_prot_dir, "cds_trimmed.fasta")
        r = subprocess.run(
            ["trimal", "-in", cds_aligned, "-out", trimmed_fasta, "-automated1"],
            capture_output=True, text=True,
        )
        if r.returncode != 0:
            print("  TrimAl failed for %s, copying untrimmed" % protein)
            shutil.copy2(cds_aligned, trimmed_fasta)
        else:
            print("  Wrote %s" % trimmed_fasta)
    try:
        shutil.rmtree(work_dir, ignore_errors=True)
    except Exception:
        pass
    return protein_dropped_ids

def split_fasta_by_species(fasta_path, alignment_dir, consensus_paths=None):
    """
    Split a combined FASTA (headers accession/species/...) into per-species files.
    consensus_paths: dict e.g. {'zaire': '/path/consensus_z.fasta'} to prepend per species.
    Returns (list of (ebolaS_species_name, path), n_assigned, n_skipped, skipped_ids).
    Records skipped: ID has no '/' or species (parts[1]) not in PROCESSING_TO_EBOLAS_SPECIES.
    """
    from Bio import SeqIO
    species_to_records = {}
    skipped_ids = []
    n_total = 0
    for rec in SeqIO.parse(fasta_path, "fasta"):
        n_total += 1
        parts = rec.id.split("/")
        if len(parts) < 2:
            skipped_ids.append(rec.id)
            continue
        proc_name = parts[1]
        ebolaS_name = PROCESSING_TO_EBOLAS_SPECIES.get(proc_name) or SPECIES_ID_ALIASES.get(proc_name)
        if ebolaS_name is None:
            # Fallback: look up species from NCBI using accession (ORGANISM field) so any header format works
            ebolaS_name = get_ebolas_species_for_accession(parts[0])
        if ebolaS_name is None:
            skipped_ids.append(rec.id)
            continue
        if ebolaS_name not in species_to_records:
            species_to_records[ebolaS_name] = []
        species_to_records[ebolaS_name].append(rec)
    n_assigned = sum(len(r) for r in species_to_records.values())
    n_skipped = n_total - n_assigned
    os.makedirs(alignment_dir, exist_ok=True)
    out_pairs = []
    for ebolaS_name, records in species_to_records.items():
        out_path = os.path.join(alignment_dir, "input_%s.fasta" % ebolaS_name)
        with open(out_path, "w") as f:
            if consensus_paths and ebolaS_name in consensus_paths:
                cons_path = consensus_paths.get(ebolaS_name)
                if cons_path and os.path.exists(cons_path):
                    with open(cons_path) as cf:
                        f.write(cf.read())
            SeqIO.write(records, f, "fasta")
        out_pairs.append((ebolaS_name, os.path.abspath(out_path)))
    return out_pairs, n_assigned, n_skipped, skipped_ids


def run_protein_pipeline(source_fasta_dir, proteins_str, do_phylogeny, virus_choices, base_name, args, run_log=None):
    """Run protein (CDS) alignment pipeline (built-in), then optionally IQTree per protein. run_log: optional list to append run summary for ebolaseq_run.log."""
    alignment_dir = "Alignment"
    os.makedirs(alignment_dir, exist_ok=True)
    # Build combined = single FASTA for alignment (all source FASTA files merged, duplicate IDs removed)
    combined = os.path.join(alignment_dir, "FASTA", "Ebola_Combined.fasta")
    os.makedirs(os.path.join(alignment_dir, "FASTA"), exist_ok=True)
    # Map consensus filename to processing name so split by species can assign them
    CONSENSUS_FILE_TO_SPECIES = {
        "consensus_zaire": "Zaire_ebolavirus", "consensus_sudan": "Sudan_ebolavirus",
        "consensus_reston": "Reston_ebolavirus", "consensus_bundibugyo": "Bundibugyo_ebolavirus",
        "consensus_tai_forest": "Tai_Forest_ebolavirus", "consensus_tai": "Tai_Forest_ebolavirus",
    }
    seen_ids = set()
    combined_records = []
    n_read = 0
    n_from_download = 0
    n_from_consensus = 0
    duplicate_ids = []
    for fname in sorted(os.listdir(source_fasta_dir)):
        if not fname.endswith(('.fasta', '.fa')) or fname == "Ebola_Combined.fasta":
            continue
        fpath = os.path.join(source_fasta_dir, fname)
        from_consensus = fname.startswith("consensus_")
        base = fname.replace(".fasta", "").replace(".fa", "")
        consensus_species = CONSENSUS_FILE_TO_SPECIES.get(base) if from_consensus else None
        consensus_seq_index = 0
        for rec in SeqIO.parse(fpath, "fasta"):
            n_read += 1
            if from_consensus:
                n_from_consensus += 1
                # Normalize consensus IDs to accession/species/... so split by species assigns them
                if consensus_species:
                    parts = rec.id.split("/")
                    if len(parts) < 2 or PROCESSING_TO_EBOLAS_SPECIES.get(parts[1]) is None:
                        base_id = (rec.id.split()[0] if rec.id else "consensus").replace("/", "_")
                        rec.id = "%s_%d/%s/consensus" % (base_id, consensus_seq_index, consensus_species)
                        consensus_seq_index += 1
            else:
                n_from_download += 1
            if rec.id in seen_ids:
                duplicate_ids.append(rec.id)
                print("Skipping duplicate ID in combined: %s" % rec.id)
                continue
            seen_ids.add(rec.id)
            combined_records.append(rec)
    with open(combined, "w") as outfile:
        SeqIO.write(combined_records, outfile, "fasta")
    n_dup = n_read - len(combined_records)
    print("Combined (for alignment): %d from download FASTA + %d from consensus = %d total; %d duplicate IDs removed → %d sequences in combined." % (
        n_from_download, n_from_consensus, n_read, n_dup, len(combined_records)))
    if run_log is not None and duplicate_ids:
        run_log.append("Duplicates dropped: %d" % len(duplicate_ids))
        for did in duplicate_ids:
            run_log.append("  duplicate: %s" % did)
    # Do not pass consensus_paths: combined FASTA already includes consensus from source_fasta_dir,
    # so prepending again would double-count (e.g. 93 in combined + 22 prepended = 115 "assigned").
    with tempfile.TemporaryDirectory(prefix="ebolaseq_by_species_") as by_species_dir:
        pairs, n_assigned, n_skipped, split_skipped_ids = split_fasta_by_species(combined, by_species_dir, consensus_paths=None)
        if not pairs:
            print("Error: no sequences could be assigned to a species for protein alignment.")
            return
        n_combined = len(combined_records)
        if n_skipped > 0:
            print("Split by species: %d in combined → %d assigned to species (%d skipped: ID not accession/species/... or species not recognized)." % (n_combined, n_assigned, n_skipped))
            for sid in split_skipped_ids:
                print("  skipped: %s" % sid)
            if run_log is not None:
                run_log.append("Split by species: %d in combined → %d assigned (%d skipped, ID or species not recognized)" % (n_combined, n_assigned, n_skipped))
                for sid in split_skipped_ids:
                    run_log.append("  split_skipped: %s" % sid)
        elif n_assigned != n_combined:
            print("Split by species: %d in combined → %d assigned to species." % (n_combined, n_assigned))
        _num_to_protein = {'1': 'L', '2': 'NP', '3': 'VP35', '4': 'VP40', '5': 'VP30', '6': 'VP24'}
        raw = [x.strip().upper() for x in proteins_str.split(",") if x.strip()] if proteins_str else ["L"]
        proteins = [_num_to_protein.get(p, p) for p in raw]
        check_protein_cds_dependencies()
        pairs_abs = [(s, os.path.abspath(p)) for s, p in pairs]
        print("\nRunning protein (CDS) alignment pipeline...")
        min_cds = getattr(args, 'min_cds_fraction', 0.5)
        nthreads = getattr(args, 'threads', 1)
        if nthreads <= 0:
            nthreads = os.cpu_count() or 1
        protein_dropped = run_protein_cds_pipeline(os.path.abspath(alignment_dir), pairs_abs, proteins, base_name=base_name, min_cds_fraction=min_cds, threads=nthreads)
        if run_log is not None and protein_dropped:
            run_log.append("Dropped (protein coverage not met):")
            for protein, ids in protein_dropped.items():
                if not ids:
                    continue
                run_log.append("  %s: %d dropped" % (protein, len(ids)))
                for sid in ids:
                    run_log.append("    %s" % sid)
    if do_phylogeny:
        check_dependencies()
        phylogeny_base = "Phylogeny"
        os.makedirs(phylogeny_base, exist_ok=True)
        trimmed_dir = os.path.join(alignment_dir, "Trimmed")
        mafft_dir = os.path.join(alignment_dir, "MAFFT")
        for protein in proteins:
            trimmed_fasta = os.path.join(trimmed_dir, protein, "cds_trimmed.fasta")
            cds_aligned = os.path.join(mafft_dir, protein, "cds_aligned.fasta")
            src = trimmed_fasta if os.path.isfile(trimmed_fasta) else cds_aligned
            if not os.path.isfile(src):
                continue
            protdir = os.path.join(phylogeny_base, protein)
            os.makedirs(protdir, exist_ok=True)
            dest = os.path.join(protdir, "cds_aligned.fasta")
            shutil.copy2(src, dest)
            cwd = os.getcwd()
            os.chdir(protdir)
            iqtree_nt = "AUTO" if (threads <= 0) else str(threads)
            r = subprocess.run(["iqtree2", "-s", "cds_aligned.fasta", "-m", "MFP", "-bb", "10000", "-nt", iqtree_nt],
                               shell=False, stderr=subprocess.PIPE, text=True)
            os.chdir(cwd)
            if r.returncode != 0:
                print("IQTree2 failed for %s: %s" % (protein, r.stderr))
            else:
                print("IQTree2 completed for %s in Phylogeny/%s/" % (protein, protein))
    print("Protein pipeline completed. Output: Alignment/FASTA, Alignment/MAFFT, Alignment/Trimmed")
    if do_phylogeny:
        print("Note: If IQ-TREE reported 'X identical sequences ignored', those sequences are still in the alignment; IQ-TREE uses only unique sequences for the tree and places identical ones at the same branch (e.g. 84 in alignment → 62 unique + 22 identical).")


def run_alignment_pipeline(source_fasta_dir, threads=1):
    """
    Build Alignment folder: Alignment/FASTA (combined), Alignment/MAFFT, Alignment/Trimmed.
    source_fasta_dir: folder with original FASTA files (e.g. 'FASTA').
    threads: for MAFFT (0 = use all CPUs).
    Returns (trimmed_fasta_path, True) on success, (None, False) on failure.
    """
    check_alignment_dependencies()
    alignment_dir = "Alignment"
    align_fasta_dir = os.path.join(alignment_dir, "FASTA")
    mafft_dir = os.path.join(alignment_dir, "MAFFT")
    trimmed_dir = os.path.join(alignment_dir, "Trimmed")
    os.makedirs(align_fasta_dir, exist_ok=True)
    os.makedirs(mafft_dir, exist_ok=True)
    os.makedirs(trimmed_dir, exist_ok=True)

    combined_fasta = os.path.join(align_fasta_dir, "Ebola_Combined.fasta")
    print("Creating combined FASTA file in Alignment/FASTA...")
    seen_ids = set()
    combined_records = []
    n_read = 0
    for fasta in os.listdir(source_fasta_dir):
        if not fasta.endswith(('.fasta', '.fa')) or fasta == "Ebola_Combined.fasta":
            continue
        fasta_path = os.path.join(source_fasta_dir, fasta)
        for rec in SeqIO.parse(fasta_path, "fasta"):
            n_read += 1
            if rec.id in seen_ids:
                print("Skipping duplicate ID in combined: %s" % rec.id)
                continue
            seen_ids.add(rec.id)
            combined_records.append(rec)
    with open(combined_fasta, "w") as outfile:
        SeqIO.write(combined_records, outfile, "fasta")
    if n_read > len(combined_records):
        print("Combined FASTA: %d unique sequences (skipped %d duplicate IDs). Source FASTA unchanged." % (len(combined_records), n_read - len(combined_records)))
    print("Combined FASTA file created")

    aligned_fasta = os.path.join(mafft_dir, "Ebola_Combined_aligned.fasta")
    print("Running MAFFT alignment...")
    mafft_thread = "-1" if (threads <= 0) else str(min(threads, 32))
    mafft_cmd = f"mafft --thread {mafft_thread} --auto {combined_fasta} > {aligned_fasta}"
    result = subprocess.run(mafft_cmd, shell=True, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Error during MAFFT alignment: {result.stderr}")
        return None, False
    print("MAFFT alignment complete")

    trimmed_fasta = os.path.join(trimmed_dir, "Ebola_trimmed.fasta")
    print("Running TrimAl...")
    trimal_cmd = f"trimal -in {aligned_fasta} -out {trimmed_fasta} -automated1"
    result = subprocess.run(trimal_cmd, shell=True, stderr=subprocess.PIPE, text=True)
    if result.returncode != 0:
        print(f"Error during TrimAl: {result.stderr}")
        return None, False
    print("TrimAl trimming complete")
    return trimmed_fasta, True


def create_alignment_only(source_fasta_dir, threads=1):
    """Create alignment in Alignment/ (FASTA, MAFFT, Trimmed) from original FASTA folder. threads: for MAFFT (0 = all CPUs)."""
    print("\nStarting alignment pipeline...")
    try:
        trimmed_fasta, ok = run_alignment_pipeline(source_fasta_dir, threads=threads)
        if not ok:
            print("Alignment pipeline completed with errors.")
            return
        print("\nAlignment pipeline completed successfully!")
        print("Output: Alignment/FASTA (combined), Alignment/MAFFT, Alignment/Trimmed")
    except Exception as e:
        print(f"Error during alignment: {str(e)}")


def create_phylogenetic_tree(source_fasta_dir, threads=1):
    """Build Alignment/ (if needed), then run IQTree2 in Phylogeny/ using trimmed alignment. threads: for MAFFT and IQTree2 (0 = all CPUs)."""
    print("\nStarting phylogenetic analysis pipeline...")
    try:
        check_dependencies()
        trimmed_fasta, ok = run_alignment_pipeline(source_fasta_dir, threads=threads)
        if not ok:
            print("Alignment step failed; phylogeny not run.")
            return
        phylogeny_dir = "Phylogeny"
        os.makedirs(phylogeny_dir, exist_ok=True)
        trimmed_copy = os.path.join(phylogeny_dir, "Ebola_trimmed.fasta")
        shutil.copy2(trimmed_fasta, trimmed_copy)
        print("Running IQTree2 in Phylogeny/...")
        current_dir = os.getcwd()
        os.chdir(phylogeny_dir)
        iqtree_nt = "AUTO" if (threads <= 0) else str(threads)
        iqtree_cmd = f"iqtree2 -s Ebola_trimmed.fasta -m MFP -bb 10000 -nt {iqtree_nt}"
        result = subprocess.run(iqtree_cmd, shell=True, stderr=subprocess.PIPE, text=True)
        os.chdir(current_dir)
        if result.returncode != 0:
            print(f"Error during IQTree2: {result.stderr}")
            return
        print("IQTree2 analysis complete")
        print("\nPhylogenetic analysis pipeline completed successfully!")
        print("Output: Alignment/ (FASTA, MAFFT, Trimmed) and Phylogeny/ (IQTree2 results)")
    except Exception as e:
        print(f"Error during phylogenetic analysis: {str(e)}")

def write_filtered_sequences(records, virus_choice_or_list, metadata_choice, output_file):
    """Write filtered sequences to FASTA file and check for duplicates. virus_choice_or_list: single species or list for pan/multi."""
    virus_choices = virus_choice_or_list if isinstance(virus_choice_or_list, list) else [virus_choice_or_list]
    seen_headers = set()
    seen_accessions = set()
    unique_records = []
    duplicates_found = False

    for record in records:
        record_species = (get_record_species_name(record) or virus_choices[0]) if len(virus_choices) > 1 else virus_choices[0]
        base_accession = record.id.split('.')[0]
        formatted_id = format_record_id(record, record_species, metadata_choice)
        
        # Skip if we've seen either the full formatted ID or the base accession
        if formatted_id in seen_headers or base_accession in seen_accessions:
            print(f"Found duplicate sequence: {formatted_id} (accession: {base_accession})")
            duplicates_found = True
            continue
            
        seen_headers.add(formatted_id)
        seen_accessions.add(base_accession)
        record.id = formatted_id
        record.description = formatted_id
        unique_records.append(record)
    
    if duplicates_found:
        print(f"Removed duplicate sequences. Keeping {len(unique_records)} unique sequences.")
    
    # Write unique sequences to file
    SeqIO.write(unique_records, output_file, "fasta")
    return len(unique_records)

def remove_duplicate_sequences(fasta_file):
    """Remove sequences with duplicate headers from a FASTA file"""
    seen_headers = set()
    unique_records = []
    duplicates_found = False
    
    # Read all records from the file
    records = list(SeqIO.parse(fasta_file, "fasta"))
    
    for record in records:
        if record.id in seen_headers:
            print(f"Found duplicate header: {record.id}")
            duplicates_found = True
            continue
        seen_headers.add(record.id)
        unique_records.append(record)
    
    if duplicates_found:
        print(f"Writing {len(unique_records)} unique sequences back to file...")
        SeqIO.write(unique_records, fasta_file, "fasta")
        print("Duplicate sequences removed")
    
    return duplicates_found

if __name__ == "__main__":
    cli_main() 
