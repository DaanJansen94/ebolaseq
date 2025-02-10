from Bio import Entrez
from Bio import SeqIO
import sys
import subprocess
import os
import re
import datetime
import time

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

def get_outgroup_reference(virus_choice):
    Entrez.email = "anonymous@example.com"
    
    # Dictionary of reference sequences for each species (RefSeq accessions)
    refseq_dict = {
        'Zaire_ebolavirus': {'outgroup': 'Sudan_ebolavirus', 'refseq': 'NC_006432'},
        'Sudan_ebolavirus': {'outgroup': 'Zaire_ebolavirus', 'refseq': 'NC_002549'},
        'Bundibugyo_ebolavirus': {'outgroup': 'Zaire_ebolavirus', 'refseq': 'NC_002549'},
        'Tai_Forest_ebolavirus': {'outgroup': 'Zaire_ebolavirus', 'refseq': 'NC_002549'},
        'Reston_ebolavirus': {'outgroup': 'Zaire_ebolavirus', 'refseq': 'NC_002549'}
    }
    
    outgroup_info = refseq_dict[virus_choice]
    refseq_id = outgroup_info['refseq']
    
    # Download the reference sequence
    handle = Entrez.efetch(db="nucleotide", id=refseq_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    
    return record, outgroup_info['outgroup']

# ... (copy ALL your original functions exactly as they were)

def main():
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
        '2': 'Chimpanzee (Pan troglodytes)',
        '3': 'Gorilla (Gorilla sp.)',
        '4': 'All hosts'
    }
    
    # Show virus menu
    print("\nAvailable Ebola virus species:")
    for key, value in virus_options.items():
        print(f"{key}. {value}")
    
    # Get all choices from user
    choice = input("\nSelect the virus type (1-5): ")
    virus_display_name = virus_options[choice]
    virus_choice = virus_processing_names[choice]
    virus_type = virus_choice.split('_')[0]
    
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
    host_choice = input("\nSelect host (1-4): ")
    
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
    
    # Get BEAST choice if applicable
    beast_choice = '1'
    if metadata_choice in ['2', '3']:
        print("\nDo you want to generate BEAST input format?")
        print("1. No")
        print("2. Yes")
        beast_choice = input("\nSelect BEAST option (1-2): ")
    
    # Download sequences
    output_file, query = download_sequences(virus_choice, genome_choice, host_choice, metadata_choice, completeness_threshold)
    
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
    
    # Get summary information
    location_counts, unknown_count, oldest_record, oldest_year, total_sequences = summarize_locations(records, completeness_threshold, virus_choice)
    
    # Get outgroup reference sequence
    outgroup_record, outgroup_species = get_outgroup_reference(virus_choice)
    
    # Create output directories and files
    os.makedirs("FASTA", exist_ok=True)
    if beast_choice == '2':
        os.makedirs("BEAST_input", exist_ok=True)
    
    # Create output filenames
    genome_type_str = "complete" if genome_choice == '1' else "all"
    host_str = host_choice if host_choice != '4' else "allhosts"
    base_filename = f"{virus_choice.lower()}_{genome_type_str}_{host_str}"
    fasta_filename = f"filtered_{base_filename}.fasta"
    outgroup_filename = f"outgroup_{outgroup_species.lower()}_{genome_type_str}_{host_str}.fasta"
    summary_filename = f"summary_{base_filename}.txt"
    
    # Write outgroup to separate FASTA file
    with open(outgroup_filename, "w") as outgroup_file:
        # Format outgroup header based on metadata choice
        if metadata_choice == '1':  # Location only
            outgroup_record.id = f"{outgroup_record.id}/{outgroup_species}/Uganda"
        elif metadata_choice == '2':  # Date only
            outgroup_record.id = f"{outgroup_record.id}/04-Aug-2004"
        elif metadata_choice == '3':  # Both
            outgroup_record.id = f"{outgroup_record.id}/{outgroup_species}/Uganda/04-Aug-2004"
        else:  # All sequences
            outgroup_record.id = f"{outgroup_record.id}/{outgroup_species}/Uganda"
        
        outgroup_record.description = ""
        SeqIO.write(outgroup_record, outgroup_file, "fasta")
    
    # Process sequences and write FASTA files
    sequences_written = 0
    fasta_path = os.path.join("FASTA", fasta_filename)
    with open(fasta_path, "w") as output:
        for record in records:
            try:
                if record.id != "MF102255.1":  # Skip known problematic sequence
                    formatted_id = format_record_id(record, virus_choice, metadata_choice)
                    if formatted_id:
                        record.id = formatted_id
                        record.description = ""
                        SeqIO.write(record, output, "fasta")
                        sequences_written += 1
            except Exception as e:
                print(f"Error processing record {record.id}: {str(e)}")
    
    # Create location file in FASTA directory
    fasta_location_file = os.path.join("FASTA", "location.txt")
    with open(fasta_location_file, "w") as loc_file:
        for record in records:
            if record.id != "MF102255.1":
                location = get_location(record)
                if location != "Unknown":
                    formatted_id = format_record_id(record, virus_choice, metadata_choice)
                    if formatted_id:
                        loc_file.write(f"{formatted_id}\t{location}\n")
        
        # Add outgroup to location file
        outgroup_location = "Uganda"
        formatted_outgroup_id = format_record_id(outgroup_record, outgroup_species, metadata_choice)
        loc_file.write(f"{formatted_outgroup_id}\t{outgroup_location}\n")
    
    # Create BEAST files if requested
    if beast_choice == '2':
        beast_fasta = os.path.join("BEAST_input", f"beast_{os.path.basename(fasta_filename)}")
        beast_outgroup = os.path.join("BEAST_input", f"beast_{os.path.basename(outgroup_filename)}")
        
        # Create BEAST format files
        with open(beast_fasta, "w") as output:
            for record in records:
                if record.id != "MF102255.1":
                    formatted_id = format_record_id(record, virus_choice, '2')
                    if formatted_id:
                        date_str = formatted_id.split('/')[-1]
                        decimal_date = convert_to_decimal_date(date_str)
                        if decimal_date:
                            beast_record = record[:]
                            beast_record.id = f"{record.id.split('/')[0]}/{decimal_date}"
                            beast_record.description = ""
                            SeqIO.write(beast_record, output, "fasta")
        
        # Create BEAST format outgroup
        with open(beast_outgroup, "w") as output:
            beast_outgroup = outgroup_record[:]
            beast_outgroup.id = f"{outgroup_record.id.split('/')[0]}/2004.591"
            beast_outgroup.description = ""
            SeqIO.write(beast_outgroup, output, "fasta")
        
        # Create location file if needed
        if metadata_choice == '3':
            location_file = os.path.join("BEAST_input", "location.txt")
            with open(location_file, "w") as loc_file:
                for record in records:
                    if record.id != "MF102255.1":
                        location = get_location(record)
                        if location != "Unknown":
                            formatted_id = format_record_id(record, virus_choice, '2')
                            if formatted_id:
                                date_str = formatted_id.split('/')[-1]
                                decimal_date = convert_to_decimal_date(date_str)
                                if decimal_date:
                                    beast_id = f"{record.id.split('/')[0]}/{decimal_date}"
                                    loc_file.write(f"{beast_id}\t{location}\n")
                
                # Add outgroup
                outgroup_date = "2004.591"
                beast_outgroup_id = f"{outgroup_record.id.split('/')[0]}/{outgroup_date}"
                loc_file.write(f"{beast_outgroup_id}\tUganda\n")
    
    # Write summary file
    with open(summary_filename, "w") as summary:
        summary.write("=== Suggested Outgroup (Reference Sequence from Sister Species) ===\n")
        summary.write(f"Species: {outgroup_species}\n")
        summary.write(f"Sequence ID: {outgroup_record.id}\n")
        summary.write("Type: Reference sequence (RefSeq)\n")
        summary.write(f"File: FASTA/{os.path.basename(outgroup_filename)}\n")
        
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
        
        # Update BEAST section if it was created
        if beast_choice == '2':
            summary.write("\n=== BEAST Input Files ===\n")
            summary.write(f"Main sequences: BEAST_input/beast_{os.path.basename(fasta_filename)}\n")
            summary.write(f"Outgroup sequence: BEAST_input/beast_{os.path.basename(outgroup_filename)}\n")
            if metadata_choice == '3':
                summary.write(f"Location data: BEAST_input/location.txt\n")
    
    try:
        # Create FASTA directory if it doesn't exist
        os.makedirs("FASTA", exist_ok=True)
        
        # Move both FASTA files to FASTA directory
        if os.path.exists(fasta_filename):
            os.rename(fasta_filename, os.path.join("FASTA", fasta_filename))
        
        if os.path.exists(outgroup_filename):
            os.rename(outgroup_filename, os.path.join("FASTA", outgroup_filename))
        
        # Create location.txt in FASTA directory only if both files exist
        fasta_path = os.path.join("FASTA", fasta_filename)
        outgroup_path = os.path.join("FASTA", outgroup_filename)
        
        if os.path.exists(fasta_path) and os.path.exists(outgroup_path):
            print("Creating location file in FASTA directory...")
            create_fasta_location_file(fasta_path, outgroup_path)
        else:
            print("Warning: Cannot create location.txt - one or both FASTA files missing")
            
    except Exception as e:
        print(f"Error during file operations: {str(e)}")
    
    # Clean up
    if os.path.exists("downloaded_genomes.gb"):
        os.remove("downloaded_genomes.gb")
    
    print("\nProcessing complete!")
    print(f"- FASTA files saved in FASTA directory")
    print(f"- Summary information saved to {summary_filename}")
    if beast_choice == '2':
        print(f"- BEAST files saved in BEAST_input directory")

def cli_main():
    """Entry point for command line interface"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Filter and analyze Ebola virus sequences')
    parser.add_argument('-output_loc', type=str, required=True, 
                       help='Output directory for results (use "." for current directory)')
    
    args = parser.parse_args()
    
    # Convert '.' to current directory absolute path
    output_dir = os.path.abspath(args.output_loc)
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Change to output directory
    original_dir = os.getcwd()
    os.chdir(output_dir)
    
    try:
        # Now run the interactive main function
        main()
    finally:
        # Change back to original directory
        os.chdir(original_dir)

def download_sequences(virus_choice, genome_choice, host_choice, metadata_choice, completeness_threshold=0):
    """Download sequences from GenBank"""
    Entrez.email = "anonymous@example.com"
    
    # Build query with all filters
    query_parts = []
    
    # Add virus species filter with variants
    if virus_choice == "Zaire_ebolavirus":
        query_parts.append('("Zaire ebolavirus"[organism] OR "Zaire Ebolavirus"[organism] OR "ZEBOV"[All Fields])')
    else:
        query_parts.append(f'"{virus_choice}"[organism]')
    
    # Add host filter with multiple possible formats
    if host_choice == '1':
        query_parts.append('("Homo sapiens"[host] OR human[host] OR "Homo sapiens"[Source Host] OR Human[Source Host])')
    elif host_choice == '2':
        query_parts.append('("Pan troglodytes"[host] OR chimpanzee[host])')
    elif host_choice == '3':
        query_parts.append('("Gorilla"[host] OR "Gorilla gorilla"[host])')
    
    # Add genome completeness filter
    if genome_choice == '1':  # Complete genomes only
        query_parts.append('("complete genome"[Title] OR "complete sequence"[Title])')
    elif genome_choice == '2' and completeness_threshold > 0:  # Partial genomes with threshold
        ref_length = get_reference_length(virus_choice.split('_')[0])
        min_length = int(ref_length * (completeness_threshold/100))
        query_parts.append(f'"{min_length}"[SLEN]:"{ref_length}"[SLEN]')
    
    # Combine all query parts
    query = " AND ".join(query_parts)
    
    # Get total count without limit
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=0)
    record = Entrez.read(handle)
    handle.close()
    total_count = int(record["Count"])
    
    print(f"Found {total_count} sequences matching criteria...")
    
    if total_count == 0:
        print("\nNo sequences found with the current filters.")
        print("Try broadening your search criteria:")
        print("- Consider selecting 'All hosts'")
        print("- Try including partial genomes")
        print("- Reduce metadata restrictions")
        return None, query
    
    # Get all IDs without limit
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=total_count)
    record = Entrez.read(handle)
    handle.close()
    
    id_list = record["IdList"]
    
    # Download sequences in batches
    batch_size = 100
    output_file = "downloaded_genomes.gb"
    
    with open(output_file, "w") as out_handle:
        for i in range(0, len(id_list), batch_size):
            batch = id_list[i:i + batch_size]
            print(f"Downloading batch {(i//batch_size)+1} of {(len(id_list)-1)//batch_size + 1}...")
            fetch_handle = Entrez.efetch(db="nucleotide", 
                                       id=batch,
                                       rettype="gb", 
                                       retmode="text")
            out_handle.write(fetch_handle.read())
            fetch_handle.close()
            time.sleep(1)  # Be nice to NCBI servers
    
    return output_file, query

def get_reference_length(virus_type):
    # Reference lengths for different Ebola virus species
    reference_lengths = {
        'Zaire': 18959,
        'Sudan': 18875,
        'Bundibugyo': 18940,
        'Tai Forest': 18935,
        'Reston': 18891
    }
    return reference_lengths.get(virus_type, 18959)

def summarize_locations(records, threshold, virus_choice):
    location_counts = {}
    unknown_count = 0
    oldest_record = None
    oldest_year = 9999
    total_sequences = 0
    
    for record in records:
        # For "All genomes" option, don't apply threshold filter
        if threshold == 0:  # Add this condition for "All genomes" option
            should_include = True
        else:
            seq_length = get_sequence_length(record)
            ref_length = get_reference_length(virus_choice.split()[0])
            completeness = (seq_length / ref_length) * 100
            should_include = completeness >= threshold
        
        if should_include:
            location = get_location(record)
            total_sequences += 1
            
            if location == "Unknown":
                unknown_count += 1
            else:
                # Count locations
                location_counts[location] = location_counts.get(location, 0) + 1
            
            # Check for oldest sequence (include unknown locations)
            year = get_collection_date(record)
            if year != 9999:  # Only consider records with valid dates
                if year < oldest_year:
                    oldest_year = year
                    oldest_record = record
                elif year == oldest_year:  # If same year, prefer DRC sequences (for 1976)
                    if get_location(record) == 'Democratic_Republic_of_the_Congo':
                        oldest_record = record
                        oldest_year = year
    
    return location_counts, unknown_count, oldest_record, oldest_year, total_sequences

def format_record_id(record, virus_choice, metadata_choice):
    """Format record ID based on metadata choice"""
    base_id = record.id
    location = get_location(record)
    date = get_collection_date(record)
    virus_name = virus_choice
    
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

def create_fasta_location_file(fasta_file, outgroup_file):
    """Create location.txt in FASTA directory based on FASTA headers."""
    # First check if both files exist
    if not (os.path.exists(fasta_file) and os.path.exists(outgroup_file)):
        print("Warning: Cannot create location.txt - one or both FASTA files missing")
        return
    
    # Get headers from main FASTA file
    headers = []
    with open(fasta_file) as f:
        for line in f:
            if line.startswith('>'):
                headers.append(line[1:].strip())  # Remove '>' and whitespace
    
    # Get header from outgroup FASTA file
    with open(outgroup_file) as f:
        for line in f:
            if line.startswith('>'):
                headers.append(line[1:].strip())
                break  # Only need first header from outgroup
    
    # Create location.txt
    fasta_dir = os.path.dirname(fasta_file)
    location_file = os.path.join(fasta_dir, "location.txt")
    
    with open(location_file, "w") as out:
        for header in headers:
            # Extract location from header (assuming format ID/species/location or ID/species/location/date)
            parts = header.split('/')
            if len(parts) >= 3:
                location = parts[2]
                out.write(f"{header}\t{location}\n")
    
    print(f"Created location.txt in {fasta_dir}")

if __name__ == "__main__":
    cli_main() 