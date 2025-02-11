def download_sequences(virus_choice, genome_choice, host_choice, metadata_choice):
    """Download sequences from GenBank"""
    Entrez.email = "anonymous@example.com"
    
    # Build query with all filters
    query_parts = []
    
    # Add virus species filter with variants
    if virus_choice == "Zaire_ebolavirus":
        query_parts.append('("Zaire ebolavirus"[organism] OR "Zaire Ebolavirus"[organism] OR "ZEBOV"[All Fields])')
    else:
        query_parts.append(f'"{virus_choice}"[organism]')
    
    # Add host filter with strict matching
    if host_choice == '1':  # Human
        query_parts.append('("Homo sapiens"[host] OR "human"[host]) NOT ("macaque"[host] OR "monkey"[host] OR "Pan troglodytes"[host] OR "Gorilla"[host] OR "bat"[host] OR "Unknown"[host] OR "cynomolgus macaque"[host] OR "cynomolgus_macaque"[host] OR "cynomolgusmacaque"[host] OR "Macaca fascicularis"[host] OR "Macaca mulatta"[host] OR "rhesus"[host])')
    elif host_choice == '2':  # Chimpanzee
        query_parts.append('"Pan troglodytes"[host]')
    elif host_choice == '3':  # Gorilla
        query_parts.append('"Gorilla"[host]')
    elif host_choice == '4':  # All hosts
        pass  # No host filter
    
    # Add other primates option (new)
    elif host_choice == '5':  # Other primates
        query_parts.append('(primate[host] OR macaque[host] OR monkey[host]) NOT ("Homo sapiens"[host] OR "human"[host] OR "Pan troglodytes"[host] OR "Gorilla"[host])')
    
    # ... rest of the function remains the same ... 