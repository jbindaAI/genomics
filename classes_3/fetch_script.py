import Bio
from Bio import Entrez, SeqIO
from typing import Union, List
from pathlib import Path

Entrez.email = "XYZ@..."

def read_organisms_from_txt(path2file:str):
    with open(path2file, "r") as f:
        spec_list = []
        for line in f:
            spec_list.append(line.strip())
    return spec_list

def is_full_sequence(seq_record:Bio.SeqRecord.SeqRecord, min_length:int=1400, max_length:int=1600):
    """Checks if a sequence record appears to be full-length based on length and title."""
    if len(seq_record.seq) >= min_length and len(seq_record.seq) <= max_length:
        return True
    return False

def download_16s_rRNA(species_list:Union[List[str], str], output_dir:str="sequences"):
    """Downloads 16S rRNA sequences from NCBI for a specified list of species.
    
    Parameters:
        species (list or str): List of species names (strings) to search for or path to file with species names (each name in new line).
        output_dir (str): Directory to save the output FASTA files.
    """
    import os
    
    # Create output directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if type(species_list) == str:
        species_list = read_organisms_from_txt(species_list)
    
    for species in species_list:
        print(f"Searching for 16S rRNA sequence of {species}...")
        
        # Search for the 16S rRNA sequences for the species in the NCBI nucleotide database
        search_term = f"({species}[Organism]) AND 16s ribosomal RNA"
        
        try:
            with Entrez.esearch(db="nucleotide", term=search_term, retmode="xml", retmax=50) as search_handle:
                search_results = Entrez.read(search_handle)
                id_list = search_results["IdList"]
            
            if not id_list:
                print(f"No 16S rRNA sequence found for {species}.")
                continue

            # Iterate over each found ID to find the first full-length sequence
            sequence_found = False
            for record_id in id_list:
                with Entrez.efetch(db="nucleotide", id=record_id, rettype="fasta", retmode="text") as fetch_handle:
                    seq_record = SeqIO.read(fetch_handle, "fasta")
                
                # Check if the sequence appears to be full-length
                if is_full_sequence(seq_record):
                    # Save the sequence to a FASTA file
                    filename = f"{output_dir}/{species.replace(' ', '_')}_16S_rRNA.fasta"
                    SeqIO.write(seq_record, filename, "fasta")
                    print(f"Downloaded full-length 16S rRNA sequence for {species} and saved as {filename}.")
                    sequence_found = True
                    break
            
            if not sequence_found:
                print(f"No full-length 16S rRNA sequence found for {species}.")
        
        except Exception as e:
            print(f"An error occurred while downloading data for {species}: {e}")

def merge_sequences_from_dir(dir_path: str, output_filename="merged.fasta"):
    dir_path = Path(dir_path)

    with open(output_filename, "w") as output_file:
        for fasta_file in dir_path.glob("*.fasta"):
            try:
                # Parse each FASTA file and write each sequence to the merged output
                records = SeqIO.parse(fasta_file, "fasta")
                SeqIO.write(records, output_file, "fasta")
                print(f"Merged: {fasta_file.name}")
            except Exception as e:
                print(f"Error processing {fasta_file.name}: {e}")
                
    print(f"All sequences merged into {output_filename}")

if __name__ == "__main__":
    species_list = "/home/jakub/genomics/classes3/seq2fetch.txt"
    #species_list = ["Halobacterium salinarum"]
    download_16s_rRNA(species_list)
    merge_sequences_from_dir("/home/jakub/genomics/classes3/sequences")
