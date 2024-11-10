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

def is_valid_sequence(seq_record:Bio.SeqRecord.SeqRecord, search_term:str, min_length:int=140, max_length:int=180):
    """Checks if a sequence record appears to be full-length based on length and title."""
    if len(seq_record.seq) >= min_length and len(seq_record.seq) <= max_length and search_term in seq_record.description and "partial" not in seq_record.description:
        return True
    return False

def download_prot_sequences(search_term:str, n:int=101, output_dir:str="sequences"):
    """Downloads protein sequences for n different species."""
    import os

    # Create output directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    search_term = "hemoglobin subunit alpha 1[Protein Name]"

    try:
        with Entrez.esearch(db="Protein", term=search_term, retmode="xml", retmax=100) as search_handle:
            search_results = Entrez.read(search_handle)
            id_list = search_results["IdList"]
        
        sequence_found = False
        for record_id in id_list:
            with Entrez.efetch(db="Protein", id=record_id, rettype="fasta", retmode="text") as fetch_handle:
                seq_record = SeqIO.read(fetch_handle, "fasta")
                
            # Check if the sequence appears to be full-length
            if is_valid_sequence(seq_record, search_term):
                # Save the sequence to a FASTA file
                filename = f"{output_dir}/{record_id}.fasta"
                SeqIO.write(seq_record, filename, "fasta")
                print(f"Downloaded full-length sequence for {record_id} and saved as {filename}.")
                sequence_found = True

    except Exception as e:
        print(f"An error occurred while downloading data for query: {search_term}: {e}")
    return


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
    #download_16s_rRNA(species_list)
    #merge_sequences_from_dir("/home/jakub/genomics/classes3/sequences")
    download_prot_sequences(search_term="")

