import os
import sys
from Bio import SeqIO
import numpy as np
import h5py
import pandas as pd

def parse_fasta(fasta_file):
    sequences = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
    return sequences

def parse_gff(gff_file):
    orf_positions = []
    with open(gff_file, "r") as file:
        for line in file:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                if len(fields) >= 7 and fields[2] == "CDS":
                    start = int(fields[3])
                    stop = int(fields[4])
                    orf_positions.append((start, stop))
    return orf_positions

def process_data(sequences, orf_positions, output_file):
    with h5py.File(output_file, "w") as hdf5_file:
        for i, sequence in enumerate(sequences):
            seq_length = len(sequence)
            labels = np.zeros(seq_length, dtype=int)
            for start, stop in orf_positions[i]:
                labels[start-1:stop] = 1

            seq_group = hdf5_file.create_group(f"sequence_{i}")
            seq_group.create_dataset("sequence", data=sequence)
            seq_group.create_dataset("orf_start_stop", data=np.array(orf_positions[i]))
            seq_group.create_dataset("labels", data=labels)

            print(f"\nProcessing sequence {i+1}:")
            print(f"Sequence length: {seq_length}")
            print(f"Number of ORFs: {len(orf_positions[i])}")
            
            print("\nGroup structure:")
            print(f"- sequence_{i}")
            print(f"  - sequence (dataset)")
            print(f"  - orf_start_stop (dataset)")
            print(f"  - labels (dataset)")

def display_sample_data(hdf5_file_path):
    with h5py.File(hdf5_file_path, 'r') as f:
        print(f"\nFinal HDF5 file structure:")
        for key in f.keys():
            group = f[key]
            print(f"Group: {key}")
            for sub_key in group.keys():
                dataset = group[sub_key]
                print(f"  Dataset: {sub_key}")
                print(f"    Type: {type(dataset)}")
                print(f"    Shape: {dataset.shape}")
                
                # Display sample data if possible
                try:
                    data = dataset[()]
                    if isinstance(data, np.ndarray):
                        if sub_key == 'orf_start_stop':
                            sample_size = min(5, len(data))  # Take up to 5 samples for orf_start_stop
                            sample_data = data[:sample_size]
                            print(f"    Sample data ({sample_size} rows):")
                            print(sample_data)
                        elif sub_key == 'labels':
                            print(f"    Unique labels: {np.unique(data)}")
                            print(f"    Label counts: {np.bincount(data)}")
                        else:
                            print("    Data type is a NumPy array.")
                    elif isinstance(data, str):
                        seq_length = len(data)
                        snippet_size = 10
                        snippet = data[:snippet_size] + "..." if seq_length > snippet_size else data
                        print(f"    Sequence length: {seq_length}")
                        print(f"    Sequence snippet: {snippet}")
                    else:
                        print("    Data type is not supported.")
                except Exception as e:
                    print(f"    Error accessing dataset value: {str(e)}")
                
                print()

def main():
    # Define the expected directory structure
    base_dir = os.path.dirname(os.getcwd())
    
    # Construct file paths
    fasta_dir = os.path.join(base_dir, "Data", "Raw", "small_data", "small_fasta")
    gff_dir = os.path.join(base_dir, "Data", "Raw", "small_data", "small_gff")
    output_file = os.path.join(base_dir, "data", "output", "processed_data.h5")

    print(f"Base Directory: {base_dir}")
    print(f"FASTA Directory: {fasta_dir}")
    print(f"GFF Directory: {gff_dir}")
    print(f"Output File: {output_file}")

    if not os.path.exists(output_file):
        try:
            # Check if directories exist
            if not os.path.exists(fasta_dir):
                raise FileNotFoundError(f"FASTA directory does not exist: {fasta_dir}")
            
            if not os.path.exists(gff_dir):
                raise FileNotFoundError(f"GFF directory does not exist: {gff_dir}")
            
            # List contents of directories for debugging
            print("\nContents of FASTA directory:")
            for item in os.listdir(fasta_dir):
                print(item)
            
            print("\nContents of GFF directory:")
            for item in os.listdir(gff_dir):
                print(item)

            fasta_files = [os.path.join(fasta_dir, file) for file in os.listdir(fasta_dir)]
            gff_files = [os.path.join(gff_dir, file) for file in os.listdir(gff_dir)]

            print("\nNumber of FASTA files:", len(fasta_files))
            print("Number of GFF files:", len(gff_files))

            sequences = []
            orf_positions = []
            for i, fasta_file in enumerate(fasta_files):
                base_name = os.path.basename(fasta_file).replace('.fna', '.gff')
                
                matching_gff_file = next((file for file in gff_files if file.endswith(base_name)), None)
                
                print(f"Processing FASTA file {i+1}: {fasta_file}")
                print(f"Matching GFF file: {matching_gff_file}")
                
                try:
                    sequence = parse_fasta(fasta_file)
                    sequences.append(sequence)
                    
                    if matching_gff_file:
                        orf_position = parse_gff(matching_gff_file)
                        orf_positions.append(orf_position)
                    else:
                        print(f"No matching GFF file found for {fasta_file}")
                except Exception as e:
                    print(f"Error processing files {i+1}: {e}")

            print(f"Number of sequences processed: {len(sequences)}")
            
            process_data(sequences, orf_positions, output_file)
            
            print("\nFinal HDF5 file structure:")
            display_sample_data(output_file)
            
            print("Data processing completed successfully.")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
    else:
        print("Output file already exists.")
        display_sample_data(output_file)

if __name__ == "__main__":
    main()