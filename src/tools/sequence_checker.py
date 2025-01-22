# src/tools/sequence_checker.py
from Bio import SeqIO
from collections import defaultdict
import os
import hashlib

def check_sequences(fasta_dir: str):
    """Check for duplicate sequences and generate sequence statistics."""
    sequence_hashes = defaultdict(list)
    sequence_stats = defaultdict(dict)
    
    # Process each FASTA file
    for fasta_file in os.listdir(fasta_dir):
        if not fasta_file.endswith(('.fa', '.fasta', '.fna')):
            continue
            
        file_path = os.path.join(fasta_dir, fasta_file)
        for record in SeqIO.parse(file_path, "fasta"):
            # Generate hash of sequence
            seq_hash = hashlib.md5(str(record.seq).encode()).hexdigest()
            sequence_hashes[seq_hash].append(record.id)
            
            # Collect statistics
            sequence_stats[record.id] = {
                'length': len(record.seq),
                'gc_content': (record.seq.count('G') + record.seq.count('C')) / len(record.seq),
                'file': fasta_file
            }
    
    # Print results
    print("\n=== Sequence Analysis ===")
    print(f"Total unique sequences: {len(sequence_hashes)}")
    
    # Check for duplicates
    duplicates = {k: v for k, v in sequence_hashes.items() if len(v) > 1}
    if duplicates:
        print("\nDuplicate sequences found:")
        for hash_val, ids in duplicates.items():
            print(f"Identical sequences: {', '.join(ids)}")
    
    # Print statistics
    print("\nSequence Statistics:")
    for seq_id, stats in sequence_stats.items():
        print(f"\n{seq_id}:")
        print(f"  Length: {stats['length']:,} bp")
        print(f"  GC Content: {stats['gc_content']:.2%}")
        print(f"  File: {stats['file']}")

if __name__ == "__main__":
    check_sequences("data/fasta")
