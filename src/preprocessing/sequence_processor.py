# src/preprocessing/sequence_processor.py
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
from pathlib import Path
import pickle
import logging
from typing import List, Dict, Tuple, Any
import gzip
from src.preprocessing.utils import setup_logger

class SequenceProcessor:
    def __init__(self):
        self.logger = setup_logger('papa_smorf')
        self.chunk_size = 5000  # Default chunk size, adjust if needed

    def parse_multiple_fasta(self, fasta_files: List[str]) -> Dict[str, SeqRecord]:
        """Parse multiple FASTA files and return sequences."""
        sequences = {}
        for fasta_file in fasta_files:
            try:
                # Handle both .gz and regular files
                if fasta_file.endswith('.gz'):
                    with gzip.open(fasta_file, 'rt') as handle:
                        for record in SeqIO.parse(handle, 'fasta'):
                            sequences[record.id] = record
                else:
                    for record in SeqIO.parse(fasta_file, 'fasta'):
                        sequences[record.id] = record
                        
                self.logger.info(f"Successfully parsed {fasta_file}")
            except Exception as e:
                self.logger.error(f"Error parsing {fasta_file}: {str(e)}")
                raise
                
        return sequences

    def parse_multiple_gff(self, gff_files: List[str]) -> Dict[str, List[Dict]]:
        """Parse multiple GFF files and return annotations."""
        annotations = {}
        for gff_file in gff_files:
            try:
                current_annotations = []
                # Handle both .gz and regular files
                opener = gzip.open if gff_file.endswith('.gz') else open
                with opener(gff_file, 'rt') as f:
                    for line in f:
                        if line.startswith('#'):
                            continue
                        parts = line.strip().split('\t')
                        if len(parts) >= 8:  # Valid GFF line
                            annotation = {
                                'seqid': parts[0],
                                'source': parts[1],
                                'type': parts[2],
                                'start': int(parts[3]),
                                'end': int(parts[4]),
                                'score': parts[5],
                                'strand': parts[6],
                                'phase': parts[7],
                                'attributes': parts[8] if len(parts) > 8 else ''
                            }
                            current_annotations.append(annotation)
                
                # Group annotations by sequence ID
                for ann in current_annotations:
                    seqid = ann['seqid']
                    if seqid not in annotations:
                        annotations[seqid] = []
                    annotations[seqid].append(ann)
                
                self.logger.info(f"Successfully parsed {gff_file}")
            except Exception as e:
                self.logger.error(f"Error parsing {gff_file}: {str(e)}")
                raise
                
        return annotations

    def find_orfs(self, sequence: str) -> Tuple[List[Tuple[int, int]], List[int]]:
        """Find all possible ORFs in a sequence."""
        start_codons = ['ATG']
        stop_codons = ['TAA', 'TAG', 'TGA']
        
        orf_positions = []
        start_positions = []
        seq_len = len(sequence)
        
        for i in range(0, seq_len - 2):
            codon = sequence[i:i+3]
            if codon in start_codons:
                start_positions.append(i)
                # Look for stop codon
                for j in range(i + 3, seq_len - 2, 3):
                    if sequence[j:j+3] in stop_codons:
                        orf_positions.append((i, j + 2))
                        break
        
        return orf_positions, start_positions

    def process_sequences(self, sequences: Dict[str, SeqRecord], 
                        annotations: Dict[str, List[Dict]]) -> List[Dict[str, Any]]:
        """Process sequences and return structured data."""
        try:
            all_results = []
            
            for seq_id, record in sequences.items():
                sequence = str(record.seq)
                
                # Find ORFs
                orf_positions, start_codon_positions = self.find_orfs(sequence)
                
                # Create chunk ID based on sequence ID
                chunk_id = f"{seq_id}_chunk"
                
                # Store processed data
                processed_data = {
                    'sequence': sequence,
                    'orf_positions': orf_positions,
                    'start_codon_positions': start_codon_positions,
                    'sample_id': seq_id,
                    'chunk_id': chunk_id
                }
                
                all_results.append(processed_data)
                
                self.logger.info(f"Processed sequence {seq_id}: "
                               f"Found {len(orf_positions)} ORFs, "
                               f"{len(start_codon_positions)} start codons")
                
            return all_results
            
        except Exception as e:
            self.logger.error(f"Error processing sequences: {str(e)}")
            raise

    def save_chunk(self, chunk_data: Dict[str, Any], output_dir: str) -> None:
        """Save processed chunk to file."""
        try:
            chunk_id = chunk_data['chunk_id']
            output_path = os.path.join(output_dir, f"{chunk_id}.pkl")
            
            with open(output_path, 'wb') as f:
                pickle.dump(chunk_data, f)
                
            self.logger.info(f"Saved chunk {chunk_id} to {output_path}")
            
        except Exception as e:
            self.logger.error(f"Error saving chunk {chunk_id}: {str(e)}")
            raise

    def load_chunk(self, chunk_id: str, processed_dir: str) -> Dict[str, Any]:
        """Load processed chunk from file."""
        try:
            chunk_path = os.path.join(processed_dir, f"{chunk_id}.pkl")
            
            with open(chunk_path, 'rb') as f:
                chunk_data = pickle.load(f)
                
            return chunk_data
            
        except Exception as e:
            self.logger.error(f"Error loading chunk {chunk_id}: {str(e)}")
            raise

    def validate_sequence(self, sequence: str) -> bool:
        """Validate sequence contains only valid nucleotides."""
        valid_nucleotides = set('ATCGN')
        return all(nuc in valid_nucleotides for nuc in sequence.upper())

    def validate_annotations(self, annotations: List[Dict], sequence_length: int) -> bool:
        """Validate annotations are within sequence bounds."""
        for ann in annotations:
            if ann['start'] < 1 or ann['end'] > sequence_length:
                return False
        return True
