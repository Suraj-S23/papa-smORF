from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import mmap
from typing import Dict, List, Tuple, Any
from pathlib import Path
import logging
from src.preprocessing.utils import setup_logger, MMapHandler


class SequenceProcessor:
    def __init__(self):
        self.logger = setup_logger('sequence_processor')
        self.mmap_handler = MMapHandler()
        self.start_codons = ['ATG']
        self.stop_codons = ['TAA', 'TAG', 'TGA']

    def parse_fasta_mmap(self, fasta_file: str) -> Dict[str, str]:
        """Parse FASTA file using memory mapping"""
        sequences = {}
        try:
# src/preprocessing/sequence_processor.py (continued)
            with open(fasta_file, 'rb') as f:
                mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
                
                current_seq = []
                current_id = None
                
                for line in iter(mm.readline, b''):
                    line = line.decode().strip()
                    if line.startswith('>'):
                        if current_id:
                            sequences[current_id] = ''.join(current_seq)
                        current_id = line[1:].split()[0]
                        current_seq = []
                    else:
                        current_seq.append(line)
                
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                
                mm.close()
                self.logger.info(f"Successfully parsed {fasta_file}")
                
            return sequences
        except Exception as e:
            self.logger.error(f"Error parsing FASTA file {fasta_file}: {str(e)}")
            raise

    def parse_gff_mmap(self, gff_file: str) -> Dict[str, List[Dict[str, Any]]]:
        """Parse GFF file using memory mapping"""
        annotations = {}
        try:
            with open(gff_file, 'rb') as f:
                mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
                
                for line in iter(mm.readline, b''):
                    line = line.decode().strip()
                    if line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) >= 8:
                        seqid = parts[0]
                        if seqid not in annotations:
                            annotations[seqid] = []
                        
                        annotations[seqid].append({
                            'start': int(parts[3]),
                            'end': int(parts[4]),
                            'strand': parts[6],
                            'type': parts[2],
                            'attributes': parts[8] if len(parts) > 8 else ''
                        })
                
                mm.close()
                self.logger.info(f"Successfully parsed {gff_file}")
                
            return annotations
        except Exception as e:
            self.logger.error(f"Error parsing GFF file {gff_file}: {str(e)}")
            raise

    def find_orfs(self, sequence: str) -> Tuple[List[Tuple[int, int]], List[int]]:
        """Find all possible ORFs in a sequence"""
        sequence = sequence.upper()
        seq_length = len(sequence)
        orf_positions = []
        start_codon_positions = []

        try:
            # Look for ORFs in all three frames
            for frame in range(3):
                frame_start = frame
                while frame_start < seq_length - 2:
                    # Check for start codon
                    codon = sequence[frame_start:frame_start + 3]
                    if codon in self.start_codons:
                        start_codon_positions.append(frame_start)
                        
                        # Look for stop codon
                        orf_start = frame_start
                        current_pos = orf_start + 3
                        
                        while current_pos < seq_length - 2:
                            current_codon = sequence[current_pos:current_pos + 3]
                            if current_codon in self.stop_codons:
                                # Found complete ORF
                                orf_positions.append((orf_start, current_pos + 3))
                                break
                            current_pos += 3
                            
                    frame_start += 3

            return orf_positions, start_codon_positions
            
        except Exception as e:
            self.logger.error(f"Error finding ORFs: {str(e)}")
            raise

    def process_sequences(self, sequences: Dict[str, str], 
                        annotations: Dict[str, List[Dict]]) -> List[Dict[str, Any]]:
        """Process sequences and store results using memory mapping"""
        try:
            all_results = []
            
            for seq_id, sequence in sequences.items():
                self.logger.info(f"Processing sequence {seq_id}")
                
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
                    'chunk_id': chunk_id,
                    'annotations': annotations.get(seq_id, [])
                }
                
                all_results.append(processed_data)
                
                self.logger.info(f"Processed sequence {seq_id}: "
                               f"Found {len(orf_positions)} ORFs, "
                               f"{len(start_codon_positions)} start codons")
                
            return all_results
            
        except Exception as e:
            self.logger.error(f"Error processing sequences: {str(e)}")
            raise

    def validate_sequence(self, sequence: str) -> bool:
        """Validate sequence content"""
        valid_bases = set('ATCGN')
        sequence = sequence.upper()
        return all(base in valid_bases for base in sequence)

    def get_sequence_stats(self, sequence: str) -> Dict[str, float]:
        """Calculate sequence statistics"""
        sequence = sequence.upper()
        length = len(sequence)
        gc_content = (sequence.count('G') + sequence.count('C')) / length * 100
        
        return {
            'length': length,
            'gc_content': gc_content,
            'n_content': sequence.count('N') / length * 100
        }

    def get_orf_stats(self, orf_positions: List[Tuple[int, int]]) -> Dict[str, float]:
        """Calculate ORF statistics"""
        orf_lengths = [end - start for start, end in orf_positions]
        
        return {
            'total_orfs': len(orf_positions),
            'mean_length': sum(orf_lengths) / len(orf_lengths) if orf_lengths else 0,
            'min_length': min(orf_lengths) if orf_lengths else 0,
            'max_length': max(orf_lengths) if orf_lengths else 0
        }
