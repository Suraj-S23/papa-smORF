import os
import struct
from collections import defaultdict
import glob
import pandas as pd
from .utils import setup_logger, load_config

logger = setup_logger()

class SequenceProcessor:
    def __init__(self, config_path='configs/config.yaml'):
        self.config = load_config(config_path)
        self.chunk_size = self.config['preprocessing']['chunk_size']
        
    def normalize_sequence_id(self, seq_id):
        """Normalize sequence ID."""
        seq_id = seq_id.split()[0]
        seq_id = seq_id.split('.')[0]
        return seq_id.strip()

    def parse_multiple_fasta(self, fasta_dir):
        """Parse multiple FASTA files."""
        sequences = {}
        fasta_files = glob.glob(os.path.join(fasta_dir, "*.fna*"))
        
        if not fasta_files:
            raise ValueError(f"No FASTA files found in {fasta_dir}")
        
        logger.info(f"Found {len(fasta_files)} FASTA files")
        
        for fasta_path in fasta_files:
            logger.info(f"Parsing FASTA file: {os.path.basename(fasta_path)}")
            with open(fasta_path, 'r') as fasta_file:
                current_id = None
                current_seq = []

                for line in fasta_file:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_id:
                            sequences[current_id] = {
                                'sequence': ''.join(current_seq),
                                'original_id': current_id,
                                'normalized_id': self.normalize_sequence_id(current_id)
                            }
                        current_id = line[1:]
                        current_seq = []
                    else:
                        current_seq.append(line)

                if current_id:
                    sequences[current_id] = {
                        'sequence': ''.join(current_seq),
                        'original_id': current_id,
                        'normalized_id': self.normalize_sequence_id(current_id)
                    }
        
        return sequences

    def parse_multiple_gff(self, gff_dir):
        """Parse multiple GFF files."""
        annotations = defaultdict(lambda: {"orf_positions": [], "start_codon_positions": []})
        gff_files = glob.glob(os.path.join(gff_dir, "*.gff*"))
        
        if not gff_files:
            raise ValueError(f"No GFF files found in {gff_dir}")
        
        logger.info(f"Found {len(gff_files)} GFF files")
        
        for gff_path in gff_files:
            logger.info(f"Parsing GFF file: {os.path.basename(gff_path)}")
            with open(gff_path, 'r') as gff_file:
                for line in gff_file:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) < 9:
                        continue

                    seq_id = parts[0]
                    normalized_id = self.normalize_sequence_id(seq_id)
                    feature_type = parts[2]
                    start = int(parts[3])
                    end = int(parts[4])

                    if feature_type == 'CDS':
                        annotations[normalized_id]["orf_positions"].append((start, end))
                        annotations[normalized_id]["start_codon_positions"].append(start)

        return annotations

    def create_target_vector(self, sequence_length, orf_positions):
        """Create target vector indicating ORF positions."""
        target = [0] * sequence_length
        for start, end in orf_positions:
            start_idx = max(0, start - 1)
            end_idx = min(sequence_length - 1, end - 1)
            
            if start_idx < sequence_length and end_idx >= 0:
                for i in range(start_idx, end_idx + 1):
                    if 0 <= i < sequence_length:
                        target[i] = 1
        
        return target

    def chunk_sequence(self, sequence):
        """Split sequence into chunks."""
        return [sequence[i:i + self.chunk_size] 
                for i in range(0, len(sequence), self.chunk_size)]

    def get_chunk_features(self, chunk_start, chunk_end, orf_positions, 
                          start_codon_positions, target_vector):
        """Get features for a specific chunk."""
        chunk_orfs = [(max(start, chunk_start+1) - chunk_start, 
                       min(end, chunk_end) - chunk_start)
                      for start, end in orf_positions
                      if start <= chunk_end and end >= chunk_start]
        
        chunk_start_codons = [pos - chunk_start 
                             for pos in start_codon_positions 
                             if chunk_start < pos <= chunk_end]
        
        chunk_target = target_vector[chunk_start:chunk_end]
        
        return chunk_orfs, chunk_start_codons, chunk_target

    def process_sequences(self, sequences, annotations, output_dir):
        """Process sequences and write to output directory."""
        os.makedirs(output_dir, exist_ok=True)
        
        # Create directories
        seq_dir = os.path.join(output_dir, 'sequences')
        meta_dir = os.path.join(output_dir, 'metadata')
        os.makedirs(seq_dir, exist_ok=True)
        os.makedirs(meta_dir, exist_ok=True)

        # Create chunk index
        chunk_index = []

        for orig_id, seq_data in sequences.items():
            normalized_id = seq_data['normalized_id']
            if normalized_id not in annotations:
                continue
                
            sequence = seq_data['sequence']
            
            try:
                sample_id = normalized_id
                safe_filename = "".join(c for c in normalized_id if c.isalnum() or c in ('_', '-'))
                
                seq_file_path = os.path.join(seq_dir, f'{safe_filename}.bin')
                meta_file_path = os.path.join(meta_dir, f'{safe_filename}.meta')

                chunks = self.chunk_sequence(sequence)
                target_vector = self.create_target_vector(
                    len(sequence), 
                    annotations[normalized_id]["orf_positions"]
                )
                
                with open(seq_file_path, 'wb') as seq_file, open(meta_file_path, 'wb') as meta_file:
                    seq_offset = seq_file.tell()
                    
                    for chunk_idx, chunk in enumerate(chunks):
                        chunk_start = chunk_idx * self.chunk_size
                        chunk_end = chunk_start + len(chunk)
                        
                        chunk_orfs, chunk_start_codons, chunk_target = self.get_chunk_features(
                            chunk_start, 
                            chunk_end,
                            annotations[normalized_id]["orf_positions"],
                            annotations[normalized_id]["start_codon_positions"],
                            target_vector
                        )
                        
                        chunk_id = f"{sample_id}_{chunk_idx}"
                        
                        chunk_index.append({
                            'sample_id': sample_id,
                            'chunk_id': chunk_id
                        })
                        
                        seq_file.write(chunk.encode('utf-8'))
                    
                    seq_file.write(b'\n')

                    # Write metadata
                    meta_file.write(struct.pack('I', seq_offset))
                    meta_file.write(struct.pack('I', len(sequence)))
                    meta_file.write(struct.pack('I', len(chunks)))
                    
                    # Write ORF positions
                    orf_positions = annotations[normalized_id]["orf_positions"]
                    meta_file.write(struct.pack('I', len(orf_positions)))
                    for start, end in orf_positions:
                        meta_file.write(struct.pack('II', start, end))
                    
                    # Write start codon positions
                    start_positions = annotations[normalized_id]["start_codon_positions"]
                    meta_file.write(struct.pack('I', len(start_positions)))
                    for pos in start_positions:
                        meta_file.write(struct.pack('I', pos))
                    
                    # Write target vector
                    meta_file.write(struct.pack('I', len(target_vector)))
                    meta_file.write(bytes(target_vector))
                
                logger.info(f"Processed: {orig_id}")
                
            except Exception as e:
                logger.error(f"Error processing sequence {orig_id}: {str(e)}")
                continue

        # Save chunk index
        chunk_index_df = pd.DataFrame(chunk_index)
        chunk_index_df.to_csv(os.path.join(output_dir, 'chunk_index.csv'), index=False)
        
        return chunk_index_df

def load_chunk(self, chunk_id, data_dir):
    """Load a specific chunk of sequence data."""
    # Read chunk index
    chunk_index_df = pd.read_csv(os.path.join(data_dir, 'chunk_index.csv'))
    
    # Find the specific chunk
    chunk_info = chunk_index_df[chunk_index_df['chunk_id'] == chunk_id].iloc[0]
    sample_id = chunk_info['sample_id']
    
    # Construct file paths
    safe_filename = "".join(c for c in sample_id if c.isalnum() or c in ('_', '-'))
    seq_file = os.path.join(data_dir, 'sequences', f'{safe_filename}.bin')
    meta_file = os.path.join(data_dir, 'metadata', f'{safe_filename}.meta')
    
    with open(seq_file, 'rb') as sf, open(meta_file, 'rb') as mf:
        # Read metadata
        seq_offset = struct.unpack('I', mf.read(4))[0]
        seq_length = struct.unpack('I', mf.read(4))[0]
        num_chunks = struct.unpack('I', mf.read(4))[0]
        
        # Get chunk number from chunk_id
        chunk_number = int(chunk_id.split('_')[-1])
        chunk_start = chunk_number * self.chunk_size
        
        # Read ORF positions
        num_orfs = struct.unpack('I', mf.read(4))[0]
        orf_positions = []
        for _ in range(num_orfs):
            start, end = struct.unpack('II', mf.read(8))
            orf_positions.append((start, end))
        
        # Read start codon positions
        num_start_codons = struct.unpack('I', mf.read(4))[0]
        start_codon_positions = []
        for _ in range(num_start_codons):
            pos = struct.unpack('I', mf.read(4))[0]
            start_codon_positions.append(pos)
        
        # Read target vector
        target_length = struct.unpack('I', mf.read(4))[0]
        target_vector = list(mf.read(target_length))
        
        # Read sequence chunk
        chunk_end = min(chunk_start + self.chunk_size, seq_length)
        sf.seek(seq_offset + chunk_start)
        chunk_sequence = sf.read(chunk_end - chunk_start).decode('utf-8')
        
        # Get chunk-specific features
        chunk_orfs, chunk_start_codons, chunk_target = self.get_chunk_features(
            chunk_start, chunk_end, orf_positions, start_codon_positions, target_vector
        )
        
        return {
            'sequence': chunk_sequence,
            'sample_id': sample_id,
            'chunk_id': chunk_id,
            'target': chunk_target,
            'orf_positions': chunk_orfs,
            'start_codon_positions': chunk_start_codons
        }

