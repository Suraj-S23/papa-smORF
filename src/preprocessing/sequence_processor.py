# src/preprocessing/sequence_processor.py
from Bio import SeqIO
from Bio.Seq import Seq
import h5py
import pandas as pd
import numpy as np
import gzip
import logging
from typing import Dict, List, Generator, Tuple, Any, Optional
from collections import defaultdict
from pathlib import Path
import yaml
from dataclasses import dataclass
import hashlib
from .utils import setup_logger

@dataclass
class ORFMetadata:
    original_id: str
    seqid: str
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive
    strand: str
    attributes: dict = None
    def uuid(self) -> str:
        """Enhanced UUID with protein_id and original sequence"""
        base = f"{self.seqid}|{self.start}|{self.end}|{self.strand}|{self.original_id}"
        return hashlib.sha256(base.encode()).hexdigest()
    def reverse_complement(self, parent_sequence: str) -> 'ProcessedORF':
        """Correct reverse complement of only the ORF subsequence"""
        orf_subseq = parent_sequence[self.start-1:self.end]  # 0-based slice
        rc_subseq = str(Seq(orf_subseq).reverse_complement())
        return ProcessedORF(
            orf_id=self.uuid(),
            sequence=rc_subseq,
            original_strand=self.strand,
            genomic_start=self.start,
            genomic_end=self.end,
            processed_start=1,
            processed_end=len(rc_subseq)
        )

@dataclass
class ProcessedORF:
    orf_id: str
    sequence: str
    original_strand: str
    genomic_start: int
    genomic_end: int
    processed_start: int
    processed_end: int

class SequenceProcessor:
    def __init__(self, config_path: str = None, chunk_size=10000):
        self.logger = setup_logger('sequence_processor')
        self.config = self.load_config(config_path)
        self.chunk_size = chunk_size

    def load_config(self, config_path: Optional[str]) -> dict:
        """Load configuration with memory-efficient YAML parsing"""
        if config_path:
            with open(config_path, "r") as f:
                return yaml.safe_load(f)
        return {}

    def parse_fasta_generator(self, fasta_file: str) -> Generator[Tuple[str, str], None, None]:
        """Memory-efficient FASTA parsing with gzip support"""
        try:
            handle = gzip.open(fasta_file, 'rt') if fasta_file.endswith('.gz') else open(fasta_file)
            for record in SeqIO.parse(handle, "fasta"):
                yield record.id, str(record.seq).upper()
            handle.close()
        except Exception as e:
            self.logger.error(f"FASTA error: {str(e)}")
            raise

    def parse_gff_generator(self, gff_file: str) -> Generator[ORFMetadata, None, None]:
        """Filter ORFs by seqid during parsing"""
        try:
            compression = 'gzip' if gff_file.endswith('.gz') else None
            for chunk in pd.read_csv(
                gff_file,
                sep='\t',
                comment='#',
                header=None,
                usecols=[0, 2, 3, 4, 6, 8],  # seqid, type, start, end, strand, attributes
                dtype={0: 'category', 2: 'category', 6: 'category'},
                chunksize=10000,
                compression=compression
            ):
                chunk.columns = ['seqid', 'type', 'start', 'end', 'strand', 'attributes']
                for _, row in chunk.iterrows():
                    if row.type.upper() != 'CDS':
                        continue
                    attributes = self._parse_attributes(row.attributes)
                    yield ORFMetadata(
                        original_id=attributes.get('ID', ''),
                        seqid=row.seqid,
                        start=int(row.start),
                        end=int(row.end),
                        strand=row.strand,
                        attributes=attributes
                    )
        except Exception as e:
            self.logger.error(f"GFF error: {str(e)}")
            raise

    def _parse_attributes(self, attributes: str) -> dict:
        """Parse GFF attributes into a dictionary"""
        if not attributes:
            return {}
        return dict(pair.strip().split('=') for pair in attributes.split(';') if '=' in pair)

    def create_orf_regions(self, sequence: str, orfs: List[ORFMetadata]) -> Tuple[np.ndarray, dict]:
        """Generate multi-channel target array: ORF presence and overlap count."""
        seq_len = len(sequence)
        target = np.zeros((seq_len, 2), dtype=np.uint8)  # Channel 0: ORF presence, Channel 1: Overlap count
        orf_map = {}

        # First pass: Mark ORF presence
        for orf in orfs:
            start = orf.start - 1  # Convert to 0-based
            end = orf.end  # 1-based inclusive
            target[start:end, 0] = 1  # Mark ORF regions in the first channel

        # Second pass: Count overlapping ORFs at each position
        for orf in orfs:
            start = orf.start - 1
            end = orf.end
            target[start:end, 1] += 1  # Increment overlap count

        # Build ORF map
        for orf in orfs:
            orf_id = orf.uuid()
            if orf.strand == '-':
                processed = orf.reverse_complement(sequence)
                orf_seq = processed.sequence
            else:
                orf_seq = sequence[orf.start-1:orf.end]
            
            orf_map[orf_id] = {
                'strand': orf.strand,
                'genomic_start': orf.start,
                'genomic_end': orf.end,
                'processed_seq': orf_seq,
                'overlap_count': np.max(target[orf.start-1:orf.end, 1]),  # Max overlap count in ORF region
                'parent_orf': self._find_parent_orf(orf, orfs),
                'chunk_id': None,
                'is_continuation': False,
                'previous_chunk_id': None
            }

        return target, orf_map


    def _find_parent_orf(self, current_orf: ORFMetadata, all_orfs: List[ORFMetadata]) -> str:
        """Identify nested ORFs using GFF hierarchy"""
        for orf in all_orfs:
            if (orf.start <= current_orf.start and 
                orf.end >= current_orf.end and 
                orf.uuid() != current_orf.uuid()):
                return orf.uuid()
        return ""

    def chunk_sequence_with_orfs(self, sequence: str, target: np.ndarray, 
                                orf_map: dict, chunk_size: int) -> List[dict]:
        """
        Improved chunk tracking with overlap awareness.
        Now, each ORF is assigned only once across all chunks by using a set
        of assigned ORF IDs.
        """
        seq_len = len(sequence)
        chunk_records = []
        assigned_orfs = set()  # Track ORFs that have already been assigned to a chunk
        previous_chunk_id = None

        for chunk_idx in range(0, seq_len, chunk_size):
            chunk_start = chunk_idx
            chunk_end = min(chunk_start + chunk_size, seq_len)
            chunk_seq = sequence[chunk_start:chunk_end].ljust(chunk_size, 'N')
            chunk_target = target[chunk_start:chunk_end]
            chunk_orfs = {}
            chunk_id = f"chunk{chunk_start//chunk_size:04d}"

            for orf_id, meta in orf_map.items():
                # Skip if this ORF has already been assigned to a previous chunk.
                if orf_id in assigned_orfs:
                    continue

                # Convert genomic start (1-based) to 0-based coordinate.
                orf_start = meta['genomic_start'] - 1  
                orf_end = meta['genomic_end']

                # Only assign the ORF to the chunk if its start is in the current chunk.
                if chunk_start <= orf_start < chunk_end:
                    assigned_orfs.add(orf_id)  # Mark this ORF as assigned.
                    rel_start = orf_start - chunk_start
                    # Even if the ORF extends beyond the current chunk, only include the portion in this chunk.
                    rel_end = min(orf_end, chunk_end) - chunk_start
                    # Mark continuation flag if the ORF extends past the current chunk.
                    is_continuation = orf_end > chunk_end

                    chunk_orfs[orf_id] = {
                        'chunk_rel_start': rel_start,
                        'chunk_rel_end': rel_end,
                        'global_start': orf_start + 1,  # Convert back to 1-based
                        'global_end': orf_end,
                        'strand': meta['strand'],
                        'processed_seq': meta['processed_seq'],
                        'overlap_count': meta['overlap_count'],
                        'parent_orf': meta['parent_orf'],
                        'segment_number': (chunk_start // chunk_size) + 1,
                        'chunk_id': chunk_id,
                        'is_continuation': is_continuation,
                        'previous_chunk_id': previous_chunk_id
                    }
            chunk_records.append({
                'sequence': chunk_seq,
                'target': chunk_target,
                'orfs': chunk_orfs,
                'chunk_start': chunk_start,
                'chunk_end': chunk_end,
                'chunk_id': chunk_id,
                'previous_chunk_id': previous_chunk_id
            })
            previous_chunk_id = chunk_id
        return chunk_records

    def _store_chunk(self, group: h5py.Group, chunk_data: dict):
        """Store multi-channel targets and hierarchy"""
        group.create_dataset('sequence', data=np.array(list(chunk_data['sequence']), dtype='S1'),
                           compression="lzf", shuffle=True)
        # Store multi-channel target
        target_ds = group.create_dataset('target', data=chunk_data['target'].astype(np.uint8),
                                       compression="lzf", shuffle=True)
        target_ds.attrs['channel_names'] = ['orf_presence', 'overlap_count']
        orf_group = group.create_group('orfs')
        for orf_id, meta in chunk_data['orfs'].items():
            orf_ds = orf_group.create_group(orf_id)
            orf_ds.attrs.update({
                'strand': meta['strand'].encode('utf-8'),  # Convert string to bytes
                'global_start': meta['global_start'],
                'global_end': meta['global_end'],
                'processed_seq': meta['processed_seq'].encode('utf-8'),  # Convert string to bytes
                'overlap_count': meta['overlap_count'],
                'parent_orf': meta['parent_orf'].encode('utf-8'),  # Convert string to bytes
                'segment_number': meta['segment_number'],
                'chunk_id': meta['chunk_id'].encode('utf-8'),  # Convert string to bytes
                'is_continuation': meta['is_continuation'],
                'previous_chunk_id': (meta['previous_chunk_id'].encode('utf-8') if meta['previous_chunk_id'] else b'')  # Convert string to bytes
            })
            orf_ds.create_dataset('chunk_coords', 
                                 data=[meta['chunk_rel_start'], meta['chunk_rel_end']],
                                 dtype=np.int32)

    def process_sequences(self, fasta_files: List[str], gff_files: List[str], output_file: str):
        """Seqid filtering and enhanced processing"""
        with h5py.File(output_file, 'w') as hdf:
            for fasta, gff in zip(fasta_files, gff_files):
                self.logger.info(f"Processing {Path(fasta).name}")
                file_stem = Path(fasta).stem
                all_orfs = list(self.parse_gff_generator(gff))
                for seq_id, sequence in self.parse_fasta_generator(fasta):
                    # Critical seqid filtering
                    filtered_orfs = [orf for orf in all_orfs if orf.seqid == seq_id]
                    unique_id = f"{file_stem}_{seq_id}"
                    target, orf_map = self.create_orf_regions(sequence, filtered_orfs)
                    chunks = self.chunk_sequence_with_orfs(sequence, target, orf_map, self.chunk_size)
                    seq_group = hdf.create_group(unique_id)
                    for i, chunk in enumerate(chunks):
                        chunk_id = f"{unique_id}_{chunk['chunk_id']}"
                        chunk_group = seq_group.create_group(chunk_id)
                        self._store_chunk(chunk_group, chunk)