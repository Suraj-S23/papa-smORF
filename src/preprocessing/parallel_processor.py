# src/preprocessing/parallel_processor.py
import h5py
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import uuid
from tqdm import tqdm
from collections import defaultdict
from typing import List, Dict, Any, Tuple
import os
from .utils import setup_logger, get_file_pairs
from .sequence_processor import SequenceProcessor, ORFMetadata

class ParallelSequenceProcessor:
    def __init__(self, n_workers=None, chunk_size=10000, batch_size=100):
        self.n_workers = n_workers or os.cpu_count()
        self.chunk_size = chunk_size
        self.batch_size = batch_size
        self.processor = SequenceProcessor(chunk_size=chunk_size)
        self.tmp_dir = "tmp_hdf5_chunks"
        self.logger = setup_logger('parallel_processor')
        Path(self.tmp_dir).mkdir(exist_ok=True)

    def _process_single_pair(self, pair: dict) -> dict:
        """Process a single FASTA-GFF pair with ORF tracking"""
        try:
            annotations = defaultdict(list)
            for ann in self.processor.parse_gff_generator(pair['gff']):
                if isinstance(ann, ORFMetadata):
                    annotations[ann.seqid].append(ann)
            chunks_data = defaultdict(list)
            file_stem = Path(pair['fasta']).stem
            for seq_id, seq in self.processor.parse_fasta_generator(pair['fasta']):
                unique_id = f"{file_stem}_{seq_id}"
                seq_ann = [a for a in annotations.get(seq_id, []) 
                           if a.start <= len(seq) and a.end <= len(seq)]
                target, orf_map = self.processor.create_orf_regions(seq, seq_ann)
                chunks = self.processor.chunk_sequence_with_orfs(seq, target, orf_map, self.chunk_size)
                for chunk in chunks:
                    chunk_id = f"{unique_id}_{chunk['chunk_id']}"
                    chunks_data[unique_id].append((chunk_id, chunk))
            return chunks_data
        except Exception as e:
            self.logger.error(f"Error processing {pair['fasta']}: {str(e)}")
            raise

    def process_batch(self, file_pairs: List[dict], output_file: str):
        """Process a batch of files with error handling"""
        batch_id = uuid.uuid4().hex[:8]
        tmp_file = Path(self.tmp_dir) / f"batch_{batch_id}.h5"
        try:
            with h5py.File(tmp_file, 'w') as hdf:
                with ProcessPoolExecutor(max_workers=self.n_workers) as executor:
                    futures = [executor.submit(self._process_single_pair, pair) 
                               for pair in file_pairs]
                    
                    for future in tqdm(futures, desc=f"Batch {batch_id}"):
                        try:
                            chunks_data = future.result()
                            for seq_id, chunks in chunks_data.items():
                                seq_group = hdf.require_group(seq_id)
                                for chunk_id, data in chunks:
                                    chunk_group = seq_group.create_group(chunk_id)
                                    self.processor._store_chunk(chunk_group, data)
                        except Exception as e:
                            self.logger.error(f"Failed to process chunk: {str(e)}")
            
            self._merge_hdf5_chunks(tmp_file, output_file)
            
        finally:
            if tmp_file.exists():
                tmp_file.unlink()

    def _merge_hdf5_chunks(self, src_file: Path, dest_file: str):
        """Atomic HDF5 merging with progress tracking"""
        try:
            with h5py.File(src_file, 'r') as src, h5py.File(dest_file, 'a') as dest:
                for seq_id in tqdm(src.keys(), desc="Merging"):
                    if seq_id in dest:
                        del dest[seq_id]
                    src.copy(seq_id, dest)
        except Exception as e:
            self.logger.error(f"Merging failed: {str(e)}")
            raise

    def process_all(self, fasta_dir: str, gff_dir: str, output_file: str):
        """Robust processing with empty file check"""
        file_pairs = get_file_pairs(fasta_dir, gff_dir)
        Path(output_file).unlink(missing_ok=True)
        
        batches = [file_pairs[i:i+self.batch_size] 
                   for i in range(0, len(file_pairs), self.batch_size)]
        
        for batch in tqdm(batches, desc="Processing batches"):
            self.process_batch(batch, output_file)
            
        # Final validation
        with h5py.File(output_file, 'r') as f:
            if not f.keys():
                raise ValueError("No sequences processed - check input files and error logs")