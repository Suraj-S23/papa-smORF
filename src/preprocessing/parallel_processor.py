# src/preprocessing/parallel_processor.py
from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from pathlib import Path
import pickle
from tqdm import tqdm
import pandas as pd
import logging
from .sequence_processor import SequenceProcessor
from .utils import setup_logger, get_file_pairs  # Added get_file_pairs import

class ParallelSequenceProcessor:
    def __init__(self, n_workers=None):
        """Initialize parallel processor"""
        self.n_workers = n_workers or multiprocessing.cpu_count()
        self.processor = SequenceProcessor()
        self.logger = setup_logger('papa_smorf')

    def process_file_pair(self, file_pair):
        """Process a single FASTA-GFF file pair"""
        try:
            fasta_file = file_pair['fasta']
            gff_file = file_pair['gff']
            
            # Process single file pair
            sequences = self.processor.parse_multiple_fasta([fasta_file])
            annotations = self.processor.parse_multiple_gff([gff_file])
            
            # Process sequences returns a list of results
            return self.processor.process_sequences(sequences, annotations)
        except Exception as e:
            return {'error': str(e), 'files': file_pair}

    def process_all(self, fasta_dir, gff_dir, output_dir):
        """Process all files in parallel"""
        file_pairs = get_file_pairs(fasta_dir, gff_dir)
        
        self.logger.info(f"Found {len(file_pairs)} FASTA-GFF pairs to process")
        if len(file_pairs) == 0:
            self.logger.warning("No matching FASTA-GFF pairs found")
            return []

        all_results = []
        with ProcessPoolExecutor(max_workers=self.n_workers) as executor:
            futures = [executor.submit(self.process_file_pair, pair) 
                      for pair in file_pairs]
            
            for future in tqdm(futures, total=len(file_pairs), 
                             desc="Processing sequences"):
                result = future.result()
                if isinstance(result, list):  # Successful processing
                    for item in result:
                        all_results.append(item)
                        self._save_result(item, output_dir)
                else:  # Error occurred
                    self.logger.error(f"Error processing {result['files']}: {result['error']}")

        # Create chunk index
        if all_results:
            chunk_index = pd.DataFrame([{
                'chunk_id': r['chunk_id'],
                'sample_id': r['sample_id'],
                'n_orfs': len(r['orf_positions']),
                'sequence_length': len(r['sequence'])
            } for r in all_results])
            
            index_path = Path(output_dir) / 'chunk_index.csv'
            chunk_index.to_csv(index_path, index=False)
            self.logger.info(f"Created chunk index with {len(chunk_index)} entries at {index_path}")

        return all_results


    def _save_result(self, result, output_dir):
        """Save processed data"""
        chunk_id = result['chunk_id']
        output_path = Path(output_dir) / f"{chunk_id}.pkl"
        with open(output_path, 'wb') as f:
            pickle.dump(result, f)

    def _log_errors(self, errors):
        """Log processing errors"""
        for error in errors:
            self.logger.error(f"Error processing {error['files']}: {error['error']}")
