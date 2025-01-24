from concurrent.futures import ProcessPoolExecutor
import multiprocessing
from pathlib import Path
from datetime import datetime
import pickle
from tqdm import tqdm
import logging
import json
import numpy as np
import mmap  # Added this import
from src.preprocessing.sequence_processor import SequenceProcessor
from src.preprocessing.utils import setup_logger, get_file_pairs, MMapHandler

class ParallelSequenceProcessor:
    def __init__(self, n_workers=None, batch_size=100):
        self.n_workers = n_workers or multiprocessing.cpu_count()
        self.batch_size = batch_size
        self.processor = SequenceProcessor()
        self.logger = setup_logger('parallel_processor')
        self.mmap_handler = MMapHandler()

    def save_result_mmap(self, result, output_dir):
        """Save processed data using memory mapping"""
        chunk_id = result['chunk_id']
        output_path = Path(output_dir) / f"{chunk_id}.mmap"
        
        try:
            # Calculate sizes and offsets
            sequence = result['sequence']
            seq_size = len(sequence)
            orf_size = len(result['orf_positions']) * 8
            start_codon_size = len(result['start_codon_positions']) * 4
            
            total_size = self.mmap_handler.header_size + seq_size + orf_size + start_codon_size
            
            metadata = {
                'chunk_id': chunk_id,
                'sample_id': result['sample_id'],
                'sequence_length': seq_size,
                'n_orfs': len(result['orf_positions']),
                'n_start_codons': len(result['start_codon_positions']),
                'offsets': {
                    'sequence': self.mmap_handler.header_size,
                    'orf_positions': self.mmap_handler.header_size + seq_size,
                    'start_codons': self.mmap_handler.header_size + seq_size + orf_size
                }
            }
            
            # Create and write to memory-mapped file
            with self.mmap_handler.create_mmap_file(output_path, total_size) as f:
                mm = mmap.mmap(f.fileno(), 0)
                
                # Write data
                self.mmap_handler.write_metadata(mm, metadata)
                self.mmap_handler.write_sequence_data(mm, metadata['offsets']['sequence'], sequence)
                self.mmap_handler.write_numeric_data(mm, metadata['offsets']['orf_positions'], result['orf_positions'])
                self.mmap_handler.write_numeric_data(mm, metadata['offsets']['start_codons'], result['start_codon_positions'])
                
                mm.flush()
                mm.close()
                
            return metadata
            
        except Exception as e:
            self.logger.error(f"Error saving memory-mapped data for {chunk_id}: {str(e)}")
            raise

    def _save_progress(self, results, output_dir, progress_info):
        """Save detailed progress information"""
        try:
            progress_data = {
                'timestamp': str(datetime.now()),
                'processed_files': progress_info['processed_files'],
                'total_files': progress_info['total_files'],
                'completion_percentage': progress_info['completion_percentage'],
                'current_batch': progress_info['current_batch'],
                'total_batches': progress_info['total_batches'],
                'total_orfs': sum(r.get('n_orfs', 0) for r in results),
                'memory_usage': psutil.Process().memory_info().rss / 1024 / 1024 / 1024,  # GB
                'cpu_usage': psutil.Process().cpu_percent()
            }
            
            progress_file = Path(output_dir) / 'progress.json'
            with open(progress_file, 'w') as f:
                json.dump(progress_data, f, indent=2)
                
        except Exception as e:
            self.logger.error(f"Error saving progress: {str(e)}")

    def process_batch(self, file_pairs, output_dir):
        """Process a batch of files"""
        batch_results = []
        with ProcessPoolExecutor(max_workers=self.n_workers) as executor:
            futures = [executor.submit(self.processor.parse_fasta_mmap, pair['fasta']) 
                      for pair in file_pairs]
            
            for future, pair in zip(futures, file_pairs):
                try:
                    sequences = future.result()
                    annotations = self.processor.parse_gff_mmap(pair['gff'])
                    
                    results = self.processor.process_sequences(sequences, annotations)
                    for result in results:
                        metadata = self.save_result_mmap(result, output_dir)
                        batch_results.append(metadata)
                        
                except Exception as e:
                    self.logger.error(f"Error processing {pair}: {str(e)}")
        
        return batch_results

    def process_all(self, fasta_dir, gff_dir, output_dir):
        """Process all files with detailed progress tracking"""
        file_pairs = get_file_pairs(fasta_dir, gff_dir)
        total_pairs = len(file_pairs)
        
        self.logger.info(f"Found {total_pairs} FASTA-GFF pairs to process")
        
        processed_count = 0
        all_results = []
        
        for i in range(0, total_pairs, self.batch_size):
            batch = file_pairs[i:i + self.batch_size]
            current_batch = i//self.batch_size + 1
            total_batches = (total_pairs - 1)//self.batch_size + 1
            
            self.logger.info(f"\n=== Processing Batch {current_batch}/{total_batches} "
                           f"({(current_batch/total_batches)*100:.1f}% complete) ===")
            
            batch_results = self.process_batch(batch, output_dir)
            all_results.extend(batch_results)
            
            processed_count += len(batch)
            
            # Save and display progress
            self._save_progress(all_results, output_dir, {
                'total_files': total_pairs,
                'processed_files': processed_count,
                'completion_percentage': (processed_count/total_pairs)*100,
                'current_batch': current_batch,
                'total_batches': total_batches
            })
            
            self.logger.info(f"Progress: {processed_count}/{total_pairs} files "
                           f"({(processed_count/total_pairs)*100:.1f}% complete)")
            self.logger.info(f"Total ORFs found so far: {sum(r.get('n_orfs', 0) for r in all_results):,}")
        
        return all_results

