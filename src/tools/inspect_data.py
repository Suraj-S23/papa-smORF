import pandas as pd
from pathlib import Path
from src.preprocessing.utils import MMapHandler
import mmap
import logging


class DataInspector:
    def __init__(self):
        self.mmap_handler = MMapHandler()

    def inspect_processed_data(self, processed_dir):
        """Inspect processed data using memory mapping"""
        processed_dir = Path(processed_dir)
        mmap_files = list(processed_dir.glob('*.mmap'))
        
        if not mmap_files:
            print(f"No processed files found in {processed_dir}")
            return
        
        total_orfs = 0
        total_sequences = 0
        sequence_stats = []
        
        for mmap_file in mmap_files:
            with open(mmap_file, 'rb') as f:
                mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
                metadata = self.mmap_handler.read_metadata(mm)
                
                sequence_stats.append({
                    'sample_id': metadata['sample_id'],
                    'sequence_length': metadata['sequence_length'],
                    'n_orfs': metadata['n_orfs'],
                    'n_start_codons': metadata['n_start_codons']
                })
                
                total_orfs += metadata['n_orfs']
                total_sequences += 1
                mm.close()
        
        stats_df = pd.DataFrame(sequence_stats)
        
        print("\nProcessed Data Summary:")
        print(f"Total sequences: {total_sequences}")
        print(f"Total ORFs: {total_orfs:,}")
        print(f"\nSequence Statistics:")
        print(f"Average sequence length: {stats_df['sequence_length'].mean():,.0f}")
        print(f"Average ORFs per sequence: {stats_df['n_orfs'].mean():,.1f}")
        
        return stats_df

    def get_sequence_region(self, mmap_file, start, end):
        """Get specific region from sequence"""
        with open(mmap_file, 'rb') as f:
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            metadata = self.mmap_handler.read_metadata(mm)
            
            sequence = self.mmap_handler.read_sequence_data(
                mm,
                metadata['offsets']['sequence'] + start,
                end - start
            )
            mm.close()
            return sequence
