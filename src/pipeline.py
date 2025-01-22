# src/pipeline.py
import argparse
import os
import sys
from pathlib import Path
import logging
from datetime import datetime

from src.preprocessing.parallel_processor import ParallelSequenceProcessor
from src.preprocessing.utils import setup_logger  # Changed from get_logger
from src.tools.sequence_checker import check_sequences
from src.tools.inspect_data import inspect_processed_data

class Pipeline:
    def __init__(self, fasta_dir: str, gff_dir: str, output_dir: str, n_workers: int = None):
        self.fasta_dir = fasta_dir
        self.gff_dir = gff_dir
        self.output_dir = output_dir
        self.n_workers = n_workers
        self.logger = setup_logger('pipeline')  # Changed from get_logger
        
        # Create output directories
        os.makedirs(output_dir, exist_ok=True)
        os.makedirs('logs', exist_ok=True)

    def run(self):
        """Run the complete pipeline."""
        start_time = datetime.now()
        self.logger.info("Starting pipeline...")

        try:
            # 1. Check input sequences
            self.logger.info("Step 1: Checking input sequences...")
            check_sequences(self.fasta_dir)

            # 2. Run parallel preprocessing
            self.logger.info("Step 2: Running parallel preprocessing...")
            processor = ParallelSequenceProcessor(n_workers=self.n_workers)
            processor.process_all(self.fasta_dir, self.gff_dir, self.output_dir)

            # 3. Inspect processed data
            self.logger.info("Step 3: Inspecting processed data...")
            inspect_processed_data(self.output_dir)

            # Calculate and log execution time
            end_time = datetime.now()
            execution_time = end_time - start_time
            self.logger.info(f"Pipeline completed successfully in {execution_time}")

        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            raise

def main():
    parser = argparse.ArgumentParser(description="Run papa-smORF pipeline")
    parser.add_argument('--fasta_dir', required=True, help="Directory containing FASTA files")
    parser.add_argument('--gff_dir', required=True, help="Directory containing GFF files")
    parser.add_argument('--output_dir', required=True, help="Output directory")
    parser.add_argument('--workers', type=int, default=None, help="Number of worker processes")
    
    args = parser.parse_args()
    
    pipeline = Pipeline(
        fasta_dir=args.fasta_dir,
        gff_dir=args.gff_dir,
        output_dir=args.output_dir,
        n_workers=args.workers
    )
    
    pipeline.run()

if __name__ == "__main__":
    main()
