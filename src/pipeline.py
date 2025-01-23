# Change from relative to absolute imports
from src.preprocessing.parallel_processor import ParallelSequenceProcessor
from src.tools.inspect_data import DataInspector
from src.preprocessing.utils import setup_logger
import argparse
from pathlib import Path
import logging
from datetime import datetime

class Pipeline:
    def __init__(self, fasta_dir, gff_dir, output_dir, n_workers=None, batch_size=100):
        self.fasta_dir = fasta_dir
        self.gff_dir = gff_dir
        self.output_dir = output_dir
        self.n_workers = n_workers
        self.batch_size = batch_size
        self.logger = setup_logger('pipeline')
        
        # Create output directories
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        Path('logs').mkdir(exist_ok=True)

    def run(self):
        """Run the complete pipeline with memory mapping"""
        start_time = datetime.now()
        self.logger.info("Starting pipeline...")

        try:
            # Initialize processor
            processor = ParallelSequenceProcessor(
                n_workers=self.n_workers,
                batch_size=self.batch_size
            )
            
            # Process sequences
            self.logger.info("Processing sequences...")
            results = processor.process_all(
                self.fasta_dir,
                self.gff_dir,
                self.output_dir
            )
            
            # Inspect results
            inspector = DataInspector()
            stats = inspector.inspect_processed_data(self.output_dir)

            # Calculate and log execution time
            end_time = datetime.now()
            execution_time = end_time - start_time
            self.logger.info(f"Pipeline completed successfully in {execution_time}")
            
            return results, stats
            
        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            raise

def main():
    parser = argparse.ArgumentParser(description="Run papa-smORF pipeline")
    parser.add_argument('--fasta_dir', required=True)
    parser.add_argument('--gff_dir', required=True)
    parser.add_argument('--output_dir', required=True)
    parser.add_argument('--workers', type=int, default=None)
    parser.add_argument('--batch_size', type=int, default=100)
    
    args = parser.parse_args()
    
    pipeline = Pipeline(
        fasta_dir=args.fasta_dir,
        gff_dir=args.gff_dir,
        output_dir=args.output_dir,
        n_workers=args.workers,
        batch_size=args.batch_size
    )
    
    pipeline.run()

if __name__ == "__main__":
    main()
