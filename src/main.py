import argparse
from preprocessing.sequence_processor import SequenceProcessor
from preprocessing.utils import setup_logger

logger = setup_logger()

def main():
    parser = argparse.ArgumentParser(description="Process DNA sequences for smORF detection")
    parser.add_argument('--fasta_dir', required=True, help="Directory containing FASTA files")
    parser.add_argument('--gff_dir', required=True, help="Directory containing GFF files")
    parser.add_argument('--output_dir', required=True, help="Directory to store processed data")
    parser.add_argument('--config', default='configs/config.yaml', help="Path to config file")
    
    args = parser.parse_args()
    
    try:
        # Initialize processor
        processor = SequenceProcessor(config_path=args.config)
        
        # Parse input files
        logger.info("Parsing FASTA files...")
        sequences = processor.parse_multiple_fasta(args.fasta_dir)
        logger.info(f"Total sequences found: {len(sequences)}")
        
        logger.info("Parsing GFF files...")
        annotations = processor.parse_multiple_gff(args.gff_dir)
        logger.info(f"Total annotations found: {len(annotations)}")
        
        # Process sequences
        logger.info("Processing sequences...")
        chunk_index = processor.process_sequences(sequences, annotations, args.output_dir)
        logger.info(f"Processing completed. Created {len(chunk_index)} chunks")
        
    except Exception as e:
        logger.error(f"Error during processing: {str(e)}")
        raise

if __name__ == '__main__':
    main()
