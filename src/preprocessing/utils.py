import logging
from pathlib import Path
import os

def setup_directories(base_dir):
    """Centralized directory creation"""
    dirs = [
        'data/fasta',
        'data/gff',
        'data/processed',
        'logs',
        'tests/test_data'
    ]
    for dir_path in dirs:
        Path(os.path.join(base_dir, dir_path)).mkdir(parents=True, exist_ok=True)

def setup_logger(name='papa_smorf'):
    """Centralized logger creation"""
    if logging.getLogger(name).handlers:
        return logging.getLogger(name)
        
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    
    # Create logs directory if it doesn't exist
    Path('logs').mkdir(exist_ok=True)
    
    # File handler
    fh = logging.FileHandler(f'logs/{name}.log')
    fh.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    ))
    logger.addHandler(fh)
    
    # Console handler
    ch = logging.StreamHandler()
    ch.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(ch)
    
    return logger

def get_logger(name='papa_smorf'):
    """Alias for setup_logger for compatibility"""
    return setup_logger(name)

def get_file_pairs(fasta_dir, gff_dir):
    """Get matching FASTA and GFF file pairs"""
    pairs = []
    for fasta_file in Path(fasta_dir).glob('*.fna*'):
        gff_file = Path(gff_dir) / f"{fasta_file.stem}.gff"
        if gff_file.exists():
            pairs.append({
                'fasta': str(fasta_file),
                'gff': str(gff_file)
            })
    return pairs

def load_config():
    """Placeholder for config loading if needed"""
    return {}
