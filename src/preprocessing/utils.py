# src/preprocessing/utils.py
import logging
import numpy as np
from pathlib import Path
from typing import Dict, List, Any, Tuple
import psutil
import os

def setup_logger(name: str = 'papa_smurf') -> logging.Logger:
    """Set up and configure logger."""
    if logging.getLogger(name).handlers:
        return logging.getLogger(name)
    
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    
    Path('logs').mkdir(exist_ok=True)
    
    fh = logging.FileHandler(f'logs/{name}.log')
    fh.setFormatter(logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    ))
    logger.addHandler(fh)
    
    ch = logging.StreamHandler()
    ch.setFormatter(logging.Formatter('%(message)s'))
    logger.addHandler(ch)
    
    return logger

def get_memory_usage() -> Tuple[float, float]:
    """Get current memory usage in MB and percent."""
    process = psutil.Process()
    memory_info = process.memory_info()
    memory_mb = memory_info.rss / (1024 * 1024)
    memory_percent = process.memory_percent()
    return memory_mb, memory_percent

def get_file_pairs(fasta_dir: str, gff_dir: str) -> List[Dict[str, str]]:
    """Find matching FASTA and GFF file pairs."""
    pairs = []
    fasta_files = list(Path(fasta_dir).glob('*.fna'))
    for fasta_file in fasta_files:
        possible_gff_files = [
            Path(gff_dir) / f"{fasta_file.stem}.gff",
            Path(gff_dir) / f"{fasta_file.stem}.gff3",
            Path(gff_dir) / f"{fasta_file.stem}.gff.gz"
        ]
        for gff_file in possible_gff_files:
            if gff_file.exists():
                pairs.append({
                    'fasta': str(fasta_file),
                    'gff': str(gff_file),
                    'fasta_size': fasta_file.stat().st_size,
                    'gff_size': gff_file.stat().st_size
                })
                break
    return pairs