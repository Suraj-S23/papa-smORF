import logging
import mmap
import json
import numpy as np
from pathlib import Path
import os
import psutil
from typing import Dict, List, Any, Union

class MMapHandler:
    """Handle memory-mapped file operations"""
    def __init__(self):
        self.header_size = 1024
        self.logger = setup_logger('mmap_handler')

    def create_mmap_file(self, file_path: Path, size: int) -> Any:
        """Create a new memory-mapped file"""
        with open(file_path, 'wb') as f:
            f.write(b'\0' * size)
        return open(file_path, 'r+b')

    def write_metadata(self, mm: mmap.mmap, metadata: Dict) -> None:
        """Write metadata to memory-mapped file"""
        metadata_bytes = json.dumps(metadata).encode()
        mm.write(metadata_bytes.ljust(self.header_size))

    def read_metadata(self, mm: mmap.mmap) -> Dict:
        """Read metadata from memory-mapped file"""
        mm.seek(0)
        header = mm.read(self.header_size).decode().strip('\0')
        return json.loads(header)

    def write_sequence_data(self, mm: mmap.mmap, offset: int, sequence: str) -> None:
        """Write sequence data to memory-mapped file"""
        mm.seek(offset)
        mm.write(sequence.encode())

    def write_numeric_data(self, mm: mmap.mmap, offset: int, data: Union[List, np.ndarray]) -> None:
        """Write numeric data to memory-mapped file"""
        mm.seek(offset)
        mm.write(np.array(data, dtype=np.int32).tobytes())

    def read_sequence_data(self, mm: mmap.mmap, offset: int, length: int) -> str:
        """Read sequence data from memory-mapped file"""
        mm.seek(offset)
        return mm.read(length).decode()

    def read_numeric_data(self, mm: mmap.mmap, offset: int, count: int, shape: tuple = None) -> np.ndarray:
        """Read numeric data from memory-mapped file"""
        mm.seek(offset)
        data = np.frombuffer(mm.read(count * 4), dtype=np.int32)
        return data.reshape(shape) if shape else data

def setup_logger(name: str = 'papa_smorf') -> logging.Logger:
    """Set up and configure logger"""
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

def get_file_pairs(fasta_dir: str, gff_dir: str) -> List[Dict[str, str]]:
    """Get matching FASTA and GFF file pairs"""
    pairs = []
    # Focus on .fna files for bacterial genomes
    fasta_files = list(Path(fasta_dir).glob('*.fna'))
    
    for fasta_file in fasta_files:
        # Try different GFF extensions
        possible_gff_files = [
            Path(gff_dir) / f"{fasta_file.stem}.gff",
            Path(gff_dir) / f"{fasta_file.stem}.gff3",
            Path(gff_dir) / f"{fasta_file.stem}.gff.gz"
        ]
        
        for gff_file in possible_gff_files:
            if gff_file.exists():
                pairs.append({
                    'fasta': str(fasta_file),
                    'gff': str(gff_file)
                })
                break
    
    return pairs

def validate_mmap_file(file_path: Union[str, Path]) -> bool:
    """Validate memory-mapped file structure"""
    try:
        with open(file_path, 'rb') as f:
            mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
            metadata = json.loads(mm.read(1024).decode().strip('\0'))
            
            # Check required fields
            required_fields = ['chunk_id', 'sample_id', 'sequence_length', 'n_orfs', 'offsets']
            for field in required_fields:
                if field not in metadata:
                    raise ValueError(f"Missing required field: {field}")
            
            # Validate file size
            expected_size = max(metadata['offsets'].values()) + \
                          metadata['n_orfs'] * 8 + \
                          metadata['n_start_codons'] * 4
            actual_size = mm.size()
            
            if actual_size < expected_size:
                raise ValueError(f"File size mismatch: expected {expected_size}, got {actual_size}")
            
            mm.close()
            return True
    except Exception as e:
        logging.error(f"Invalid mmap file {file_path}: {str(e)}")
        return False

def check_resources() -> Dict[str, float]:
    """Check system resources"""
    cpu_percent = psutil.cpu_percent(interval=1)
    memory = psutil.virtual_memory()
    disk = psutil.disk_usage('/')
    
    return {
        'cpu_percent': cpu_percent,
        'memory_percent': memory.percent,
        'disk_percent': disk.percent,
        'memory_available_gb': memory.available / (1024**3),
        'disk_available_gb': disk.free / (1024**3)
    }

def setup_directories(base_dir: Union[str, Path]) -> None:
    """Create project directories"""
    dirs = [
        'data/fasta',
        'data/gff',
        'data/processed',
        'logs',
        'tests/test_data'
    ]
    for dir_path in dirs:
        Path(os.path.join(base_dir, dir_path)).mkdir(parents=True, exist_ok=True)

def get_sequence_region(mmap_file: Union[str, Path], start: int, end: int) -> str:
    """Get specific region from sequence in memory-mapped file"""
    handler = MMapHandler()
    with open(mmap_file, 'rb') as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        metadata = handler.read_metadata(mm)
        sequence = handler.read_sequence_data(
            mm,
            metadata['offsets']['sequence'] + start,
            end - start
        )
        mm.close()
        return sequence
