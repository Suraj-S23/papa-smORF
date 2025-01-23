import pytest
from src.preprocessing.utils import MMapHandler
import numpy as np
import os

@pytest.fixture
def mmap_handler():
    return MMapHandler()

def test_mmap_write_read(mmap_handler, tmp_path):
    # Test data
    test_data = {
        'sequence': 'ATCG',
        'orf_positions': [(0, 3)],
        'start_codon_positions': [0]
    }
    
    # Create test file
    test_file = tmp_path / "test.mmap"
    metadata = {
        'chunk_id': 'test_chunk',
        'sample_id': 'test_sample',
        'sequence_length': len(test_data['sequence']),
        'n_orfs': len(test_data['orf_positions']),
        'n_start_codons': len(test_data['start_codon_positions']),
        'offsets': {
            'sequence': 1024,
            'orf_positions': 1024 + len(test_data['sequence']),
            'start_codons': 1024 + len(test_data['sequence']) + len(test_data['orf_positions']) * 8
        }
    }
    
    # Write test data
    with mmap_handler.create_mmap_file(test_file, 2048) as f:
        mm = mmap.mmap(f.fileno(), 0)
        mmap_handler.write_metadata(mm, metadata)
        mmap_handler.write_sequence_data(mm, metadata['offsets']['sequence'], test_data['sequence'])
        mmap_handler.write_numeric_data(mm, metadata['offsets']['orf_positions'], test_data['orf_positions'])
        mm.close()
    
    # Read and verify
    with open(test_file, 'rb') as f:
        mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        read_metadata = mmap_handler.read_metadata(mm)
        assert read_metadata['chunk_id'] == metadata['chunk_id']
        
        read_sequence = mmap_handler.read_sequence_data(
            mm, 
            metadata['offsets']['sequence'], 
            metadata['sequence_length']
        )
        assert read_sequence == test_data['sequence']
        
        mm.close()
