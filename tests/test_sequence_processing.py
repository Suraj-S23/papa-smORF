# tests/test_sequence_processor.py
import unittest
from src.preprocessing.sequence_processor import SequenceProcessor

class TestSequenceProcessor(unittest.TestCase):
    def setUp(self):
        """This runs before each test."""
        self.processor = SequenceProcessor()
        
    def test_chunk_sequence(self):
        """Test if sequence chunking works correctly."""
        test_sequence = "ATGCATGCATGC"  # 12 bases
        chunks = self.processor.chunk_sequence(test_sequence)
        
        # Test if chunks are correct size
        self.assertEqual(len(chunks[0]), self.processor.chunk_size)
        
        # Test if joining chunks gives original sequence
        self.assertEqual(''.join(chunks), test_sequence)

    def test_create_target_vector(self):
        """Test if target vector creation works correctly."""
        sequence_length = 10
        orf_positions = [(1, 5)]  # ORF from position 1 to 5
        
        target = self.processor.create_target_vector(sequence_length, orf_positions)
        
        # First 5 positions should be 1, rest should be 0
        expected = [1, 1, 1, 1, 1, 0, 0, 0, 0, 0]
        self.assertEqual(target, expected)
