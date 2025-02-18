# src/analysis/start_codon_analysis.py
import h5py
from collections import defaultdict
from pathlib import Path
from Bio.Seq import Seq
import logging

class StartCodonAnalyzer:
    def __init__(self, hdf5_path: str):
        self.logger = logging.getLogger('start_codon_analyzer')
        self.hdf5_path = Path(hdf5_path)
        self.codon_counts = defaultdict(int)
        self.total_codons = 0

    def analyze(self) -> dict:
        self._validate_input()
        self._process_file()
        return self._calculate_percentages()

    def _validate_input(self):
        """Ensure the HDF5 file is non-empty and contains ORF data."""
        with h5py.File(self.hdf5_path, 'r') as f:
            if not f.keys():
                raise ValueError("Empty HDF5 file")
            orf_count = 0
            for seq in f.values():
                for chunk in seq.values():
                    if 'orfs' in chunk:
                        orf_count += len(chunk['orfs'])
            if orf_count == 0:
                raise ValueError("No ORF data found in HDF5 file")

    def _process_file(self):
        """Process the HDF5 file and deduplicate ORFs based on their IDs."""
        processed_ids = set()  # Deduplication set
        with h5py.File(self.hdf5_path, 'r') as f:
            for seq_id in f:
                for chunk_id in f[seq_id]:
                    chunk = f[seq_id][chunk_id]
                    self._process_chunk(chunk, processed_ids)

    def _process_chunk(self, chunk, processed_ids: set):
        """Process a single chunk and count start codons,
           skipping ORFs that have already been processed.
        """
        if 'orfs' not in chunk:
            return
        for orf_id in chunk['orfs']:
            # Do not process the same ORF twice.
            if orf_id in processed_ids:
                continue
            processed_ids.add(orf_id)
            orf = chunk['orfs'][orf_id]
            try:
                strand = orf.attrs['strand']
                seq = orf.attrs['processed_seq']
                if len(seq) < 3:
                    continue

                # For both strands, the ORF's sequence is already in the desired orientation.
                codon = seq[:3]
                self.codon_counts[codon.upper()] += 1
                self.total_codons += 1
            except Exception as e:
                self.logger.warning(f"Error processing {orf_id}: {str(e)}")

    def _calculate_percentages(self) -> dict:
        """Calculate the percentage of each start codon."""
        if self.total_codons == 0:
            raise ValueError("No codons found")
        return {codon: count / self.total_codons for codon, count in self.codon_counts.items()}

    def generate_report(self, percentages: dict) -> str:
        """Generate a human-readable report."""
        report_lines = ["Start Codon Distribution:"]
        for codon, freq in sorted(percentages.items(), key=lambda x: (-x[1], x[0])):
            report_lines.append(f"{codon}: {freq:.2%}")
        return "\n".join(report_lines)

    def validate_reverse_strand_orfs(self, sample_size=100):
        """Return a sample of reverse strand orf end codons for validation."""
        reverse_codons = []
        with h5py.File(self.hdf5_path, 'r') as f:
            for seq in f.values():
                for chunk in seq.values():
                    if 'orfs' not in chunk:
                        continue
                    for orf_id, orf in chunk['orfs'].items():
                        if orf.attrs['strand'] == '-':
                            seq_part = orf.attrs['processed_seq']
                            if len(seq_part) >= 3:
                                # For reverse strand ORFs, we reverse-complement the last 3 bases
                                rev_codon = str(Seq(seq_part[-3:]).reverse_complement()).upper()
                                reverse_codons.append(rev_codon)
                            if len(reverse_codons) >= sample_size:
                                return reverse_codons
        return reverse_codons

# For standalone usage, you can include a simple CLI interface:
if __name__ == "__main__":
    import argparse
    logging.basicConfig(level=logging.INFO)
    
    parser = argparse.ArgumentParser(description="Analyze start codon distribution from processed HDF5 file")
    parser.add_argument('--input', required=True, help='Path to processed HDF5 file')
    args = parser.parse_args()
    
    analyzer = StartCodonAnalyzer(args.input)
    percentages = analyzer.analyze()
    
    print("\nStart Codon Distribution (Percentages):")
    for codon, freq in sorted(percentages.items(), key=lambda x: (-x[1], x[0])):
        print(f"{codon}: {freq:.2%}")
    
    # Optionally, print the deduplicated reverse strand codon sample.
    reverse_samples = analyzer.validate_reverse_strand_orfs()
    print("\nReverse strand sample (deduplicated):", reverse_samples)