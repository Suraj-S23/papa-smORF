# src/tools/visualize_data.py
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import SeqIO
import os
from collections import defaultdict
import numpy as np
from src.preprocessing.sequence_processor import SequenceProcessor
from tqdm import tqdm
import concurrent.futures

class SequenceVisualizer:
    def __init__(self, processed_dir, fasta_dir, gff_dir):
        self.processed_dir = processed_dir
        self.fasta_dir = fasta_dir
        self.gff_dir = gff_dir
        self.processor = SequenceProcessor()
        
        # Set standard plotting parameters
        plt.style.use('seaborn-v0_8-paper')  # Clean, professional style
        plt.rcParams['font.family'] = 'DejaVu Serif'
        plt.rcParams['font.size'] = 8
        plt.rcParams['axes.titlesize'] = 10
        plt.rcParams['axes.labelsize'] = 8
        plt.rcParams['xtick.labelsize'] = 7
        plt.rcParams['ytick.labelsize'] = 7
        
        # Create plots directory if it doesn't exist
        os.makedirs('figure_files', exist_ok=True)
        
        # Cache for storing processed data
        self.chunk_df = pd.read_csv(os.path.join(self.processed_dir, 'chunk_index.csv'))
        self.orf_data = defaultdict(list)
        self.start_codons = defaultdict(int)
        
    def preprocess_data(self):
        """Preprocess all data once to avoid multiple reads."""
        print("Preprocessing data...")
        
        # Process chunks in parallel
        with concurrent.futures.ProcessPoolExecutor() as executor:
            futures = []
            for _, row in self.chunk_df.iterrows():
                futures.append(
                    executor.submit(self.processor.load_chunk, 
                                  row['chunk_id'], 
                                  self.processed_dir)
                )
            
            # Collect results
            for future in tqdm(concurrent.futures.as_completed(futures), 
                             total=len(futures), 
                             desc="Processing chunks"):
                chunk_data = future.result()
                sample_id = chunk_data['sample_id']
                chunk_id = chunk_data['chunk_id']
                chunk_num = int(chunk_id.split('_')[-1])
                chunk_start = chunk_num * self.processor.chunk_size
                
                # Store ORF positions
                for start, end in chunk_data['orf_positions']:
                    self.orf_data[sample_id].append(
                        (chunk_start + start, chunk_start + end)
                    )
                
                # Store start codons
                sequence = chunk_data['sequence']
                for pos in chunk_data['start_codon_positions']:
                    if pos + 3 <= len(sequence):
                        codon = sequence[pos:pos+3]
                        self.start_codons[codon] += 1

    def plot_orf_positions(self, sample_id):
        """Plot ORF positions along the genome."""
        orfs = self.orf_data[sample_id]
        if not orfs:
            return
        
        fig, ax = plt.subplots(figsize=(10, 3))
        
        # Plot ORFs more efficiently
        starts, ends = zip(*orfs)
        y = np.ones_like(starts)
        ax.hlines(y=y, xmin=starts, xmax=ends, colors='blue', alpha=0.5, linewidth=0.5)
        
        ax.set_title(f'ORF Positions in {sample_id}')
        ax.set_xlabel('Genome Position (bp)')
        ax.set_ylabel('ORFs')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig(f"figure_files/orf_positions_{sample_id}.pdf", dpi=300, bbox_inches='tight')
        plt.close()

    def plot_orf_length_distribution(self):
        """Plot distribution of ORF lengths."""
        orf_lengths = []
        for orfs in self.orf_data.values():
            orf_lengths.extend(end - start + 1 for start, end in orfs)
        
        fig, ax = plt.subplots(figsize=(6, 4))
        sns.histplot(data=orf_lengths, bins=50, ax=ax)
        
        ax.set_title('Distribution of ORF Lengths')
        ax.set_xlabel('Length (bp)')
        ax.set_ylabel('Count')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        plt.tight_layout()
        plt.savefig("figure_files/orf_length_distribution.pdf", dpi=300, bbox_inches='tight')
        plt.close()

    def plot_start_codons(self):
        """Plot start codon usage."""
        if not self.start_codons:
            return
            
        fig, ax = plt.subplots(figsize=(6, 6))
        
        codons = list(self.start_codons.keys())
        counts = list(self.start_codons.values())
        total = sum(counts)
        percentages = [count/total * 100 for count in counts]
        
        wedges, texts, autotexts = ax.pie(percentages, 
                                         labels=codons,
                                         autopct='%1.1f%%',
                                         textprops={'fontsize': 8})
        
        ax.set_title('Start Codon Usage')
        
        plt.tight_layout()
        plt.savefig("figure_files/start_codon_distribution.pdf", dpi=300, bbox_inches='tight')
        plt.close()

    def generate_all_plots(self):
        """Generate all visualizations."""
        try:
            # Preprocess all data once
            self.preprocess_data()
            
            print("Generating plots...")
            # Generate plots for each sample
            for sample_id in tqdm(self.chunk_df['sample_id'].unique(), desc="Plotting ORF positions"):
                self.plot_orf_positions(sample_id)

            # Generate overall statistics plots
            self.plot_orf_length_distribution()
            self.plot_start_codons()

        except Exception as e:
            print(f"Error generating plots: {str(e)}")
            raise

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Generate visualizations for processed data")
    parser.add_argument('--processed_dir', required=True, help="Directory containing processed data")
    parser.add_argument('--fasta_dir', required=True, help="Directory containing FASTA files")
    parser.add_argument('--gff_dir', required=True, help="Directory containing GFF files")
    
    args = parser.parse_args()
    
    visualizer = SequenceVisualizer(
        processed_dir=args.processed_dir,
        fasta_dir=args.fasta_dir,
        gff_dir=args.gff_dir
    )
    
    visualizer.generate_all_plots()

if __name__ == "__main__":
    main()
