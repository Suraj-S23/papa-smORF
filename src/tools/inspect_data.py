import pandas as pd
import os
import sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
from src.preprocessing.sequence_processor import SequenceProcessor

def inspect_processed_data(processed_dir):
    """Inspect processed data with better error handling"""
    try:
        index_file = os.path.join(processed_dir, 'chunk_index.csv')
        if not os.path.exists(index_file):
            print(f"\nNo chunk index found at {index_file}")
            # Try to find any processed files
            pkl_files = list(Path(processed_dir).glob('*.pkl'))
            if pkl_files:
                print(f"Found {len(pkl_files)} processed files but no index.")
                return
            else:
                print("No processed files found.")
                return
            
        chunk_df = pd.read_csv(index_file)
        print("\nProcessed Data Summary:")
        print(f"Total chunks: {len(chunk_df)}")
        print(f"Total samples: {chunk_df['sample_id'].nunique()}")
        print(f"Total ORFs: {chunk_df['n_orfs'].sum():,}")
        print("\nPer-sample statistics:")
        print(chunk_df.groupby('sample_id')['n_orfs'].sum().to_string())
        
    except Exception as e:
        print(f"Error inspecting data: {str(e)}")

if __name__ == "__main__":
    inspect_processed_data("data/processed")
