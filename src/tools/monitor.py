import psutil
import time
from pathlib import Path
import json

def monitor_resources(output_dir):
    """Monitor system resources during processing"""
    while True:
        try:
            # Get resource usage
            cpu_percent = psutil.cpu_percent(interval=1)
            memory = psutil.virtual_memory()
            disk = psutil.disk_usage('/')
            
            # Get processing progress
            mmap_files = list(Path(output_dir).glob('*.mmap'))
            
            print(f"\nResource Usage:")
            print(f"CPU: {cpu_percent}%")
            print(f"Memory: {memory.percent}%")
            print(f"Disk: {disk.percent}%")
            print(f"Processed files: {len(mmap_files)}")
            
            time.sleep(60)
            
        except KeyboardInterrupt:
            break
