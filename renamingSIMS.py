import os
from pathlib import Path
import shutil

def rename_sample_folders(base_path='CAMISIM/out'):
    # Create raw_reads directory if it doesn't exist
    raw_reads_dir = os.path.join(os.path.dirname(base_path), 'raw_reads')
    os.makedirs(raw_reads_dir, exist_ok=True)
    
    # Get all directories in the base path
    dirs = [d for d in os.listdir(base_path) if os.path.isdir(os.path.join(base_path, d)) and 'sample' in d]
    
    # Sort directories to ensure consistent ordering
    dirs.sort()
    
    # Rename folders and their contents
    for index, dir_name in enumerate(dirs, start=1):
        old_path = os.path.join(base_path, dir_name)
        new_name = f'sample_{index}'
        new_path = os.path.join(base_path, new_name)
        
        # First rename and move the reads file
        reads_path = os.path.join(old_path, 'reads', 'anonymous_reads.fq.gz')
        if os.path.exists(reads_path):
            new_reads_name = f'sample_{index}.fastq.gz'
            temp_reads_path = os.path.join(old_path, 'reads', new_reads_name)
            final_reads_path = os.path.join(raw_reads_dir, new_reads_name)
            try:
                # First rename
                os.rename(reads_path, temp_reads_path)
                # Then move to raw_reads
                shutil.move(temp_reads_path, final_reads_path)
                print(f'Moved and renamed reads: anonymous_reads.fq.gz → {final_reads_path}')
            except OSError as e:
                print(f'Error processing reads in {dir_name}: {e}')
        
        # Then rename the sample folder
        try:
            os.rename(old_path, new_path)
            print(f'Renamed folder: {dir_name} → {new_name}')
        except OSError as e:
            print(f'Error renaming folder {dir_name}: {e}')

if __name__ == '__main__':
    rename_sample_folders()