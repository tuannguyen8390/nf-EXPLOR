import os
import sys
from json import loads

def get_snf_version(filename: str) -> str:
    try:
        with open(filename, 'rb') as f:
            header = loads(f.readline())
        return f'{header["config"]["version"]} {header["config"]["build"]}'
    except Exception as e:
        return f"Error processing {filename}: {str(e)}"

def process_directory_recursively(directory_path: str):
    if not os.path.isdir(directory_path):
        print(f"Error: {directory_path} is not a valid directory")
        return
    
    found_files = False
    
    for root, dirs, files in os.walk(directory_path):
        for filename in files:
            if filename.endswith('.snf'):
                found_files = True
                full_path = os.path.join(root, filename)
                relative_path = os.path.relpath(full_path, directory_path)
                version = get_snf_version(full_path)
                print(f"{relative_path}: {version}")
    
    if not found_files:
        print(f"No .snf files found in {directory_path} or its subdirectories")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        directory = sys.argv[1]
    else:
        directory = "."  # Current directory if none specified
    
    process_directory_recursively(directory)
