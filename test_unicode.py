# %%
from __future__ import annotations
from pathlib import Path
import argparse
# %%
parser = argparse.ArgumentParser(description="Check which file and line has invalid unicode characters.")
parser.add_argument("path", type=Path, help="Path to the file to check.")
args = parser.parse_args()
# %%
paths = list(args.path.glob('*.f90')) + list(args.path.glob('*.f'))
if not paths:
    print(f"No files found in {args.path} with the specified extensions.")
    exit(1)
# %%
for path in paths:
    with open(path, 'r', encoding='utf-8') as file:
        for line_number, line in enumerate(file, start=1):
            try:
                line.encode('ascii')
            except UnicodeEncodeError:
                print(f"Invalid unicode character found in {path} at line {line_number}:")
                print(line.strip())