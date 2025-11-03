# %%
from __future__ import annotations
from pathlib import Path
import argparse
# %%
parser = argparse.ArgumentParser(description="Check which file and line has invalid unicode characters.")
parser.add_argument("path", type=Path, help="Path to the file to check.")
args = parser.parse_args()
# %%
paths = list(args.path.glob('*.f90'))
if not paths:
    print(f"No files found in {args.path} with the specified extensions.")
    exit(1)
# %%
for path in paths:
    lines = open(path, 'r', encoding='utf-8').readlines()
    for line in lines:
        parts = line.split('!', maxsplit=1)
        if len(parts) > 1:
            parts[0] = parts[0].lower()
        line = '!'.join(parts)
    with open(path, 'w', encoding='utf-8') as file:
        file.writelines(lines)
        