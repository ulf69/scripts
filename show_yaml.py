#!/usr/bin/env python3
"""
Convert gorder YAML output to XMGrace/XVG-ready file
Usage:
    python3 gorder_to_grace.py input.dat output.xvg
"""

import yaml
import sys
import os
import re

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <input_file.yaml> <output_file.xvg>")
    sys.exit(1)

input_file = sys.argv[1]
output_file = sys.argv[2]

if not os.path.isfile(input_file):
    print(f"Error: input file '{input_file}' does not exist")
    sys.exit(1)

# Load YAML (skip comment lines)
with open(input_file) as f:
    content = "".join(line for line in f if not line.strip().startswith("#"))
    data = yaml.safe_load(content)

# Write XMGrace/XVG file
with open(output_file, 'w') as out:
    out.write(f"# Carbon\tOrder\n")
    for carbon_key, carbon_data in data['POPC']['order parameters'].items():
        # Remove the number in parentheses
        clean_name = re.sub(r"\s*\(\d+\)$", "", carbon_key)  # C211 (66) -> C211
        clean_name = re.sub(r"^POPC\s+", "", clean_name)        # remove POPC
        out.write(f"{clean_name}\t{carbon_data['total']:.5f}\n")

