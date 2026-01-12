#!/usr/bin/env python3
"""
Apply drag multiplier: 1x below 10 deg, 2x at/above 10 deg.
Reads from data/v3/2D_polars_CFD/ and writes to data/v3/2D_2x_drag_polars/
"""

import os
import csv

input_dir = "data/v3/2D_polars_CFD"
output_dir = "data/v3/2D_2x_drag_polars"

os.makedirs(output_dir, exist_ok=True)

for filename in os.listdir(input_dir):
    if not filename.endswith('.csv'):
        continue

    filepath = os.path.join(input_dir, filename)

    # Read CSV
    with open(filepath, 'r') as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames

    # Apply drag multiplier: 1x below 10, 2x at/above 10
    for row in rows:
        alpha = float(row['alpha'])
        multiplier = 2.0 if alpha >= 10.0 else 1.0
        row['Cd'] = float(row['Cd']) * multiplier

    # Write output
    outpath = os.path.join(output_dir, filename)
    with open(outpath, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Processed {filename}")

print(f"Done. Output in {output_dir}")
