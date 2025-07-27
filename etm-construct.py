#!/usr/bin/env python3
"""
Generate MOL2 files for visualization of molecular properties on van der Waals surface.
Each MOL2 file represents one property, with surface points as atoms and delta values as charges.
Follows the FMN MOL2 template format.
"""

import pandas as pd
import numpy as np
import os
import sys
import shutil
import argparse
from pathlib import Path

def load_surface_coordinates(surface_file):
    """Load surface coordinates from text file."""
    try:
        # Read surface coordinates (assuming format: x y z)
        coords = []
        with open(surface_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                            coords.append([x, y, z])
                        except ValueError:
                            continue
        
        if not coords:
            raise ValueError("No valid coordinates found in surface file")
            
        return np.array(coords)
    except Exception as e:
        print(f"Error loading surface coordinates: {e}")
        return None

def create_mol2_file(property_name, surface_coords, delta_values, output_dir):
    """
    Create a single MOL2 file for a property following FMN template format.
    Each surface point becomes an atom with the delta value as charge.
    Format: index H x y z H1 1 XXX129 charge_value
    """
    if len(surface_coords) != len(delta_values):
        print(f"Warning: Coordinate count ({len(surface_coords)}) != delta count ({len(delta_values)}) for {property_name}")
        min_len = min(len(surface_coords), len(delta_values))
        surface_coords = surface_coords[:min_len]
        delta_values = delta_values[:min_len]
    
    # Create output file
    mol2_file = output_dir / f"{property_name}_tuning_map.mol2"
    
    with open(mol2_file, 'w') as f:
        # Write header following FMN template
        f.write("@<TRIPOS>MOLECULE\n")
        f.write(f"{property_name}_tuning_map.mol2\n")
        f.write(f"     {len(surface_coords)} 0 0 0\n")
        f.write("SMALL\n")
        f.write("GASTEIGER\n")
        f.write("@<TRIPOS>ATOM       \n")
        
        # Write atoms (surface points) following FMN format
        # Format: index H x y z H1 1 XXX129 charge_value
        for i, (coord, delta) in enumerate(zip(surface_coords, delta_values), 1):
            x, y, z = coord
            f.write(f"{i}    H  {x:.4f}  {y:.4f}  {z:.4f}  H1  1  XXX129  {delta:.6f}\n")
    
    print(f"Created MOL2 file: {mol2_file}")
    return mol2_file

def cleanup_files(molecule_name, keep_inputs=False):
    """Clean up unnecessary files after MOL2 generation."""
    files_to_remove = [
        "etm.tex",
        f"{molecule_name}_props.csv", 
        "properties copy.csv",
        "properties.csv",
        "etm-extract copy.py",
        "test.py",
        "timer.dat"
    ]
    
    dirs_to_remove = [
        "TCE-Files",
        f"{molecule_name}_psi4_inputs copy 2",
        "__pycache__",
        "mol2_files"  # Old MOL2 directory
    ]
    
    if not keep_inputs:
        dirs_to_remove.extend([
            f"{molecule_name}_psi4_inputs",
            f"{molecule_name}_redox_inputs"  # Also clean up redox inputs
        ])
    
    # Remove files
    for file_path in files_to_remove:
        if os.path.exists(file_path):
            try:
                os.remove(file_path)
                print(f"Removed: {file_path}")
            except Exception as e:
                print(f"Could not remove {file_path}: {e}")
    
    # Remove directories
    for dir_path in dirs_to_remove:
        if os.path.exists(dir_path):
            try:
                shutil.rmtree(dir_path)
                print(f"Removed directory: {dir_path}")
            except Exception as e:
                print(f"Could not remove directory {dir_path}: {e}")

def detect_molecule_name():
    """Auto-detect molecule name from available files."""
    # Look for XYZ files (excluding surface files)
    xyz_files = [f for f in os.listdir('.') if f.endswith('.xyz') and 'surface' not in f.lower()]
    
    if not xyz_files:
        print("Error: No molecule XYZ file found!")
        print("Please ensure you have a molecule XYZ file (e.g., water.xyz, benzene.xyz)")
        sys.exit(1)
    
    if len(xyz_files) > 1:
        print(f"Warning: Multiple XYZ files found: {xyz_files}")
        print(f"Using: {xyz_files[0]}")
    
    molecule_name = xyz_files[0].replace('.xyz', '')
    return molecule_name

def main():
    parser = argparse.ArgumentParser(description="Generate MOL2 files for ETM visualization")
    parser.add_argument('-ki', '--keep-inputs', action='store_true', 
                        help='Keep input files (molecule_psi4_inputs directory)')
    parser.add_argument('--excited-state', '-exs', type=int, default=1,
                        help='Number of excited states to calculate (tdscf_states)')
    args = parser.parse_args()
    
    # Try to get molecule name from config first, then auto-detect
    molecule_name = None
    config_file = 'psi4_params.json'
    if os.path.exists(config_file):
        try:
            import json
            with open(config_file) as f:
                config = json.load(f)
            if 'xyz' in config:
                molecule_name = config['xyz'].replace('.xyz', '')
                print(f"Using molecule from config: {molecule_name}")
        except:
            pass
    
    if not molecule_name:
        molecule_name = detect_molecule_name()
        print(f"Detected molecule: {molecule_name}")
    
    # File paths
    properties_file = f"{molecule_name}_etm_properties.csv"
    surface_file = f"{molecule_name}_vdw_surface.txt"
    output_dir = Path(f"{molecule_name}_etm_results")
    
    # Create output directory
    output_dir.mkdir(exist_ok=True)
    
    # Load data
    print("Loading surface coordinates...")
    surface_coords = load_surface_coordinates(surface_file)
    if surface_coords is None:
        print("Failed to load surface coordinates")
        return
    
    print(f"Loaded {len(surface_coords)} surface coordinates")
    
    print("Loading properties data...")
    try:
        df = pd.read_csv(properties_file)
        print(f"Loaded properties data with {len(df)} rows")
        print(f"Columns: {list(df.columns)}")
        
        # Remove the reference row (first row with ref='Ref')
        df = df[df['ref'] != 'Ref'].reset_index(drop=True)
        print(f"After removing reference row: {len(df)} rows")
        
    except Exception as e:
        print(f"Error loading properties file: {e}")
        return
    
    # Get delta columns only (these contain the values we want to visualize)
    delta_cols = [col for col in df.columns if col.startswith('delta_')]
    
    print(f"Found delta columns: {delta_cols}")
    
    # Generate MOL2 files for each property
    for delta_col in delta_cols:
        # Extract property name from delta column
        property_name = delta_col.replace('delta_', '')
        
        print(f"\nProcessing property: {property_name}")
        
        # Get delta values for all surface points
        delta_values = df[delta_col].values
        
        # Check if we have the right number of values
        if len(delta_values) != len(surface_coords):
            print(f"Warning: {len(delta_values)} delta values but {len(surface_coords)} surface points")
            # Use minimum to avoid index errors
            min_len = min(len(delta_values), len(surface_coords))
            delta_values = delta_values[:min_len]
            coords_to_use = surface_coords[:min_len]
        else:
            coords_to_use = surface_coords
        
        # Create MOL2 file for this property
        create_mol2_file(property_name, coords_to_use, delta_values, output_dir)
    
    print(f"\nMOL2 files saved in: {output_dir}")
    
    # Clean up unnecessary files
    print("\nCleaning up unnecessary files...")
    cleanup_files(molecule_name, keep_inputs=args.keep_inputs)
    
    print("\nETM construction complete!")
    print(f"Generated {len(delta_cols)} MOL2 files in {output_dir}/")
    if args.keep_inputs:
        print("Input files preserved (--keep-inputs flag used)")
    else:
        print("Input files cleaned up")

if __name__ == "__main__":
    main()
