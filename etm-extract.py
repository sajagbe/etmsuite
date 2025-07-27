# --- Property Calculation Functions (fix for NameError) ---
def electron_affinity(eneutral, eanion):
    return eneutral - eanion

def ionization_energy(ecation, eneutral):
    return ecation - eneutral

def chemical_potential(I, A):
    return -0.5 * (I + A)

def electronegativity(I, A):
    return 0.5 * (I + A)

def hardness(I, A):
    return 0.5 * (I - A)

def electrophilicity(mu, eta):
    return mu**2 / (2 * eta) if eta != 0 else None

def nucleophilicity(mu, eta):
    omega = electrophilicity(mu, eta)
    return 2 * eta / mu**2 if mu != 0 else None
import os
import re
import sys
import json
import argparse
import pandas as pd
import numpy as np
import glob
from collections import defaultdict

try:
    from cclib.io import ccread
    CCLIB_AVAILABLE = True
except ImportError:
    CCLIB_AVAILABLE = False
    print("Warning: cclib not available, using regex parsing only")

# --- Property Calculation Functions ---
def ground_state_energy(f):
    if CCLIB_AVAILABLE:
        try:
            return ccread(f).scfenergies[-1]
        except:
            pass
    
    # Fallback regex parsing
    pattern = re.compile(r'@DF-RKS Final Energy:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)')
    try:
        with open(f, 'r') as infile:
            for line in infile:
                m = pattern.search(line)
                if m:
                    return float(m.group(1))
    except Exception as e:
        print(f"Error reading ground state energy from {f}: {e}")
    return None

def homo(f):
    if CCLIB_AVAILABLE:
        try:
            d = ccread(f)
            return d.moenergies[0][d.homos[0]]
        except:
            pass
    
    # Fallback regex parsing
    au_to_ev = 27.2114
    pattern = re.compile(r'(\d+)A\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)')
    try:
        with open(f, 'r') as infile:
            energies = []
            in_orbital_section = False
            for line in infile:
                if "Orbital Energies [Eh]" in line:
                    in_orbital_section = True
                    continue
                if in_orbital_section and "Doubly Occupied:" in line:
                    continue
                if in_orbital_section and "Virtual:" in line:
                    break
                if in_orbital_section:
                    m = pattern.search(line)
                    if m:
                        energies.append(float(m.group(2)) * au_to_ev)
            if energies:
                return energies[-1]  # HOMO is the last occupied orbital
    except Exception as e:
        print(f"Error reading HOMO from {f}: {e}")
    return None

def lumo(f):
    if CCLIB_AVAILABLE:
        try:
            d = ccread(f)
            return d.moenergies[0][d.homos[0] + 1]
        except:
            pass
    
    # Fallback regex parsing
    au_to_ev = 27.2114
    pattern = re.compile(r'(\d+)A\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)')
    try:
        with open(f, 'r') as infile:
            in_virtual_section = False
            for line in infile:
                if "Virtual:" in line:
                    in_virtual_section = True
                    continue
                if in_virtual_section:
                    m = pattern.search(line)
                    if m:
                        return float(m.group(2)) * au_to_ev  # First virtual orbital is LUMO
    except Exception as e:
        print(f"Error reading LUMO from {f}: {e}")
    return None

def dipole_moment(f):
    if CCLIB_AVAILABLE:
        try:
            return ccread(f).moments[1][-1]
        except:
            pass
    
    # Fallback regex parsing
    pattern = re.compile(r'Magnitude\s*:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)')
    try:
        with open(f, 'r') as infile:
            for line in infile:
                m = pattern.search(line)
                if m:
                    return float(m.group(1))
    except Exception as e:
        print(f"Error reading dipole moment from {f}: {e}")
    return None

def extract_excitation_energy(f, state_number):
    """Extract excitation energy for the specified excited state using cclib if available, fallback to regex parsing"""
    if CCLIB_AVAILABLE:
        try:
            data = ccread(f)
            if hasattr(data, 'etenergies') and data.etenergies is not None and len(data.etenergies) >= state_number:
                return data.etenergies[state_number-1]
        except:
            pass
    # Fallback to regex parsing
    au_to_ev = 27.2114
    pattern = re.compile(rf"Excited State\s+{state_number}\s+\([^)]+\):\s+([\d\.Ee+-]+)\s+au")
    try:
        with open(f, 'r') as infile:
            for line in infile:
                m = pattern.search(line)
                if m:
                    return float(m.group(1)) * au_to_ev
    except Exception as e:
        print(f"Error reading S{state_number} excitation energy from {f}: {e}")
    return None

def extract_oscillator_strength(f, state_number):
    """Extract oscillator strength for the specified excited state using cclib if available, fallback to regex parsing"""
    if CCLIB_AVAILABLE:
        try:
            data = ccread(f)
            if hasattr(data, 'etoscs') and data.etoscs is not None and len(data.etoscs) >= state_number:
                return data.etoscs[state_number-1]
        except:
            pass
    # Fallback to regex parsing
    pattern = re.compile(rf"Excited State\s+{state_number}\s+\([^)]+\):.*f\s*=\s*([\d\.Ee+-]+)")
    try:
        with open(f, 'r') as infile:
            for line in infile:
                m = pattern.search(line)
                if m:
                    return float(m.group(1))
    except Exception as e:
        print(f"Error reading S{state_number} oscillator strength from {f}: {e}")
    return None

# Generate dynamic excited state functions
def create_excited_state_functions():
    """Dynamically create excited state extraction functions"""
    globals_dict = globals()
    
    for i in range(1, 11):  # s1 through s10
        # Create excitation energy function
        def make_excitation_func(state_num):
            def func(f):
                return extract_excitation_energy(f, state_num)
            func.__name__ = f's{state_num}_excitation_energy'
            func.__doc__ = f'Extract S{state_num} excitation energy'
            return func
        
        # Create oscillator strength function  
        def make_oscillator_func(state_num):
            def func(f):
                return extract_oscillator_strength(f, state_num)
            func.__name__ = f's{state_num}_oscillator_strength'
            func.__doc__ = f'Extract S{state_num} oscillator strength'
            return func
        
        # Add functions to global namespace
        globals_dict[f's{i}_excitation_energy'] = make_excitation_func(i)
        globals_dict[f's{i}_oscillator_strength'] = make_oscillator_func(i)

# Generate the functions
create_excited_state_functions()

def oscillator_strength(f):
    """Extract oscillator strength for the specified excited state"""
    excited_state = getattr(oscillator_strength, 'excited_state', 1)
    if CCLIB_AVAILABLE:
        try:
            data = ccread(f)
            if hasattr(data, 'etoscs') and data.etoscs is not None and len(data.etoscs) >= excited_state:
                return data.etoscs[excited_state-1]
        except:
            pass
    # Fallback to regex parsing
    pattern = re.compile(rf"Excited State\s+{excited_state}\s+\([^)]+\):.*f\s*=\s*([\d\.Ee+-]+)")
    try:
        with open(f, 'r') as infile:
            for line in infile:
                m = pattern.search(line)
                if m:
                    return float(m.group(1))
    except Exception as e:
        print(f"Error reading S{excited_state} oscillator strength from {f}: {e}")
    return None

def excitation_energy(f):
    """Extract excitation energy for the specified excited state"""
    excited_state = getattr(excitation_energy, 'excited_state', 1)
    au_to_ev = 27.2114
    if CCLIB_AVAILABLE:
        try:
            data = ccread(f)
            if hasattr(data, 'etenergies') and data.etenergies is not None and len(data.etenergies) >= excited_state:
                return data.etenergies[excited_state-1]
        except:
            pass
    # Fallback to regex parsing
    pattern = re.compile(rf"Excited State\s+{excited_state}\s+\([^)]+\):\s*([\d\.Ee+-]+)\s*au")
    try:
        with open(f, 'r') as infile:
            for line in infile:
                m = pattern.search(line)
                if m:
                    return float(m.group(1)) * au_to_ev
    except Exception as e:
        print(f"Error reading S{excited_state} excitation energy from {f}: {e}")
    return None

def gibbs_correction(f):
    """Extract Gibbs free energy correction from frequency calculation output"""
    try:
        with open(f, 'r') as infile:
            for line in infile:
                if "Correction G" in line and "[Eh]" in line:
                    # Extract the last number in the line (in Hartree/Eh)
                    parts = line.strip().split()
                    return float(parts[-2])  # Second to last element (before [Eh])
    except Exception as e:
        print(f"Error reading Gibbs correction from {f}: {e}")
    return None

def calculate_redox_properties(neutral_sscf_file, neutral_gfec_file, anion_sscf_file, anion_gfec_file):
    """Calculate redox properties from the 4 calculation files"""
    
    # Extract energies (no conversion - keep in atomic units)
    E_neutral_sscf = ground_state_energy(neutral_sscf_file)  # Solvated SCF energy
    G_corr_neutral = gibbs_correction(neutral_gfec_file)     # Gas-phase Gibbs correction
    E_anion_sscf = ground_state_energy(anion_sscf_file)      # Solvated SCF energy  
    G_corr_anion = gibbs_correction(anion_gfec_file)         # Gas-phase Gibbs correction

    # Always convert extracted SCF energies from eV to au (corrections are already in au)
    E_neutral_sscf = E_neutral_sscf / 27.2114
    E_anion_sscf = E_anion_sscf / 27.2114

    # Check if all extractions succeeded
    if None in [E_neutral_sscf, G_corr_neutral, E_anion_sscf, G_corr_anion]:
        return None, None

    # Calculate total Gibbs free energies (solvated SCF + gas-phase correction)
    gibbs_neutral = E_neutral_sscf + G_corr_neutral
    gibbs_anion = E_anion_sscf + G_corr_anion

    # Calculate reduction free energy (in atomic units)
    reduction_gfe = gibbs_anion - gibbs_neutral

    # Convert to redox potential (eV vs NHE)
    reduction_gfe_eV = reduction_gfe * 27.2114  # Convert au to eV
    redox_potential = (-1 * reduction_gfe_eV) - 4.43  # Apply transformation as requested


    return reduction_gfe, redox_potential  # reduction_gfe in au, redox_potential in eV
def build_file_map_for_index(base_dir, target_index, molecule_name):
    """
    Scan the base directory subfolders for files matching {molecule}_{species}_{charge|Ref}.out
    and collect paths ONLY for the given target_index.
    Returns dict: {'neutral': path or None, 'anion': path or None, 'cation': path or None}
    """
    file_map = {'neutral': None, 'anion': None, 'cation': None}
    
    # Enhanced pattern to handle multiple naming conventions:
    # 1. {molecule}_{species}_{charge_####|Ref}[_attempt#].out (local mode with retries)
    # 2. {molecule}_{species}.out (local mode only)
    # 3. {molecule}_{species}_{charge_####|Ref}.out (standard format)
    pattern = re.compile(
        rf'{molecule_name}_(neutral|anion|cation)(?:_(?:charge_(\d{{4}})|Ref))?(?:_attempt\d+)?\.out$'
    )

    for root, _, files in os.walk(base_dir):
        for file in files:
            if not file.endswith(".out"):
                continue
            match = pattern.match(file)
            if not match:
                continue
            species, index = match.groups()
            
            # Priority system for file selection:
            # 1. If target_index is 'Ref', prefer files with no index (local style) or explicit 'Ref'
            # 2. If target_index is a charge, match exact charge
            # 3. For files without charge suffix, treat as 'Ref' equivalent
            
            if target_index == 'Ref':
                if index is None or index == 'Ref':
                    # Only overwrite if we don't have a file yet, or if this is explicitly 'Ref'
                    if file_map[species] is None:
                        file_map[species] = os.path.join(root, file)
                        # Debug output can be enabled by uncommenting the next line
                        # print(f"[DEBUG] Mapped {species} -> {file} for index {target_index}")
            elif index == target_index:
                file_map[species] = os.path.join(root, file)
                # Debug output can be enabled by uncommenting the next line
                # print(f"[DEBUG] Mapped {species} -> {file} for index {target_index}")
    
    return file_map

def build_redox_file_map_for_index(base_dir, target_index, molecule_name):
    """
    Scan the flat redox directory structure for files matching the target index.
    Expected structure: molecule_redox_inputs/{anion_gfec, anion_sscf, neutral_gfec, neutral_sscf}/molecule_{charge|species}_{index}.out
    Returns dict with paths to the 4 required redox calculation files.
    """
    file_map = {
        'neutral_gfec': None,
        'neutral_sscf': None, 
        'anion_gfec': None,
        'anion_sscf': None
    }


    # Define the redox directory structure (flat)
    # Only join molecule_name_redox_inputs if not already present
    if base_dir.endswith(f"{molecule_name}_redox_inputs"):
        redox_base = base_dir
    else:
        redox_base = os.path.join(base_dir, f"{molecule_name}_redox_inputs")

    # Map the file types to their flat subfolders and naming patterns
    file_types = {
        'neutral_gfec': ('neutral_gfec', f'{molecule_name}_neutral_gfec'),
        'neutral_sscf': ('neutral_sscf', f'{molecule_name}_neutral_sscf'),
        'anion_gfec': ('anion_gfec', f'{molecule_name}_anion_gfec'),
        'anion_sscf': ('anion_sscf', f'{molecule_name}_anion_sscf')
    }


    for file_type, (subdir, name_prefix) in file_types.items():
        target_dir = os.path.join(redox_base, subdir)
        if os.path.exists(target_dir):
            # Look for files matching the target index
            if target_index == 'Ref':
                pattern = f"{name_prefix}_Ref*.out"
            else:
                pattern = f"{name_prefix}_charge_{target_index}*.out"

            full_pattern = os.path.join(target_dir, pattern)
            matching_files = glob.glob(full_pattern)
            if matching_files:
                # Sort files by attempt number (if present), highest first
                def attempt_sort_key(f):
                    import re
                    m = re.search(r'_attempt(\\d+)\\.out$', f)
                    return int(m.group(1)) if m else 0
                matching_files.sort(key=attempt_sort_key, reverse=True)
                file_map[file_type] = matching_files[0]

    return file_map

def compute_properties_for_index(index, file_map, props_requested, base_dir=None, molecule_name=None):
    """
    Given an index and its file paths (file_map), compute the requested properties.
    Returns a dict with results.
    """
    row = {}

    try:
        if 'ground_state_energy' in props_requested and file_map['neutral']:
            eneutral = ground_state_energy(file_map['neutral'])
            row['ground_state_energy'] = eneutral
        else:
            eneutral = None

        if 'homo' in props_requested and file_map['neutral']:
            ehomo = homo(file_map['neutral'])
            row['homo'] = ehomo
        else:
            ehomo = None

        if 'lumo' in props_requested and file_map['neutral']:
            elumo = lumo(file_map['neutral'])
            row['lumo'] = elumo

        if 'gap' in props_requested and 'homo' in row and 'lumo' in row:
            row['gap'] = row['lumo'] - row['homo']

        if 'dipole_moment' in props_requested and file_map['neutral']:
            row['dipole_moment'] = dipole_moment(file_map['neutral'])

        # Handle dynamic excited state properties
        for prop in props_requested:
            if prop.startswith('s') and ('_excitation_energy' in prop or '_oscillator_strength' in prop) and file_map['neutral']:
                if prop in globals():
                    row[prop] = globals()[prop](file_map['neutral'])

        # Generic excited state properties
        if 'excitation_energy' in props_requested and file_map['neutral']:
            row['excitation_energy'] = excitation_energy(file_map['neutral'])

        if 'oscillator_strength' in props_requested and file_map['neutral']:
            row['oscillator_strength'] = oscillator_strength(file_map['neutral'])

        if {'electron_affinity', 'chemical_potential', 'electronegativity', 'hardness', 'electrophilicity'} & props_requested:
            if file_map['neutral'] and file_map['anion']:
                if eneutral is None:
                    eneutral = ground_state_energy(file_map['neutral'])
                eanion = ground_state_energy(file_map['anion'])
                A = electron_affinity(eneutral, eanion)
                row['electron_affinity'] = A

        if {'ionization_energy', 'chemical_potential', 'electronegativity', 'hardness', 'electrophilicity'} & props_requested:
            if file_map['neutral'] and file_map['cation']:
                if eneutral is None:
                    eneutral = ground_state_energy(file_map['neutral'])
                ecation = ground_state_energy(file_map['cation'])
                I = ionization_energy(ecation, eneutral)
                row['ionization_energy'] = I

        if {'chemical_potential', 'electronegativity', 'hardness', 'electrophilicity', 'nucleophilicity'} & props_requested:
            I = row.get('ionization_energy')
            A = row.get('electron_affinity')
            if I is not None and A is not None:
                mu = chemical_potential(I, A)
                eta = hardness(I, A)
                omega = electrophilicity(mu, eta)
                nu = nucleophilicity(mu, eta)
                row.update({
                    'chemical_potential': mu,
                    'electronegativity': electronegativity(I, A),
                    'hardness': eta,
                    'electrophilicity': omega,
                    'nucleophilicity': nu
                })

        # Remove the old nucleophilicity calculation based on HOMO

        # Handle redox properties
        if 'redox' in props_requested:
            # Robustly detect the redox input directory for any invocation context
            redox_dir_name = f"{molecule_name}_redox_inputs"
            redox_base_dir = None
            if base_dir is not None:
                if base_dir.endswith(redox_dir_name):
                    redox_base_dir = base_dir
                elif base_dir.endswith(f"{molecule_name}_psi4_inputs"):
                    parent = os.path.dirname(base_dir)
                    candidate = os.path.join(parent, redox_dir_name)
                    if os.path.exists(candidate):
                        redox_base_dir = candidate
                else:
                    candidate = os.path.join(base_dir, redox_dir_name)
                    if os.path.exists(candidate):
                        redox_base_dir = candidate
            if redox_base_dir is None:
                # Fallback: look in current working directory
                candidate = os.path.join(os.getcwd(), redox_dir_name)
                if os.path.exists(candidate):
                    redox_base_dir = candidate
            if redox_base_dir is None:
                print(f"[DEBUG] Could not find redox input directory for molecule '{molecule_name}' (tried base_dir: {base_dir})")
                redox_file_map = {k: None for k in ['neutral_gfec','neutral_sscf','anion_gfec','anion_sscf']}
            else:
                redox_file_map = build_redox_file_map_for_index(redox_base_dir, index, molecule_name)
            # Debug: print missing files if not all found
            if redox_file_map and all(redox_file_map.values()):
                reduction_gfe, redox_potential = calculate_redox_properties(
                    redox_file_map['neutral_sscf'],
                    redox_file_map['neutral_gfec'],
                    redox_file_map['anion_sscf'],
                    redox_file_map['anion_gfec']
                )
                if reduction_gfe is not None:
                    row['reduction_gfe'] = reduction_gfe
                if redox_potential is not None:
                    row['redoxpotential'] = redox_potential

        if row:
            row['ref'] = index

    except Exception as e:
        print(f"Error computing properties for index {index}: {e}")

    return row

def load_surface_points(surface_file):
    """Load surface points from the VDW surface file"""
    surface_points = []
    try:
        with open(surface_file, 'r') as f:
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 3:
                        x, y, z = float(parts[0]), float(parts[1]), float(parts[2])
                        surface_points.append([x, y, z])
    except Exception as e:
        print(f"Error reading surface file {surface_file}: {e}")
    return surface_points

def extract_all_properties_with_deltas(config_file, surface_file, base_dir, output_file, molecule_name, requested_props=None):
    """Extract properties for all indices and calculate deltas relative to reference"""
    
    # Check if base_dir is valid
    if base_dir is None:
        print(f"Error: base_dir is None. Cannot proceed with extraction.")
        return
    
    if not os.path.exists(base_dir):
        print(f"Error: base_dir '{base_dir}' does not exist.")
        return
    
    # Load configuration
    with open(config_file, 'r') as f:
        config = json.load(f)
    
    # Check for both standard and redox directories
    redox_dir = os.path.join(os.path.dirname(base_dir), f"{molecule_name}_redox_inputs")
    
    # Use requested properties parameter if provided, otherwise fall back to config
    if requested_props:
        props_requested = set(requested_props)
        print(f"Using requested properties from parameter: {sorted(props_requested)}")
    else:
        # Fall back to config-based detection
        config_props = set(config.get('properties', {}).keys()) if config.get('properties') else set()
        
        # If config only has redox but standard directory exists, we need to infer all properties
        if config_props == {'redox'} and os.path.exists(base_dir) and os.path.exists(redox_dir):
            print("Detected merged standard+redox workflow - extracting all standard properties plus redox")
            # This happens when both directories exist but config was overwritten by redox step
            props_requested = {
                'ground_state_energy', 'homo', 'lumo', 'gap', 'electron_affinity',
                'ionization_energy', 'electronegativity', 'hardness', 'electrophilicity',
                'nucleophilicity', 'chemical_potential', 'dipole_moment', 'redox'
            }
            # Add excited states that actually have outputs (check for TDDFT files)
            tddft_files = []
            for root, dirs, files in os.walk(base_dir):
                tddft_files.extend([f for f in files if 'neutral' in f and f.endswith('.out')])
            
            if tddft_files:
                # Conservatively add s1 only unless we detect higher states in outputs
                props_requested.add('s1_excitation_energy')
                props_requested.add('s1_oscillator_strength')
                
        else:
            # Use config properties as-is
            props_requested = config_props
        print(f"Properties to extract: {sorted(props_requested)}")
    
    # Load surface points
    surface_points = load_surface_points(surface_file)
    
    # Prepare results list
    results = []
    
    # Extract properties for reference (Ref)
    print("Extracting properties for reference (Ref)...")
    ref_file_map = build_file_map_for_index(base_dir, 'Ref', molecule_name)
    ref_props = compute_properties_for_index('Ref', ref_file_map, props_requested, base_dir, molecule_name)
    
    if ref_props:
        ref_row = {
            'ref': 'Ref',
            'potential_x': 0.0,
            'potential_y': 0.0,
            'potential_z': 0.0,
            **{k: v for k, v in ref_props.items() if k != 'ref'}
        }
        results.append(ref_row)
    else:
        print("Could not extract reference properties!")
        return
    
    # Extract properties for each surface point
    print(f"Extracting properties for {len(surface_points)} surface points...")
    for i, (x, y, z) in enumerate(surface_points):
        index_str = f"{i:04d}"
        print(f"Processing index {index_str}...")
        
        file_map = build_file_map_for_index(base_dir, index_str, molecule_name)
        props = compute_properties_for_index(index_str, file_map, props_requested, base_dir, molecule_name)
        
        if props:
            row = {
                'ref': index_str,
                'potential_x': x,
                'potential_y': y,
                'potential_z': z,
                **{k: v for k, v in props.items() if k != 'ref'}
            }
            results.append(row)
        else:
            print(f"Skipping index {index_str} due to missing data")
    
    # Create DataFrame
    if results:
        df = pd.DataFrame(results)
        
        # Reorder columns to have ref and potential coordinates first
        base_cols = ['ref', 'potential_x', 'potential_y', 'potential_z']
        prop_cols = [col for col in df.columns if col not in base_cols]
        df = df[base_cols + sorted(prop_cols)]
        
        # Get reference values for delta calculations
        ref_row = df[df['ref'] == 'Ref'].iloc[0]
        
        # Add delta columns for each property
        print("Calculating delta columns...")
        for prop_col in prop_cols:
            delta_col = f'delta_{prop_col}'
            ref_value = ref_row[prop_col]
            
            if pd.notna(ref_value):
                # Calculate delta: current_value - reference_value
                df[delta_col] = df[prop_col] - ref_value
                # Set delta for reference row to 0
                df.loc[df['ref'] == 'Ref', delta_col] = 0.0
            else:
                # If reference value is NaN, set all deltas to NaN
                df[delta_col] = np.nan
        
        # Reorder columns to interleave properties with their deltas
        final_cols = base_cols.copy()
        for prop_col in sorted(prop_cols):
            final_cols.append(prop_col)
            final_cols.append(f'delta_{prop_col}')
        
        df = df[final_cols]
        
        # Save to CSV
        df.to_csv(output_file, index=False)
        print(f"\nResults saved to {output_file}")
        print(f"Extracted properties for {len(results)} indices")
        print(f"Total columns: {len(df.columns)} (including {len(prop_cols)} delta columns)")
        
        # Display summary
        print(f"\nSummary:")
        print(f"- Reference (Ref): 1 row")
        print(f"- Surface points: {len(results)-1} rows")
        print(f"- Properties: {len(prop_cols)} columns")
        print(f"- Delta columns: {len(prop_cols)} columns")
        
    else:
        print("No properties extracted!")

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
    parser = argparse.ArgumentParser(description="Calculate molecular properties.")
    parser.add_argument('--index', '-idx', help="Charge index or 'Ref' to process (enables single index mode)")
    parser.add_argument('--outfolder', '-outf', help="Base folder containing neutral/anion/cation subfolders (auto-detected if not provided)")
    parser.add_argument('--config', '-cfg', default='psi4_params.json', help="JSON config path with properties")
    parser.add_argument('--outfile', '-of', default='properties.csv', help="CSV output filename (for single index mode)")
    parser.add_argument('--surface-file', '-sf', help="VDW surface file for all-indices extraction (auto-detected if not provided)")
    parser.add_argument('--excited-state', '-exs', type=int, default=1,
                        help='Number of excited states to extract (affects which excited state properties are extracted)')
    parser.add_argument('--requested-properties', '-req', nargs='+',
                        help='Specific properties requested by user (overrides config file)')
    args = parser.parse_args()
    
    # Propagate excited_state to extraction functions (no longer needed for dynamic functions)
    excitation_energy.excited_state = args.excited_state
    oscillator_strength.excited_state = args.excited_state

    # Try to get molecule name from config first, then auto-detect
    molecule_name = None
    if os.path.exists(args.config):
        try:
            with open(args.config) as f:
                config = json.load(f)
            if 'xyz' in config:
                molecule_name = config['xyz'].replace('.xyz', '')
                print(f"Using molecule from config: {molecule_name}")
        except:
            pass
    
    if not molecule_name:
        molecule_name = detect_molecule_name()
        print(f"Detected molecule: {molecule_name}")
    
    # Set defaults based on molecule name if not provided
    if not args.outfolder:
        args.outfolder = f"{molecule_name}_psi4_inputs"
    if not args.surface_file:
        args.surface_file = f"{molecule_name}_vdw_surface.txt"

    # If --index is provided, use single index mode (original functionality)
    if args.index:
        print(f"Single index mode: extracting properties for index {args.index}")
        
        with open(args.config) as f:
            config = json.load(f)

        props_requested = set(config['properties'].keys())

        file_map = build_file_map_for_index(args.outfolder, args.index, molecule_name)
        result_row = compute_properties_for_index(args.index, file_map, props_requested, args.outfolder, molecule_name)

        if not result_row:
            print(f"No properties computed for index {args.index}. Check if files exist.")
            return

        df = pd.DataFrame([result_row])

        # Only include reduction_gfe and redoxpotential if present in config properties
        ordered_cols = ['ref'] + [p for p in config['properties'].keys() if p in df.columns]
        df = df[ordered_cols]

        df.to_csv(args.outfile, index=False)
        print(f"Results saved to {args.outfile}")
        return

    # Default mode: Extract all properties with deltas
    print("Default mode: extracting properties for all indices with delta calculations")
    
    # Check if both standard and redox directories exist
    standard_dir = f"{molecule_name}_psi4_inputs"
    redox_dir = f"{molecule_name}_redox_inputs"
    
    if os.path.exists(standard_dir) and os.path.exists(redox_dir):
        print(f"Found both {standard_dir} and {redox_dir} - will extract from both")
        # Use standard directory as primary, extraction will automatically handle redox
        output_file = f'{molecule_name}_etm_properties.csv'
        extract_all_properties_with_deltas(args.config, args.surface_file, standard_dir, output_file, molecule_name, getattr(args, 'requested_properties', None))
    else:
        output_file = f'{molecule_name}_etm_properties.csv'
        extract_all_properties_with_deltas(args.config, args.surface_file, args.outfolder, output_file, molecule_name, getattr(args, 'requested_properties', None))

if __name__ == "__main__":
    main()
