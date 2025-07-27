import os
import sys
import json
import argparse
import subprocess
from tqdm import tqdm
from rdkit import Chem
from colorama import Fore, Style, init

init(autoreset=True)

DEFAULTS = {
    'vdw_scale': 1.0,
    'density': 1.0,
    'basis': '6-31G*',
    'method': 'b3lyp',
    'charge': 0,
    'multiplicity': 1,
    'properties': [
        'ground_state_energy', 'homo', 'lumo', 'gap', 'electron_affinity',
        'ionization_energy', 'electronegativity', 'hardness',
        'electrophilicity', 'nucleophilicity', 'chemical_potential',
        'dipole_moment', 'excitation_energy', 'oscillator_strength', 'redox'
    ] + [f's{i}_excitation_energy' for i in range(1, 11)] + [f's{i}_oscillator_strength' for i in range(1, 11)],
    'xyz': None,
    'solvent': {
        'enabled': False,
        'name': 'Water',
        'solver_type': 'CPCM',
        'radii_set': 'UFF',
        'cavity_type': 'GePol',
        'scaling': False,
        'area': 0.3,
        'mode': 'Implicit'
    },
    'redox_solvent': {
        'name': 'Acetonitrile',
        'solver_type': 'CPCM',
        'radii_set': 'UFF',
        'cavity_type': 'GePol',
        'scaling': False,
        'area': 0.3,
        'mode': 'Implicit'
    }
}

PROPERTIES = DEFAULTS['properties']

property_folders = {
    'ground_state_energy': {'neutral'},
    'homo': {'neutral'},
    'lumo': {'neutral'},
    'gap': {'neutral'},
    'electron_affinity': {'neutral', 'anion'},
    'ionization_energy': {'neutral', 'cation'},
    'chemical_potential': {'neutral', 'anion', 'cation'},
    'electronegativity': {'neutral', 'anion', 'cation'},
    'hardness': {'neutral', 'anion', 'cation'},
    'electrophilicity': {'neutral', 'anion', 'cation'},
    'nucleophilicity': {'neutral', 'anion', 'cation'},
    'dipole_moment': {'neutral'},
    'excitation_energy': {'neutral'},
    'oscillator_strength': {'neutral'},
    'redox': {'neutral_gfec', 'neutral_sscf', 'anion_gfec', 'anion_sscf'},
}

# Add dynamic excited state properties
for i in range(1, 11):
    property_folders[f's{i}_excitation_energy'] = {'neutral'}
    property_folders[f's{i}_oscillator_strength'] = {'neutral'}


def generate_vdw_surface(xyz_file, scale, density):
    print(Fore.YELLOW + f"\U0001F527 Running vsg: scale={scale}, density={density}")
    subprocess.run(['vsg', xyz_file, '--scale', str(scale), '--density', str(density), '--txt'], check=True)

def read_surface_points(txt_file):
    points = []
    with open(txt_file) as f:
        for line in f:
            parts = line.split()
            if len(parts) == 3:
                try:
                    points.append(tuple(map(float, parts)))
                except ValueError:
                    continue
    return points

def extract_atoms_from_xyz(xyz_file):
    mol = Chem.MolFromXYZFile(xyz_file)
    if mol is None:
        raise ValueError(f"RDKit could not parse {xyz_file}")
    conf = mol.GetConformer()
    return [(atom.GetSymbol(), conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z)
            for i, atom in enumerate(mol.GetAtoms())]

def construct_molecule_block(atoms, charge, multiplicity):
    lines = [f"molecule mol {{", f"  {charge} {multiplicity}"]
    lines += [f"  {s} {x:.6f} {y:.6f} {z:.6f}" for s, x, y, z in atoms]
    lines += ["  symmetry c1", "  no_reorient", "  no_com", "}\n"]
    return "\n".join(lines)

def construct_settings_block(basis, method, multiplicity, use_tddft=False):
    reference = 'uhf' if multiplicity != 1 else 'rks'
    if use_tddft:
        def settings_block_with_states(tdscf_states):
            return f"""set_options({{
  \"basis\": \"{basis}\",
  \"tdscf_states\": {tdscf_states},
  \"tdscf_maxiter\": 200,
  \"maxiter\": 200,
  \"diis\": True,
  \"diis_start\": 2,
  \"diis_max_vecs\": 10,
  \"level_shift\": 0.5,
  \"guess\": \"sad\",
  \"reference\": \"{reference}\"
}})"""
        return settings_block_with_states
    else:
        return f"""set {{
  basis {basis}
  scf_type df
  d_convergence 12
  maxiter 200
  diis true
  diis_start 2
  diis_max_vecs 10
  level_shift 0.5
  guess sad
  reference {reference}
}}"""

def construct_external_block(point):
    q, x, y, z = 1.0, *point
    return f"""external_potentials = np.array([
    [{q:.6f}, {x:.6f}, {y:.6f}, {z:.6f}]
])
external_potentials[:, 1:] /= constants.bohr2angstroms
"""

def construct_pcm_block(solvent_config):
    """Generate PCM solvation block"""
    if not solvent_config:
        return ""
    
    name = solvent_config.get('name', 'Water')
    solver_type = solvent_config.get('solver_type', 'CPCM')
    radii_set = solvent_config.get('radii_set', 'UFF')
    cavity_type = solvent_config.get('cavity_type', 'GePol')
    scaling = str(solvent_config.get('scaling', False)).lower()
    area = solvent_config.get('area', 0.3)
    mode = solvent_config.get('mode', 'Implicit')
    
    return f"""set_options({{
  'pcm': True,
  'pcm_scf_type': 'total',
}})

pcm_helper(\"\"\"
   Units = Angstrom
   Medium {{
   SolverType = {solver_type}
   Solvent = {name}
   }}
   Cavity {{
   RadiiSet = {radii_set}
   Type = {cavity_type}
   Scaling = {scaling}
   Area = {area}
   Mode = {mode}
   }}
\"\"\")
"""

def construct_property_block(properties, method, use_tddft=False, has_external=False):
    blocks = []
    ext_param = ", external_potentials=external_potentials" if has_external else ""
    
    if use_tddft:
        blocks.append(f"energy('td-{method}', molecule=mol{ext_param}, return_wfn=True)")
    else:
        blocks.append(f"energy('{method}', molecule=mol{ext_param}, return_wfn=True)")

        if 'dipole_moment' in properties:
            if has_external:
                blocks.append("dipole = properties('scf', properties=['dipole'], external_potentials=external_potentials)")
            else:
                blocks.append("dipole = properties('scf', properties=['dipole'])")

        if any(p in properties for p in ['homo', 'lumo', 'gap']):
            blocks.append(f"""
e, wfn = energy('{method}'{ext_param}, return_wfn=True)
homo = wfn.epsilon_a_subset('AO', 'ALL').np[wfn.nalpha() - 1]
lumo = wfn.epsilon_a_subset('AO', 'ALL').np[wfn.nalpha()]
gap = lumo - homo
""")

    return "\n".join(blocks)

def construct_frequency_block(properties, method, has_external=False):
    """Generate frequency calculation block for GFEC calculations"""
    blocks = []
    ext_param = ", external_potentials=external_potentials" if has_external else ""
    
    # First do SCF energy calculation
    blocks.append(f"energy('{method}', molecule=mol{ext_param}, return_wfn=True)")
    
    # Then add frequency calculation for Gibbs correction
    blocks.append(f"frequency('{method}/6-31g*', molecule=mol, return_wfn=True, dertype='gradient')")
    
    return "\n".join(blocks)

def write_psi4_input(path, mol_block, settings_block, ext_block, prop_block, pcm_block=None):
    with open(path, 'w') as f:
        f.write("import numpy as np\nfrom psi4 import constants\n\n")
        f.write(mol_block + "\n")
        if pcm_block:
            f.write(pcm_block + "\n")
        f.write(settings_block + "\n")
        if ext_block:
            f.write(ext_block + "\n")
        f.write(prop_block + "\n")

def write_all_inputs(atoms, mol_name, output_dir, config, surface_points):
    requested_props = set(config.get('properties', []))

    has_redox = 'redox' in requested_props
    non_redox_props = requested_props - {'redox'}
    has_non_redox = bool(non_redox_props)

    # If both redox and non-redox properties are present, generate both folders
    if has_redox:
        write_redox_inputs(atoms, mol_name, output_dir, config, surface_points)
    if has_non_redox:
        write_standard_inputs(atoms, mol_name, output_dir, config, surface_points)
    if not has_redox and not has_non_redox:
        print(Fore.RED + '\u274c No valid properties specified for input generation.')
        sys.exit(1)

def write_standard_inputs(atoms, mol_name, output_dir, config, surface_points):
    """Write standard ETM input files"""
    os.makedirs(output_dir, exist_ok=True)

    charge_anion = config.get('charge_anion', config['charge'] - 1)
    multiplicity_anion = config.get('multiplicity_anion', config['multiplicity'] + 1)
    charge_cation = config.get('charge_cation', config['charge'] + 1)
    multiplicity_cation = config.get('multiplicity_cation', config['multiplicity'] + 1)

    job_map = {
        'neutral': {
            'charge': config['charge'],
            'multiplicity': config['multiplicity'],
            'props': [
                'ground_state_energy', 'homo', 'lumo', 'gap',
                'dipole_moment', 'excitation_energy', 'oscillator_strength',
                'electrophilicity', 'nucleophilicity', 'electronegativity', 'hardness', 'chemical_potential'
            ] + [f's{i}_excitation_energy' for i in range(1, 11)] + [f's{i}_oscillator_strength' for i in range(1, 11)],
            'desc': 'neutral molecule'
        },
        'anion': {
            'charge': charge_anion,
            'multiplicity': multiplicity_anion,
            'props': ['electron_affinity', 'electronegativity', 'hardness', 'electrophilicity', 'chemical_potential'],
            'desc': 'anion'
        },
        'cation': {
            'charge': charge_cation,
            'multiplicity': multiplicity_cation,
            'props': ['ionization_energy', 'electronegativity', 'hardness', 'electrophilicity', 'chemical_potential'],
            'desc': 'cation'
        }
    }

    requested_props = set(config.get('properties', []))
    job_types = [job for job, data in job_map.items() if not requested_props or any(p in requested_props for p in data['props'])]

    if not job_types:
        print(Fore.RED + "\u274c No valid job types based on requested properties.")
        sys.exit(1)

    # Check if optional solvent is enabled
    solvent_config = config.get('solvent', {})
    use_solvent = solvent_config.get('enabled', False)
    pcm_block = construct_pcm_block(solvent_config) if use_solvent else None

    for job_type in job_types:
        job = job_map[job_type]
        folder = os.path.join(output_dir, job_type)
        os.makedirs(folder, exist_ok=True)

        props = [p for p in requested_props if p in job['props']] if requested_props else job['props']
        if not props:
            print(Fore.YELLOW + f"\u26a0\ufe0f Skipping {job_type}; no relevant properties requested.")
            continue

        excited_state_props = [f's{i}_excitation_energy' for i in range(1, 11)] + [f's{i}_oscillator_strength' for i in range(1, 11)]
        use_tddft = job_type == 'neutral' and any(p in props for p in ['excitation_energy', 'oscillator_strength'] + excited_state_props)
        input_type = 'TDDFT' if use_tddft else 'SCF'
        
        solvent_msg = " (with solvent)" if use_solvent else ""
        print(Fore.YELLOW + f"\u270d\ufe0f Writing {input_type} inputs for {job['desc']} ({', '.join(props)}){solvent_msg}")
        
        mol_block = construct_molecule_block(atoms, job['charge'], job['multiplicity'])
        if use_tddft:
            settings_block_fn = construct_settings_block(config['basis'], config['method'], job['multiplicity'], use_tddft)
            settings_block = settings_block_fn(config.get('excited_state', 1))
        else:
            settings_block = construct_settings_block(config['basis'], config['method'], job['multiplicity'], use_tddft)
        prop_block = construct_property_block(props, config['method'], use_tddft, has_external=False)

        suffix = 'Ref'
        write_psi4_input(os.path.join(folder, f"{mol_name}_{job_type}_{suffix}.in"),
                         mol_block, settings_block, None, prop_block, pcm_block)

        for idx, point in enumerate(tqdm(surface_points, desc=f"Generating {job['desc']} charge inputs")):
            ext_block = construct_external_block(point)
            prop_block_ext = construct_property_block(props, config['method'], use_tddft, has_external=True)
            charge_path = os.path.join(folder, f"{mol_name}_{job_type}_charge_{idx:04d}.in")
            if use_tddft:
                settings_block = settings_block_fn(config.get('excited_state', 1))
            else:
                settings_block = construct_settings_block(config['basis'], config['method'], job['multiplicity'], use_tddft)
            write_psi4_input(charge_path, mol_block, settings_block, ext_block, prop_block_ext, pcm_block)

    print(Fore.GREEN + f"\u2705 All Psi4 input files saved to: {output_dir}")

def write_redox_inputs(atoms, mol_name, output_dir, config, surface_points):
    """Write redox calculation input files"""
    # Use redox-specific directory name
    redox_output_dir = output_dir.replace('_psi4_inputs', '_redox_inputs')
    os.makedirs(redox_output_dir, exist_ok=True)
    
    charge_anion = config.get('charge_anion', config['charge'] - 1)
    multiplicity_anion = config.get('multiplicity_anion', config['multiplicity'] + 1)
    
    # Get redox solvent configuration (mandatory for SSCF)
    redox_solvent_config = config.get('redox_solvent', {'name': 'Acetonitrile', 'solver_type': 'CPCM'})
    pcm_block = construct_pcm_block(redox_solvent_config)
    

    # Define new flat folder structure for redox jobs
    redox_jobs = {
        'neutral_gfec': {
            'charge': config['charge'],
            'multiplicity': config['multiplicity'],
            'calc_type': 'gfec',
            'desc': 'neutral gas-phase frequency'
        },
        'neutral_sscf': {
            'charge': config['charge'],
            'multiplicity': config['multiplicity'],
            'calc_type': 'sscf',
            'desc': 'neutral solvated SCF'
        },
        'anion_gfec': {
            'charge': charge_anion,
            'multiplicity': multiplicity_anion,
            'calc_type': 'gfec',
            'desc': 'anion gas-phase frequency'
        },
        'anion_sscf': {
            'charge': charge_anion,
            'multiplicity': multiplicity_anion,
            'calc_type': 'sscf',
            'desc': 'anion solvated SCF'
        }
    }

    print(Fore.CYAN + f"\u2699\ufe0f  Starting redox calculation input generation")
    print(Fore.CYAN + f"Redox solvent: {redox_solvent_config['name']} ({redox_solvent_config['solver_type']})")

    for job_name, job in redox_jobs.items():
        # Create flat folder structure: anion_gfec, anion_sscf, neutral_gfec, neutral_sscf
        folder = os.path.join(redox_output_dir, job_name)
        os.makedirs(folder, exist_ok=True)

        print(Fore.YELLOW + f"\u270d\ufe0f Writing {job['desc']} inputs")

        mol_block = construct_molecule_block(atoms, job['charge'], job['multiplicity'])
        settings_block = construct_settings_block(config['basis'], config['method'], job['multiplicity'], use_tddft=False)

        # Choose calculation type
        if job['calc_type'] == 'gfec':
            # Gas-phase frequency calculation (no PCM)
            prop_block = construct_frequency_block(['redox'], config['method'], has_external=False)
            job_pcm_block = None
        else:  # sscf
            # Solvated SCF calculation (with PCM)
            prop_block = construct_property_block(['redox'], config['method'], use_tddft=False, has_external=False)
            job_pcm_block = pcm_block

        # Reference calculation (no external charge)
        ref_file = os.path.join(folder, f"{mol_name}_{job_name}_Ref.in")
        write_psi4_input(ref_file, mol_block, settings_block, None, prop_block, job_pcm_block)

        # Surface point calculations
        for idx, point in enumerate(tqdm(surface_points, desc=f"Generating {job['desc']} charge inputs")):
            ext_block = construct_external_block(point)

            if job['calc_type'] == 'gfec':
                prop_block_ext = construct_frequency_block(['redox'], config['method'], has_external=True)
                job_pcm_block = None
            else:  # sscf
                prop_block_ext = construct_property_block(['redox'], config['method'], use_tddft=False, has_external=True)
                job_pcm_block = pcm_block

            charge_path = os.path.join(folder, f"{mol_name}_{job_name}_charge_{idx:04d}.in")
            write_psi4_input(charge_path, mol_block, settings_block, ext_block, prop_block_ext, job_pcm_block)

    print(Fore.GREEN + f"\u2705 All redox input files saved to: {redox_output_dir}")

def build_property_folder_mapping(requested_props):
    result = {}
    for prop in requested_props:
        folders = property_folders.get(prop, set())
        folder_str = ", ".join(sorted(folders))
        result[prop] = folder_str
    return result

def save_config(path, config):
    requested_props = set(config.get('properties', []))

    config_to_save = dict(config)
    config_to_save['properties'] = build_property_folder_mapping(requested_props)
    
    with open(path, 'w') as f:
        json.dump(config_to_save, f, indent=2)
    print(Fore.GREEN + f"\U0001F4BE Saved config to {path}")

def load_or_create_config(path, cli_args):
    if os.path.exists(path):
        try:
            with open(path, 'r') as f:
                existing_config = json.load(f)
            print(Fore.YELLOW + f"\u26a0\ufe0f Updating existing config: {path}")
        except:
            print(Fore.YELLOW + f"\u26a0\ufe0f Replacing corrupted config: {path}")

    config = DEFAULTS.copy()

    if cli_args.xyz:
        config['xyz'] = cli_args.xyz

    for key in ['vdw_scale', 'density', 'basis', 'method', 'charge', 'multiplicity',
                'charge_anion', 'multiplicity_anion', 'charge_cation', 'multiplicity_cation']:
        val = getattr(cli_args, key, None)
        if val is not None:
            config[key] = val

    if cli_args.properties:
        config['properties'] = cli_args.properties

    # Handle standard ETM solvent configuration
    if cli_args.solvent:
        config['solvent']['enabled'] = True
        config['solvent']['name'] = cli_args.solvent
        if hasattr(cli_args, 'solvent_type') and cli_args.solvent_type:
            config['solvent']['solver_type'] = cli_args.solvent_type  
        if hasattr(cli_args, 'solvent_radii') and cli_args.solvent_radii:
            config['solvent']['radii_set'] = cli_args.solvent_radii

    # Handle redox solvent configuration  
    if hasattr(cli_args, 'redox_solvent') and cli_args.redox_solvent:
        config['redox_solvent']['name'] = cli_args.redox_solvent
    if hasattr(cli_args, 'redox_solvent_type') and cli_args.redox_solvent_type:
        config['redox_solvent']['solver_type'] = cli_args.redox_solvent_type
    if hasattr(cli_args, 'redox_solvent_radii') and cli_args.redox_solvent_radii:
        config['redox_solvent']['radii_set'] = cli_args.redox_solvent_radii

    if not config.get('xyz'):
        print(Fore.RED + "\u274c No input .xyz file provided.")
        sys.exit(1)

    save_config(path, config)
    return config

def parse_args():
    parser = argparse.ArgumentParser(description="Generate Psi4 inputs with VDW point charges.")
    parser.add_argument('--config', '-cfg', type=str, default="psi4_params.json", help="JSON config path.")
    parser.add_argument('xyz', help="Input .xyz file.")
    parser.add_argument('--vdw_scale', '-vds', type=float)
    parser.add_argument('--density', '-dns', type=float)
    parser.add_argument('--basis', '-bss', type=str)
    parser.add_argument('--method', '-mtd', type=str)
    parser.add_argument('--charge', '-chg', type=int)
    parser.add_argument('--multiplicity', '-mult', type=int)
    parser.add_argument('--properties', '-prop', nargs='+', choices=PROPERTIES)
    parser.add_argument('--charge_anion', '-chga', type=int)
    parser.add_argument('--multiplicity_anion', '-multa', type=int)
    parser.add_argument('--charge_cation', '-chgc', type=int)
    parser.add_argument('--multiplicity_cation', '-multc', type=int)
    parser.add_argument('--excited-state', '-exs', type=int, default=1,
                        help='Number of excited states to calculate (tdscf_states)')
    parser.add_argument('--solvent', '-sol', type=str, help='Enable solvent for ETM calculations (e.g., Water, Acetonitrile)')
    parser.add_argument('--solvent-type', '-solt', type=str, default='CPCM', help='Solvent model type (CPCM, IEFPCM)')
    parser.add_argument('--solvent-radii', '-solr', type=str, default='UFF', help='Radii set for solvent (UFF, BONDI)')
    parser.add_argument('--redox-solvent', '-rsol', type=str, help='Solvent for redox SSCF calculations (default: Acetonitrile)')
    parser.add_argument('--redox-solvent-type', '-rsolt', type=str, default='CPCM', help='Redox solvent model type')
    parser.add_argument('--redox-solvent-radii', '-rsolr', type=str, default='UFF', help='Redox solvent radii set')
    
    return parser.parse_args()

def detect_molecule_name():
    """Auto-detect molecule name from available files."""
    # Look for XYZ files (excluding surface files)
    xyz_files = [f for f in os.listdir('.') if f.endswith('.xyz') and 'surface' not in f.lower()]
    
    if not xyz_files:
        print(Fore.RED + "Error: No molecule XYZ file found!")
        print("Please ensure you have a molecule XYZ file (e.g., water.xyz, benzene.xyz)")
        sys.exit(1)
    
    if len(xyz_files) > 1:
        print(Fore.YELLOW + f"Warning: Multiple XYZ files found: {xyz_files}")
        print(f"Using: {xyz_files[0]}")
    
    molecule_name = xyz_files[0].replace('.xyz', '')
    return molecule_name

def main():
    args = parse_args()
    
    # If no XYZ file provided as argument, try to detect it
    if not args.xyz:
        try:
            molecule_name = detect_molecule_name()
            args.xyz = f"{molecule_name}.xyz"
            print(Fore.CYAN + f"Auto-detected molecule: {molecule_name}")
        except SystemExit:
            # detect_molecule_name calls sys.exit(1) if no XYZ found
            raise
    
    config = load_or_create_config(args.config, args)
    config['excited_state'] = args.excited_state
    mol_name = os.path.splitext(os.path.basename(config['xyz']))[0]

    print(Fore.BLUE + f"\U0001F4C4 Molecule: {mol_name}")
    atoms = extract_atoms_from_xyz(config['xyz'])

    print(Fore.YELLOW + "\U0001F527 Generating VDW surface ...")
    generate_vdw_surface(config['xyz'], config['vdw_scale'], config['density'])

    surface_file = f"{mol_name}_vdw_surface.txt"
    print(Fore.MAGENTA + f"\U0001F4E5 Reading surface points from: {surface_file}")
    surface_points = read_surface_points(surface_file)

    # Always use the standard folder as the base output_dir for write_all_inputs
    # The function itself will create both folders if needed
    output_dir = f"{mol_name}_psi4_inputs"
    write_all_inputs(atoms, mol_name, output_dir, config, surface_points)

if __name__ == "__main__":
    main()
