#!/usr/bin/env python3
"""
ETM Suite Main Pipeline
======================

Automated workflow for running the complete ETM Suite (Electrostatic Tuning Map Suite) analysis.
This script orchestrates all four steps of the ETM Suite workflow:

1. Generate Psi4 input files with external point charges
2. Submit and run quantum chemistry calculations  
3. Extract molecular properties from output files
4. Generate MOL2 visualization files

Usage:
    python main.py [molecule.xyz] [options]

Examples:
    python main.py                           # Auto-detect molecule from .xyz files
    python main.py benzene.xyz               # Use specific molecule
    python main.py water.xyz --threads 8     # Use 8 threads total
    python main.py --help                    # Show all options
"""

import os
import sys
import time
import json
import argparse
import subprocess
from pathlib import Path
from colorama import Fore, Style, init

# Initialize colorama for colored output
init(autoreset=True)

def print_header(text):
    """Print a styled header"""
    print(f"\n{Fore.CYAN}{'='*60}")
    print(f"{Fore.CYAN}{text:^60}")
    print(f"{Fore.CYAN}{'='*60}{Style.RESET_ALL}")

def print_step(step_num, description):
    """Print a step header"""
    print(f"\n{Fore.YELLOW}🚀 Step {step_num}: {description}{Style.RESET_ALL}")
    print(f"{Fore.YELLOW}{'-'*50}{Style.RESET_ALL}")

def print_success(message):
    """Print a success message"""
    print(f"{Fore.GREEN}✅ {message}{Style.RESET_ALL}")

def print_error(message):
    """Print an error message"""
    print(f"{Fore.RED}❌ {message}{Style.RESET_ALL}")

def print_warning(message):
    """Print a warning message"""
    print(f"{Fore.YELLOW}⚠️  {message}{Style.RESET_ALL}")

def run_command(cmd, description, capture_output=False):
    """Run a shell command with error handling"""
    print(f"{Fore.BLUE}💻 Running: {' '.join(cmd)}{Style.RESET_ALL}")
    
    try:
        if capture_output:
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)
            return result.stdout, result.stderr
        else:
            subprocess.run(cmd, check=True)
            return None, None
    except subprocess.CalledProcessError as e:
        print_error(f"Command failed with exit code {e.returncode}")
        if hasattr(e, 'stderr') and e.stderr:
            print(f"{Fore.RED}Error output: {e.stderr}{Style.RESET_ALL}")
        return None, e
    except FileNotFoundError:
        print_error(f"Command not found: {cmd[0]}")
        print_warning("Make sure the ETM Suite scripts are in the current directory")
        return None, FileNotFoundError(f"Command not found: {cmd[0]}")

def load_local_config(config_path='psi4_params.json'):
    """Load local execution settings from JSON config file"""
    if not config_path or not os.path.exists(config_path):
        return {}
    
    try:
        with open(config_path) as f:
            config = json.load(f)
        return config.get("local", {})
    except Exception as e:
        print_warning(f"Could not load local config from {config_path}: {e}")
        return {}

def check_prerequisites():
    """Check if required files and dependencies are available"""
    print_step("0", "Checking Prerequisites")
    
    # Check for required scripts
    required_scripts = [
        "etm-inputs.py",
        "etm-submit.py", 
        "etm-extract.py",
        "etm-construct.py"
    ]
    
    missing_scripts = []
    for script in required_scripts:
        if not os.path.exists(script):
            missing_scripts.append(script)
    
    if missing_scripts:
        print_error(f"Missing required scripts: {', '.join(missing_scripts)}")
        return False
    
    # Check for psi4_params.json and create default if missing
    if not os.path.exists('psi4_params.json'):
        print(f"{Fore.YELLOW}⚠️  psi4_params.json not found{Style.RESET_ALL}")
        print(f"{Fore.BLUE}📋 Creating default configuration...{Style.RESET_ALL}")
        if not create_default_psi4_params():
            print_warning("Failed to create default psi4_params.json")
            print_warning("Some downstream scripts may fail without this configuration")
        else:
            print(f"{Fore.GREEN}✅ Default psi4_params.json created successfully{Style.RESET_ALL}")
    else:
        print(f"{Fore.GREEN}✅ Found psi4_params.json configuration{Style.RESET_ALL}")
    
    print_success("All prerequisites satisfied")
    return True

def step1_generate_inputs(molecule_xyz=None, args=None):
    """Step 1: Generate Psi4 input files"""
    print_step("1", "Generate Psi4 Input Files")
    
    cmd = ["python", "etm-inputs.py"]
    if molecule_xyz:
        cmd.append(molecule_xyz)
    if args:
        cmd.extend(build_inputs_args(args))
    
    stdout, stderr = run_command(cmd, "Generating Psi4 input files")
    
    if stderr:
        return False
    
    print_success("Psi4 input files generated successfully")
    return True

def step2_submit_jobs(molecule_name, mode="local", threads_per_job=2, total_threads=4, 
                     clean_outputs=True, args=None):
    """Step 2: Submit and run quantum chemistry calculations"""
    print_step("2", "Submit and Run Quantum Chemistry Calculations")
    
    # Submit jobs for both input directories if present
    redox_input_dir = f"{molecule_name}_redox_inputs"
    standard_input_dir = f"{molecule_name}_psi4_inputs"
    input_dirs = []
    if os.path.exists(redox_input_dir):
        input_dirs.append((redox_input_dir, "redox"))
    if os.path.exists(standard_input_dir):
        input_dirs.append((standard_input_dir, "standard"))

    if not input_dirs:
        print_error(f"No input directory found: neither {redox_input_dir} nor {standard_input_dir}")
        return False

    all_success = True
    for input_dir, label in input_dirs:
        if label == "redox":
            print(f"{Fore.BLUE}📋 Using redox input directory: {input_dir}{Style.RESET_ALL}")
        else:
            print(f"{Fore.BLUE}📋 Using standard input directory: {input_dir}{Style.RESET_ALL}")
        cmd = [
            "python", "etm-submit.py",
            "--input", input_dir,
            "--threads_per_job", str(threads_per_job),
            "--total_threads", str(total_threads)
        ]
        if clean_outputs:
            cmd.append("--clean_outputs")
        if args and hasattr(args, 'force') and args.force:
            cmd.append("--force")
        stdout, stderr = run_command(cmd, f"Running quantum chemistry calculations for {label} inputs")
        if stderr:
            all_success = False
    if all_success:
        print_success("Quantum chemistry calculations completed successfully")
        return True
    else:
        return False

def step3_extract_properties(args=None, molecule_name=None):
    """Step 3: Extract molecular properties from output files"""
    print_step("3", "Extract Molecular Properties")
    
    cmd = ["python", "etm-extract.py"]
    if args:
        cmd.extend(build_extract_args(args))
    
    stdout, stderr = run_command(cmd, "Extracting molecular properties")
    
    if stderr:
        return False
    
    # Check if properties file was created
    csv_file = f"{molecule_name}_etm_properties.csv" if molecule_name else "etm_properties.csv"
    if os.path.exists(csv_file):
        print_success("Molecular properties extracted successfully")
        print(f"{Fore.BLUE}📊 Properties saved to: {csv_file}{Style.RESET_ALL}")
        return True
    else:
        print_error(f"Properties file not created: {csv_file}")
        return False

def step4_generate_mol2(keep_inputs=False, custom_args=None, molecule_name=None, args=None):
    """Step 4: Generate MOL2 visualization files"""
    print_step("4", "Generate MOL2 Visualization Files")
    
    cmd = ["python", "etm-construct.py"]
    if keep_inputs:
        cmd.extend(["-ki", "--keep-inputs"])
    # Add excited state arguments for construction
    if args and hasattr(args, 'requested_excited_states') and args.requested_excited_states:
        # Pass specific excited states requested by user for construction
        for state in args.requested_excited_states:
            cmd.extend(['--excited-state', str(state)])
    elif args and hasattr(args, 'excited_state'):
        # Fallback to original behavior if no specific states requested
        cmd.extend(['--excited-state', str(args.excited_state)])
    if custom_args:
        cmd.extend(custom_args)
    
    stdout, stderr = run_command(cmd, "Generating MOL2 files")
    
    if stderr:
        return False
    
    # Check if MOL2 files were created
    results_dir = f"{molecule_name}_etm_results" if molecule_name else "etm-results"
    if os.path.exists(results_dir) and os.listdir(results_dir):
        mol2_files = [f for f in os.listdir(results_dir) if f.endswith(".mol2")]
        print_success(f"MOL2 files generated successfully ({len(mol2_files)} files)")
        print(f"{Fore.BLUE}📁 MOL2 files saved to: {results_dir}/{Style.RESET_ALL}")
        return True
    else:
        print_error(f"MOL2 files not created in: {results_dir}")
        return False

def print_summary(molecule_name, start_time, keep_inputs=False):
    """Print a summary of the completed workflow"""
    end_time = time.time()
    duration = end_time - start_time
    
    print_header("ETM Suite Workflow Summary")
    print(f"{Fore.GREEN}🎉 ETM Suite workflow completed successfully!{Style.RESET_ALL}")
    print(f"{Fore.BLUE}📋 Molecule: {molecule_name}{Style.RESET_ALL}")
    print(f"{Fore.BLUE}⏱️  Total time: {duration:.1f} seconds{Style.RESET_ALL}")
    
    print(f"\n{Fore.CYAN}📂 Generated Files:{Style.RESET_ALL}")
    
    # List generated files
    files_to_check = [
        (f"{molecule_name}_vdw_surface.txt", "VDW surface coordinates"),
        (f"{molecule_name}_etm_properties.csv", "Extracted molecular properties"),
        (f"{molecule_name}_etm_results/", "MOL2 visualization files directory")
    ]
    
    # Handle input directory based on keep_inputs setting
    inputs_dir = f"{molecule_name}_psi4_inputs/"
    if keep_inputs:
        files_to_check.insert(1, (inputs_dir, "Psi4 input files directory"))
    else:
        # Check if directory exists (might not have been cleaned up yet)
        if os.path.exists(inputs_dir):
            files_to_check.insert(1, (inputs_dir, "Psi4 input files directory"))
        else:
            # Directory was cleaned up - show positive message
            print(f"  ✅ Psi4 input files directory - removed as requested")
            files_to_check.insert(1, (None, None))  # Placeholder to maintain order
    
    for filepath, description in files_to_check:
        if filepath is None:  # Skip placeholder
            continue
        if os.path.exists(filepath):
            if os.path.isdir(filepath):
                count = len([f for f in os.listdir(filepath) if not f.startswith('.')])
                print(f"  ✅ {filepath} ({count} files) - {description}")
            else:
                print(f"  ✅ {filepath} - {description}")
        else:
            print(f"  ❌ {filepath} - {description} (not found)")
    
    print(f"\n{Fore.YELLOW}🚀 Next Steps:{Style.RESET_ALL}")
    print(f"  1. Open the MOL2 files in {molecule_name}_etm_results/ with a molecular viewer (VMD, ChimeraX, etc.)")
    print(f"  2. Analyze the {molecule_name}_etm_properties.csv file for quantitative data")
    print("  3. Visualize property variations across the molecular surface")

def create_full_workflow_script(molecule_name, args, script_name="etm_workflow.sh"):
    """Create a script that runs the complete workflow including job submission and waiting"""
    
    # Remove any existing workflow scripts to avoid confusion
    existing_workflows = ["etm_workflow.sh", "etm_local_workflow.sh"]
    for old_script in existing_workflows:
        if os.path.exists(old_script):
            os.remove(old_script)
    
    # Build the etm-submit.py command
    submit_cmd = [
        "python", "etm-submit.py",
        "--input", f"{molecule_name}_psi4_inputs",
        "--threads_per_job", str(args.threads_per_job),
        "--total_threads", str(args.threads)
    ]
    
    if args.clean_outputs:
        submit_cmd.append("--clean_outputs")
    
    # Add local arguments
    if args.force:
        submit_cmd.append('--force')
    
    script_content = f"""#!/bin/bash
# ETM Suite Full Workflow Script (Detached Mode)
# Auto-generated on {time.strftime('%Y-%m-%d %H:%M:%S')}

# Set up signal handling to prevent early termination
set -euo pipefail
trap 'echo "Workflow interrupted at $(date)"' INT TERM

# Redirect all output to log file while still showing progress
exec > >(tee -a etm_workflow.log) 2>&1

echo "=== ETM Suite Full Workflow for {molecule_name} ==="
echo "Start time: $(date)"
echo "Process ID: $$"
echo "Working directory: $(pwd)"
echo ""

# Step 2: Submit and wait for calculations to complete
echo "Step 2: Running quantum chemistry calculations..."
{' '.join(submit_cmd)}

EXIT_CODE=$?
if [ $EXIT_CODE -ne 0 ]; then
    echo "WARNING: Some calculations failed (exit code: $EXIT_CODE)"
    if [ -f "failed_jobs.csv" ]; then
        FAILED_COUNT=$(tail -n +2 failed_jobs.csv | wc -l)
        echo "Number of failed jobs: $FAILED_COUNT"
        echo "Failed jobs logged in failed_jobs.csv"
        
        # Check if we have any successful output files to continue with
        TOTAL_OUTPUTS=$(find {molecule_name}_psi4_inputs -name "*.out" 2>/dev/null | wc -l)
        if [ $TOTAL_OUTPUTS -gt 10 ]; then
            echo "Found $TOTAL_OUTPUTS successful outputs - continuing with extraction..."
        else
            echo "ERROR: Too few successful calculations ($TOTAL_OUTPUTS outputs)"
            echo "Cannot proceed with meaningful analysis"
            exit 1
        fi
    else
        echo "ERROR: All calculations failed"
        echo "Check the output above for details"
        exit 1
    fi
else
    echo "All calculations completed successfully"
fi

echo ""
echo "Step 2 completed successfully at $(date)"
echo ""

# Step 3: Extract properties
echo "Step 3: Extracting molecular properties..."
python etm-extract.py
if [ $? -ne 0 ]; then
    echo "ERROR: Property extraction failed"
    echo "Check the output above for details"
    exit 1
fi

# Verify extraction completed
if [ ! -f "{molecule_name}_etm_properties.csv" ]; then
    echo "ERROR: {molecule_name}_etm_properties.csv not created"
    exit 1
fi
echo "✅ Properties extracted successfully: {molecule_name}_etm_properties.csv"

echo "Step 3 completed successfully at $(date)"
echo ""

# Step 4: Generate MOL2 files
echo "Step 4: Generating MOL2 visualization files..."
python etm-construct.py {'-ki --keep-inputs' if args.keep_inputs else ''}
if [ $? -ne 0 ]; then
    echo "ERROR: MOL2 generation failed"
    echo "Check the output above for details"
    exit 1
fi

# Verify MOL2 files were created
if [ ! -d "{molecule_name}_etm_results" ]; then
    echo "ERROR: {molecule_name}_etm_results directory not created"
    exit 1
fi
echo "✅ MOL2 files generated successfully in {molecule_name}_etm_results/"

echo "Step 4 completed successfully at $(date)"
echo ""

# Final cleanup: Remove Psi4 cleanup manifest files
echo "Cleaning up Psi4 temporary files..."
rm -f psi.*.clean 2>/dev/null || true
echo "Cleanup completed"
echo ""

echo "=== ETM Suite workflow completed! ==="
echo "Results available in {molecule_name}_etm_results/ directory"
echo "Properties saved to {molecule_name}_etm_properties.csv"

# Show summary of any failed jobs
if [ -f "failed_jobs.csv" ]; then
    FAILED_COUNT=$(tail -n +2 failed_jobs.csv | wc -l)
    if [ $FAILED_COUNT -gt 0 ]; then
        echo "⚠️  $FAILED_COUNT jobs failed due to convergence issues (see failed_jobs.csv)"
        echo "Consider using fix_convergence.py for problematic surface points"
    fi
fi

echo "Completion timestamp: $(date)"
"""
    
    with open(script_name, 'w') as f:
        f.write(script_content)
    
    os.chmod(script_name, 0o755)  # Make executable
    print(f"{Fore.BLUE}📝 Created full workflow script: {script_name}{Style.RESET_ALL}")
    return script_name

def step2_submit_jobs_detached(molecule_name, args):
    """Submit jobs in detached mode"""
    print_step("2", "Submit and Run Quantum Chemistry Calculations (DETACHED)")
    
    input_dir = f"{molecule_name}_psi4_inputs"
    
    if not os.path.exists(input_dir):
        print_error(f"Input directory not found: {input_dir}")
        return False
    
    if args.mode == 'local':
        # For local mode: create and run full workflow script in background
        workflow_script = create_full_workflow_script(molecule_name, args, "etm_workflow.sh")
        
        cmd = ["nohup", f"./{workflow_script}"]
        
        print(f"{Fore.BLUE}💻 Starting detached local workflow...{Style.RESET_ALL}")
        print(f"{Fore.YELLOW}📋 Progress will be logged to: nohup.out and etm_workflow.log{Style.RESET_ALL}")
        print(f"{Fore.YELLOW}📋 Monitor with: tail -f nohup.out or tail -f etm_workflow.log{Style.RESET_ALL}")
        print(f"{Fore.YELLOW}📋 Check if running: ps aux | grep etm{Style.RESET_ALL}")
        print(f"{Fore.YELLOW}📋 Cancel if needed: pkill -f etm{Style.RESET_ALL}")
        print(f"{Fore.YELLOW}📋 Check completion: ls -la {molecule_name}_etm_results/ {molecule_name}_etm_properties.csv{Style.RESET_ALL}")
        
        try:
            # Start the background process with proper nohup redirection
            with open('nohup.out', 'w') as f:
                subprocess.Popen(cmd, stdout=f, stderr=subprocess.STDOUT, 
                               preexec_fn=os.setsid)
            print_success("Background workflow started - you can safely close this terminal")
            return True
        except Exception as e:
            print_error(f"Failed to start background process: {e}")
            return False
    
    return False

def create_default_psi4_params(config_path='psi4_params.json', molecule_xyz=None):
    """Create a default psi4_params.json file with sensible defaults"""
    
    # Auto-detect molecule file if not provided
    if not molecule_xyz:
        xyz_files = [f for f in os.listdir('.') if f.endswith('.xyz') and 'surface' not in f.lower()]
        if xyz_files:
            # Prefer LF.xyz if available, otherwise use first available
            if 'LF.xyz' in xyz_files:
                molecule_xyz = 'LF.xyz'
            else:
                molecule_xyz = xyz_files[0]
        else:
            molecule_xyz = "molecule.xyz"  # Default fallback
    
    default_config = {
        "vdw_scale": 1.0,
        "density": 1.0,
        "basis": "6-31G*",
        "method": "b3lyp",
        "charge": 0,
        "multiplicity": 1,
        "properties": {
            "homo": "neutral",
            "ground_state_energy": "neutral",
            "electronegativity": "anion, cation, neutral",
            "electrophilicity": "anion, cation, neutral",
            "gap": "neutral",
            "dipole_moment": "neutral",
            "electron_affinity": "anion, neutral",
            "excitation_energy": "neutral",
            "oscillator_strength": "neutral",
            "hardness": "anion, cation, neutral",
            "chemical_potential": "anion, cation, neutral",
            "lumo": "neutral",
            "ionization_energy": "cation, neutral",
            "nucleophilicity": "anion, cation, neutral"
        },
        "xyz": molecule_xyz,
        "local": {
            "threads_per_job": 2,
            "total_threads": 4
        }
    }
    
    try:
        with open(config_path, 'w') as f:
            json.dump(default_config, f, indent=2)
        print(f"{Fore.GREEN}✅ Created default configuration: {config_path}{Style.RESET_ALL}")
        print(f"{Fore.BLUE}📋 Using molecule file: {molecule_xyz}{Style.RESET_ALL}")
        print(f"{Fore.BLUE}📋 Local execution: threads_per_job=2, total_threads=4{Style.RESET_ALL}")
        return True
    except Exception as e:
        print_warning(f"Could not create default config {config_path}: {e}")
        return False

def build_inputs_args(args):
    """Use optimized argument builder"""
    from argument_builder import build_etm_inputs_args
    return build_etm_inputs_args(args)

def build_submit_args(args):
    """Build command-line arguments for etm-submit.py (local mode only)"""
    cmd_args = []
    
    # For local mode, only add force if specified
    if args.force:
        cmd_args.append('--force')
    
    return cmd_args

def build_extract_args(args):
    """Use optimized argument builder"""
    from argument_builder import build_etm_extract_args
    return build_etm_extract_args(args)

def extract_max_excited_state_from_properties(properties):
    """Extract the maximum excited state number from property requests like s2o, s5e, etc."""
    max_state = 1  # Default minimum
    
    for prop in properties:
        # Look for patterns like s1e, s2o, s5_excitation_energy, etc.
        if prop.startswith('s') and ('e' in prop or 'o' in prop or '_' in prop):
            # Extract number from patterns like s2o, s5e, s1_excitation_energy
            if '_' in prop:
                # Handle s1_excitation_energy, s2_oscillator_strength format
                state_part = prop.split('_')[0]
                if state_part.startswith('s') and state_part[1:].isdigit():
                    state_num = int(state_part[1:])
                    max_state = max(max_state, state_num)
            else:
                # Handle s2o, s5e format
                for i, char in enumerate(prop[1:], 1):
                    if not char.isdigit():
                        if i > 1:  # We found at least one digit
                            state_num = int(prop[1:i])
                            max_state = max(max_state, state_num)
                        break
    
    return max_state

def extract_requested_excited_states(properties):
    """Extract which specific excited states the user requested for construction in step 4."""
    requested_states = set()
    
    for prop in properties:
        if prop.startswith('s') and ('e' in prop or 'o' in prop or '_' in prop):
            if '_' in prop:
                # Handle s1_excitation_energy, s2_oscillator_strength format
                state_part = prop.split('_')[0]
                if state_part.startswith('s') and state_part[1:].isdigit():
                    state_num = int(state_part[1:])
                    requested_states.add(state_num)
            else:
                # Handle s2o, s5e format
                for i, char in enumerate(prop[1:], 1):
                    if not char.isdigit():
                        if i > 1:  # We found at least one digit
                            state_num = int(prop[1:i])
                            requested_states.add(state_num)
                        break
    
    return sorted(list(requested_states))

def main():
    parser = argparse.ArgumentParser(
        description="ETM Suite - Complete workflow for electrostatic tuning map analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python main.py benzene.xyz --properties all         # Compute all properties for benzene
  python main.py water.xyz --properties homo lumo gap # Compute only HOMO, LUMO, and gap  
  python main.py H2O2.xyz --properties all --threads 8 # Use 8 total threads
  python main.py molecule.xyz --properties electronegativity --detached # Run in background
  python main.py benzene.xyz --properties all --keep-inputs # Preserve input files
        """
    )
    
    # Molecule specification (REQUIRED)
    parser.add_argument('molecule', help='REQUIRED: Input molecule .xyz file')
    # Workflow control
    parser.add_argument('--skip-step1', '-ss1', action='store_true', help='Skip input generation (use existing inputs)')
    parser.add_argument('--skip-step2', '-ss2', action='store_true', help='Skip calculations (use existing outputs)')
    parser.add_argument('--skip-step3', '-ss3', action='store_true', help='Skip property extraction (use existing CSV)')
    parser.add_argument('--skip-step4', '-ss4', action='store_true', help='Skip MOL2 generation')
    # Computational parameters
    parser.add_argument('--threads', '-ttt', '--total-threads', type=int, default=4,
                       help='Total number of threads to use (default: 4)')
    parser.add_argument('--threads-per-job', '-tpj', type=int, default=2,
                       help='Threads per individual job (default: 2)')
    parser.add_argument('--detached', '-det', action='store_true',
                       help='Run in detached mode - jobs run in background, terminal can be closed')
    # Molecular parameters (override JSON config)
    parser.add_argument('--vdw-scale', '-vds', type=float, help='VDW surface scaling factor (default from JSON)')
    parser.add_argument('--density', '-dns', type=float, help='Surface point density (default from JSON)')
    parser.add_argument('--basis', '-bss', help='Quantum chemistry basis set (default from JSON)')
    parser.add_argument('--method', '-mtd', help='QM method, e.g., B3LYP, M06-2X (default from JSON)')
    parser.add_argument('--charge', '-chg', type=int, help='Molecular charge (default from JSON)')
    parser.add_argument('--multiplicity', '-mult', type=int, help='Spin multiplicity (default from JSON)')
    parser.add_argument('--charge-anion', '-chga', type=int, help='Anion charge (default: charge-1)')
    parser.add_argument('--multiplicity-anion', '-multa', type=int, help='Anion multiplicity (auto-determined)')
    parser.add_argument('--charge-cation', '-chgc', type=int, help='Cation charge (default: charge+1)')
    parser.add_argument('--multiplicity-cation', '-multc', type=int, help='Cation multiplicity (auto-determined)')
    parser.add_argument('--solvent', '-sol', type=str, help='Solvent to use for calculation (optional)')
    # Local execution options
    parser.add_argument('--force', '-frs', action='store_true', help='Force rerun of completed jobs')
    # File management
    parser.add_argument('--keep-inputs', '-ki', action='store_true', 
                       help='Keep input files after completion')
    parser.add_argument('--clean-outputs', '-co', action='store_true', default=True,
                       help='Clean output files after extraction (default: True)')
    # Property selection (REQUIRED)
    # Accept both long and short forms for all property choices
    # Generate property choices including excited states up to s10
    base_properties = [
        'all', 'ground_state_energy', 'gse', 'homo', 'hom', 'lumo', 'lum', 'gap',
        'electron_affinity', 'eaf', 'ionization_energy', 'ioe', 'electronegativity', 'eng',
        'hardness', 'hrd', 'electrophilicity', 'efl', 'nucleophilicity', 'nfl',
        'chemical_potential', 'cpt', 'dipole_moment', 'dpm', 'excitation_energy', 'exe',
        'oscillator_strength', 'osc', 'redox', 'rdx'
    ]
    
    # Add excited state properties up to s10 (should be sufficient for most cases)
    excited_state_properties = []
    for i in range(1, 11):  # s1 through s10
        excited_state_properties.extend([
            f's{i}_excitation_energy', f's{i}e',
            f's{i}_oscillator_strength', f's{i}o'
        ])
    
    property_choices = base_properties + excited_state_properties
    parser.add_argument('--properties', '-prop', nargs='+', required=True,
                       choices=property_choices,
                       help='REQUIRED: Select properties to compute. Use "all" for all properties or specify individual ones. Short and long forms accepted. "excitation_energy" and "oscillator_strength" refer to the excited state specified by --excited-state.')


    # Add short forms for s1/s2 excited state properties
    # Excited state configuration
    parser.add_argument('--excited-state', '-exs', type=int, default=1,
                       help='Number of excited states to calculate (tdscf_states, default: 1)')
    # Output control
    parser.add_argument('--csv-path', '-csvp', default='properties.csv', help='Properties CSV output path (default: properties.csv)')
    parser.add_argument('--extract-csv', '-xcsv', default='etm_properties.csv', help='Final ETM properties CSV (default: etm_properties.csv)')
    # Debug options
    parser.add_argument('--single-index', '-si', help='Extract single charge index for debugging (e.g., "0000" or "Ref")')

    args = parser.parse_args()
    
    # Load local execution settings from JSON config and apply as defaults
    local_config = load_local_config()
    
    # Apply local settings from JSON as defaults
    if 'threads_per_job' in local_config and args.threads_per_job == 2:  # Use JSON default if not specified
        args.threads_per_job = local_config['threads_per_job']
    if 'total_threads' in local_config and args.threads == 4:  # Use JSON default if not specified
        args.threads = local_config['total_threads']
    
    # Show applied local settings
    if local_config:
        settings_str = ', '.join(f'{k}={v}' for k, v in local_config.items())
        print(f"{Fore.BLUE}📋 Applied local settings: {settings_str}{Style.RESET_ALL}")
    
    # Set up quiet mode
    if hasattr(args, 'quiet') and args.quiet:
        # Redirect stdout for subprocesses but keep our prints
        pass
    
    start_time = time.time()
    
    print_header("       Electrostatic Tuning Map Suite (ETM Suite)\n                  by Stephen O. Ajagbe")

    # Check prerequisites
    if not check_prerequisites():
        sys.exit(1)
    
    # Determine molecule name
    if args.molecule:
        # Handle both 'molecule' and 'molecule.xyz' formats
        if args.molecule.endswith('.xyz'):
            molecule_xyz = args.molecule
            molecule_name = os.path.splitext(args.molecule)[0]
        else:
            molecule_xyz = f"{args.molecule}.xyz"
            molecule_name = args.molecule
        print(f"{Fore.BLUE}📋 Using specified molecule: {molecule_name}{Style.RESET_ALL}")
    else:
        print_error("No molecule XYZ file specified. Please provide a molecule file.")
        print(f"{Fore.YELLOW}� Example: python main.py benzene.xyz --properties all{Style.RESET_ALL}")
        sys.exit(1)
    
    # Force local mode only
    args.mode = 'local'
    
    # Map short forms to long forms for all properties
    property_longform_map = {
        'gse': 'ground_state_energy',
        'hom': 'homo',
        'lum': 'lumo',
        'gap': 'gap',
        'eaf': 'electron_affinity',
        'ioe': 'ionization_energy',
        'eng': 'electronegativity',
        'hrd': 'hardness',
        'efl': 'electrophilicity',
        'nfl': 'nucleophilicity',
        'cpt': 'chemical_potential',
        'dpm': 'dipole_moment',
        'exe': 'excitation_energy',
        'osc': 'oscillator_strength',
        'rdx': 'redox',
    }
    
    # Add excited state mappings dynamically
    for i in range(1, 11):  # s1 through s10
        property_longform_map[f's{i}e'] = f's{i}_excitation_energy'
        property_longform_map[f's{i}o'] = f's{i}_oscillator_strength'
    # First, determine the maximum excited state needed from user requests
    max_excited_state_needed = extract_max_excited_state_from_properties(args.properties)
    
    # Override excited_state with the dynamically determined value
    args.excited_state = max(max_excited_state_needed, getattr(args, 'excited_state', 1))
    
    # Extract which specific states the user requested for step 4 construction
    requested_excited_states = extract_requested_excited_states(args.properties)
    args.requested_excited_states = requested_excited_states
    
    # Handle "all" properties option
    user_selected_all = 'all' in args.properties
    if user_selected_all:
        all_properties = [
            'ground_state_energy', 'homo', 'lumo', 'gap', 'electron_affinity',
            'ionization_energy', 'electronegativity', 'hardness',
            'electrophilicity', 'nucleophilicity', 'chemical_potential',
            'dipole_moment', 'redox'
        ]
        # Add all excitation energies and oscillator strengths up to max needed state
        for i in range(1, args.excited_state+1):
            all_properties.append(f's{i}_excitation_energy')
            all_properties.append(f's{i}_oscillator_strength')
        args.properties = all_properties
        print(f"{Fore.BLUE}📋 Computing all properties up to excited state {args.excited_state}{Style.RESET_ALL}")
        print(f"{Fore.BLUE}📋 Will construct MOL2 files for states: {requested_excited_states if requested_excited_states else [1]}{Style.RESET_ALL}")
    else:
        mapped_properties = []
        for p in args.properties:
            longform = property_longform_map.get(p, p)
            # If user requests excitation_energy or oscillator_strength, expand for all states up to max needed
            if longform == 'excitation_energy':
                for i in range(1, args.excited_state+1):
                    mapped_properties.append(f's{i}_excitation_energy')
            elif longform == 'oscillator_strength':
                for i in range(1, args.excited_state+1):
                    mapped_properties.append(f's{i}_oscillator_strength')
            else:
                mapped_properties.append(longform)
        args.properties = mapped_properties
        print(f"{Fore.BLUE}📋 Will calculate up to excited state {args.excited_state} (auto-detected from: {', '.join([p for p in args.properties if 's' in p])}){Style.RESET_ALL}")
        print(f"{Fore.BLUE}📋 Will construct MOL2 files for states: {requested_excited_states if requested_excited_states else [1]}{Style.RESET_ALL}")
        print(f"{Fore.BLUE}📋 Selected properties: {', '.join(args.properties)}{Style.RESET_ALL}")
    
    success = True
    # Step 1: Generate inputs
    if not args.skip_step1:
        # Define folder requirements for each property
        property_folders = {
            "redox": "anion_gfec, anion_sscf, neutral_gfec, neutral_sscf",
            "ground_state_energy": "neutral",
            "homo": "neutral",
            "lumo": "neutral",
            "gap": "neutral",
            "electron_affinity": "anion, neutral",
            "ionization_energy": "cation, neutral",
            "electronegativity": "anion, cation, neutral",
            "hardness": "anion, cation, neutral",
            "electrophilicity": "anion, cation, neutral",
            "nucleophilicity": "anion, cation, neutral",
            "chemical_potential": "anion, cation, neutral",
            "dipole_moment": "neutral",
        }
        # Add excited state properties
        for i in range(1, args.excited_state + 1):
            property_folders[f"s{i}_excitation_energy"] = "neutral"
            property_folders[f"s{i}_oscillator_strength"] = "neutral"

        # Build the properties block for JSON
        def build_properties_block(props):
            return {p: property_folders.get(p, "neutral") for p in props}

        # The properties to use are already mapped in args.properties
        properties_block = build_properties_block(args.properties)

        # Write the JSON file (or pass to your input generation function)
        psi4_params = {
            "vdw_scale": getattr(args, "vdw_scale", 1.0),
            "density": getattr(args, "density", 1.0),
            "basis": getattr(args, "basis", "6-31G*"),
            "method": getattr(args, "method", "b3lyp"),
            "charge": getattr(args, "charge", 0),
            "multiplicity": getattr(args, "multiplicity", 1),
            "properties": properties_block,
            "xyz": molecule_xyz,
            "local": {
                "threads_per_job": args.threads_per_job,
                "total_threads": args.threads
            }
        }
        with open('psi4_params.json', 'w') as f:
            json.dump(psi4_params, f, indent=2)
        print(f"{Fore.GREEN}✅ psi4_params.json written with properties: {list(properties_block.keys())}{Style.RESET_ALL}")

        # Now generate inputs as usual
        success = step1_generate_inputs(molecule_xyz, args)
        if not success:
            print_error("Step 1 failed")
            sys.exit(1)
    else:
        print_step("1", "Generate Psi4 Input Files (SKIPPED)")
    
    # Step 2: Submit jobs
    if not args.skip_step2:
        if args.detached:
            # Detached mode - submit jobs in background
            success = step2_submit_jobs_detached(molecule_name, args)
            if not success:
                print_error("Step 2 (detached) failed")
                sys.exit(1)
            
            # In detached mode, we exit here - continuation will happen automatically
            print(f"\n{Fore.GREEN}🎉 ETM Suite started in detached mode!{Style.RESET_ALL}")
            print(f"{Fore.YELLOW}📋 Jobs are running in the background{Style.RESET_ALL}")
            print(f"{Fore.YELLOW}📋 Steps 3 and 4 will run automatically after jobs complete{Style.RESET_ALL}")
            print(f"{Fore.YELLOW}📋 You can safely close this terminal{Style.RESET_ALL}")
            return 0
        else:
            # Regular mode - wait for jobs to complete
            success = step2_submit_jobs(
                molecule_name, 
                mode=args.mode,
                threads_per_job=args.threads_per_job,
                total_threads=args.threads,
                clean_outputs=args.clean_outputs,
                args=args
            )
            if not success:
                print_error("Step 2 failed")
                sys.exit(1)
    else:
        print_step("2", "Submit and Run Quantum Chemistry Calculations (SKIPPED)")
    
    # In detached mode, these steps are handled by the continuation script
    # In regular mode, we continue with steps 3 and 4
    
    # Step 3: Extract properties
    if not args.skip_step3:
        success = step3_extract_properties(args, molecule_name)
        if not success:
            print_error("Step 3 failed")
            sys.exit(1)
    else:
        print_step("3", "Extract Molecular Properties (SKIPPED)")
    
    # Step 4: Generate MOL2 files
    if not args.skip_step4:
        success = step4_generate_mol2(keep_inputs=args.keep_inputs, molecule_name=molecule_name, args=args)
        if not success:
            print_error("Step 4 failed")
            sys.exit(1)
    else:
        print_step("4", "Generate MOL2 Visualization Files (SKIPPED)")
    
    # Print summary
    if success:
        print_summary(molecule_name, start_time, keep_inputs=args.keep_inputs)
    
    return 0

if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt:
        print(f"\n{Fore.YELLOW}⚠️  Workflow interrupted by user{Style.RESET_ALL}")
        sys.exit(1)
    except Exception as e:
        print_error(f"Unexpected error: {e}")
        sys.exit(1)
