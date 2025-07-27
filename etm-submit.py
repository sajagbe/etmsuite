import os
import argparse
import csv
import shutil
import subprocess
import tempfile
import json
import time
import re
import glob
from concurrent.futures import ThreadPoolExecutor, as_completed
from cclib import io as cclibio

STATES = ['neutral', 'cation', 'anion']
EXT_KEEP = {'.in', '.out'}
HARTREE_TO_EV = 27.2114
HARTREE_TO_KCAL = 627.509
EV_CONV = 0.0367493

def cleanup_psi4_files(working_dir='.'):
    """
    Clean up Psi4 cleanup manifest files (psi.*.clean) left behind after calculations.
    These files are created by Psi4 to track temporary files but are not automatically removed.
    """
    clean_files = glob.glob(os.path.join(working_dir, 'psi.*.clean'))
    if clean_files:
        print(f"[CLEANUP] Removing {len(clean_files)} Psi4 cleanup manifest files...")
        for clean_file in clean_files:
            try:
                os.remove(clean_file)
            except OSError as e:
                print(f"[WARN] Could not remove {clean_file}: {e}")
        print(f"[CLEANUP] Psi4 cleanup files removed")
    else:
        print(f"[CLEANUP] No Psi4 cleanup files found to remove")

JOB_MAP = {
    'neutral': {
        'props': [
            'ground_state_energy', 'homo', 'lumo', 'gap',
            'dipole_moment', 's1_excitation_energy', 's1_oscillator_strength'
        ]
    },
    'anion': {
        'props': ['electron_affinity']
    },
    'cation': {
        'props': ['ionization_energy']
    }
}

def is_job_completed(input_file):
    """
    Check if a job has already completed successfully by looking for existing output files.
    Returns True if a completed output file exists, False otherwise.
    """
    base_name = input_file.replace('.in', '')
    
    # Check for output files with attempt numbers
    for attempt in range(1, 4):  # Check attempts 1-3
        output_file = f"{base_name}_attempt{attempt}.out"
        if os.path.exists(output_file):
            try:
                with open(output_file, 'r') as f:
                    content = f.read()
                if 'Psi4 exiting successfully' in content:
                    print(f"[SKIP] {input_file} already completed (found {output_file})")
                    return True
            except Exception as e:
                print(f"[WARN] Could not read {output_file}: {e}")
                continue
    
    # Only check for .out files from local jobs
    array_output = input_file.replace('.in', '.out')
    if os.path.exists(array_output):
        try:
            with open(array_output, 'r') as f:
                content = f.read()
            if 'Psi4 exiting successfully' in content:
                print(f"[SKIP] {input_file} already completed (found {array_output})")
                return True
        except Exception as e:
            print(f"[WARN] Could not read {array_output}: {e}")
    
    return False

def collect_in_files(root_folder, target_states, force_rerun=False):
    """Collect input files and filter out already completed jobs"""
    files = []
    
    # Check if this is a redox calculation directory (ends with _redox_inputs)
    if root_folder.endswith('_redox_inputs'):
        # New flat redox folder structure: molecule_redox_inputs/{anion_gfec, anion_sscf, neutral_gfec, neutral_sscf}/
        for subdir in os.listdir(root_folder):
            subdir_path = os.path.join(root_folder, subdir)
            if os.path.isdir(subdir_path):
                files.extend(
                    sorted(os.path.join(subdir_path, f) for f in os.listdir(subdir_path) if f.endswith('.in'))
                )
    else:
        # Handle standard directory structure
        for state in target_states:
            subdir = os.path.join(root_folder, state)
            if os.path.isdir(subdir):
                files.extend(
                    sorted(os.path.join(subdir, f) for f in os.listdir(subdir) if f.endswith('.in'))
                )
    
    if force_rerun:
        print(f"[INFO] Force mode enabled - will rerun all {len(files)} jobs")
        return files
    
    # Filter out completed jobs
    uncompleted_files = []
    completed_count = 0
    
    for file in files:
        if is_job_completed(file):
            completed_count += 1
        else:
            uncompleted_files.append(file)
    
    if completed_count > 0:
        print(f"[INFO] Found {completed_count} already completed jobs, skipping them")
        print(f"[INFO] {len(uncompleted_files)} jobs need to be run")
    
    return uncompleted_files

def run_single_psi4(file, threads, attempt=1):
    with tempfile.TemporaryDirectory() as tmp_dir:
        base = os.path.basename(file).replace('.in', '')
        tmp_input = os.path.join(tmp_dir, f'{base}.in')
        tmp_output = os.path.join(tmp_dir, f'{base}.out')
        shutil.copy(file, tmp_input)

        result = subprocess.run(
            ['psi4', tmp_input], cwd=tmp_dir,
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )

        # Check if Psi4 output file was created
        if not os.path.exists(tmp_output):
            print(f"[PSI4 ERROR] Output file missing: {tmp_output} (attempt {attempt})")
            print(f"Return code: {result.returncode}")
            print(f"Stdout + Stderr:\n{result.stdout}\n{result.stderr}")
            return False

        with open(tmp_output, 'r') as f:
            out_content = f.read()

        if 'Psi4 exiting successfully' not in out_content:
            print(f"[PSI4 ERROR] {file} (attempt {attempt}) - success string not found in output file")
            tail_lines = out_content.strip().split('\n')[-20:]
            print("\n".join(tail_lines))
            return False

        # Save output file back with attempt number
        dest_out = file.replace('.in', f'_attempt{attempt}.out')
        shutil.copy(tmp_output, dest_out)
        print(f"[SUCCESS] {file} (attempt {attempt})")
        return True

def run_local_with_retries_parallel(in_files, threads_per_job, max_workers, log_path):
    def wrapper(infile):
        return None if any(run_single_psi4(infile, threads_per_job, a) for a in range(1, 3)) else infile

    failed = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(wrapper, f) for f in in_files]
        for future in as_completed(futures):
            result = future.result()
            if result: failed.append(result)

    # Clean up Psi4 files after local execution
    cleanup_psi4_files()

    if failed:
        with open(log_path, 'w', newline='') as f:
            csv.writer(f).writerows([['Filename']] + [[j] for j in failed])
        print(f"\n[LOGGED] {len(failed)} failed jobs in {log_path}")
        return False
    else:
        print("\n[INFO] All jobs completed successfully.")
        return True

# Main execution
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Submit Psi4 jobs locally')
    parser.add_argument('--input', '-inp', required=True,
                       help='Input directory containing .in files')
    parser.add_argument('--threads_per_job', '-tpj', type=int, default=2,
                       help='Number of threads per Psi4 job')
    parser.add_argument('--total_threads', '--threads', '-ttt', type=int, default=4,
                       help='Total number of threads to use')
    parser.add_argument('--clean_outputs', '-co', action='store_true',
                       help='Clean up output files after successful completion')
    parser.add_argument('--force', '-frs', action='store_true',
                       help='Force rerun of completed jobs')
    parser.add_argument('--excited-state', '-exs', type=int, default=1,
                       help='Number of excited states (for consistency with other scripts)')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.input):
        print(f"[ERROR] Input directory does not exist: {args.input}")
        exit(1)
    
    # Collect input files
    in_files = collect_in_files(args.input, STATES, force_rerun=args.force)
    
    if not in_files:
        print("[INFO] No input files to process")
        cleanup_psi4_files()  # Clean up any leftover files
        exit(0)
    
    print(f"[INFO] Found {len(in_files)} input files to process")
    
    # Execute based on mode
    max_workers = max(1, args.total_threads // args.threads_per_job)
    print(f"[LOCAL] Using {max_workers} parallel workers with {args.threads_per_job} threads each")
    
    success = run_local_with_retries_parallel(
        in_files, args.threads_per_job, max_workers, 'failed_jobs.csv'
    )
    
    # Clean up after local execution
    cleanup_psi4_files()
    
    # Final status
    if success:
        print("[FINAL] All jobs completed successfully!")
        if args.clean_outputs:
            print("[CLEANUP] Cleaning up output files...")
            # Additional cleanup logic if needed
    else:
        print("[FINAL] Some jobs failed after all retry attempts")
        exit(1)
