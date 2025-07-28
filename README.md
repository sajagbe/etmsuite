# ETM Suite (Electrostatic Tuning Map Suite)

**A comprehensive quantum chemistry workflow for computing positional electrostatic effects on molecular properties.**

## 🎯 Overview

The ETM Suite is a complete workflow that generates electrostatic tuning maps by:
1. Creating van der Waals surfaces around molecules
2. Placing point charges at surface positions
3. Computing molecular properties with external electrostatic perturbations
4. Extracting property variations across the surface
5. Generating MOL2 visualization files for molecular viewers

## 🚀 Quick Start

1. **Clone GitHub**
   ```bash
   git clone https://github.com/sajagbe/etmsuite.git

   cd etmsuite

   ```

2. **Create the Conda environment:**

   ```bash
   conda env create -f etmsuite.yml
   ```

3. **Activate the environment:**

   ```bash
   conda activate etmsuite
   ```

4. **Run the workflow using SLURM:**

   - Edit `etmsuite.sh` and replace `molecule.xyz` with your input file (e.g., `water.xyz`).
   - Add calculation options as needed using descriptions below. 
   - Submit the job:

     ```bash
     sbatch etmsuite.sh
     ```
   

   OR:

4. **Quick Run**
```bash
# Compute all properties for water
python main.py water.xyz -prop all 

# Compute specific properties with solvent
python main.py benzene.xyz -prop homo lumo gap -det

# Compute redox and excited state properties (up to S2)
python main.py molecule.xyz -prop redox exe osc -exs 2
```


## 📋 Molecular Properties Available

### Standard Electronic Properties (12)
| Short | Long | Description |
|-------|------|-------------|
| `gse` | ground_state_energy | Ground state SCF energy |
| `hom` | homo | Highest occupied molecular orbital energy |
| `lum` | lumo | Lowest unoccupied molecular orbital energy |
| `gap` | gap | HOMO-LUMO gap |
| `dpm` | dipole_moment | Electric dipole moment |
| `eaf` | electron_affinity | Electron affinity (neutral + anion) |
| `ioe` | ionization_energy | Ionization energy (neutral + cation) |
| `cpt` | chemical_potential | Electronic chemical potential |
| `eng` | electronegativity | Mulliken electronegativity |
| `hrd` | hardness | Chemical hardness |
| `efl` | electrophilicity | Electrophilic character |
| `nfl` | nucleophilicity | Nucleophilic character |

### Excited State Properties (2)
| Short | Long | Description |
|-------|------|-------------|
| `exe` | excitation_energy | Electronic excitation energy (TDDFT, expands to all states up to -exs) |
| `osc` | oscillator_strength | Oscillator strength (TDDFT, expands to all states up to -exs) |

### Redox Properties (1 requestable property → 4 computed values)
| Short | Long | Description |
|-------|------|-------------|
| `rdx` | redox | Redox potential (outputs: reduction_gfe, delta_reduction_gfe, redox_potential, delta_redox_potential) |

## ⚗️ Solvent Models

### Non-Redox Calculations
```bash
# Gas phase (default)
python main.py molecule.xyz -prop homo

# Water solvent (default when -sol specified)
python main.py molecule.xyz -prop homo -sol Water

# Custom solvent
python main.py molecule.xyz -prop homo -sol Acetonitrile
python main.py molecule.xyz -prop homo -sol Benzene
python main.py molecule.xyz -prop homo -sol DMSO
```

### Redox Calculations
```bash
# Default: Acetonitrile for solvated steps, gas phase for thermal corrections
python main.py molecule.xyz -prop redox

# Custom redox solvent
python main.py molecule.xyz -prop redox -rsol DMSO -rsolt IEFPCM -rsolr BONDI

```

### Solvent Options
**Note:** Advanced solvent options are configured through the JSON configuration file, not command-line flags.

| Parameter | Options | Default |
|-----------|---------|---------|
| solver_type | CPCM, IEFPCM | CPCM |
| radii_set | UFF, BONDI | UFF |
| cavity_type | GePol | GePol |

## 🎛️ Computational Parameters

### Surface Generation
```bash
-vds 1.2    # VDW surface scaling factor (default: 1.0)
-dns 0.5    # Surface point density (default: 1.0)
```

### Quantum Chemistry
```bash
-bss cc-pVDZ      # Basis set (default: 6-31G*)
-mtd m06-2x       # QM method (default: b3lyp)
-chg 1            # Molecular charge (default: 0)
-mult 2           # Spin multiplicity (default: 1)
```

### Charge States
```bash
-chga -2       # Anion charge (default: charge-1)
-multa 1       # Anion multiplicity (auto-determined)
-chgc 2        # Cation charge (default: charge+1)
-multc 1       # Cation multiplicity (auto-determined)
```

### Excited States
```bash
-exs 3            # Number of excited states to calculate (default: 1)
```

## 🎮 Workflow Control

### Step Control
```bash
-ss1    # Skip input generation (use existing)
-ss2    # Skip calculations (use existing outputs)
-ss3    # Skip property extraction (use existing CSV)
-ss4    # Skip MOL2 generation
```

### Execution Control
```bash
-ttt 8     # Total threads to use (default: 4)
-tpj 4     # Threads per individual job (default: 2)
-det       # Run in background
-frs       # Force rerun completed jobs
```

### File Management
```bash
-ki         # Preserve input files after completion
-co         # Clean output files after extraction (default: True)
-csvp path  # Custom CSV output path (default: properties.csv)
-xcsv path  # Final ETM CSV path (default: etm_properties.csv)
```

### Debug and Analysis Options
```bash
-si 0000   # Extract single surface point for debugging
-si Ref    # Extract only reference calculation
```

## 📊 Output Files

### Generated Files
```
molecule_vdw_surface.txt              # VDW surface coordinates
molecule_psi4_inputs/                 # Standard QM input files
molecule_redox_inputs/                # Redox QM input files (if applicable)
molecule_etm_properties.csv           # Extracted properties with deltas
molecule_etm_results/                 # MOL2 visualization files
  ├── property1_tuning_map.mol2
  ├── property2_tuning_map.mol2
  └── ...
```

### CSV Output Format
```csv
ref,potential_x,potential_y,potential_z,property,delta_property
0000,1.234,5.678,9.012,0.123,-0.045
0001,2.345,6.789,0.123,0.156,-0.012
Ref,0.000,0.000,0.000,0.168,0.000
...
```

### Advanced Configuration Features

#### JSON Configuration File (`psi4_params.json`)
The ETM Suite uses a JSON configuration file for default settings:

```json
{
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
  "xyz": "water.xyz",
  "local": {
    "threads_per_job": 2,
    "total_threads": 4
  }
}
```

#### Auto-Detection Features
- **Directory Detection**: Auto-detects `molecule_psi4_inputs/` and `molecule_redox_inputs/` folders
- **Config Creation**: Creates default `psi4_params.json` if none exists.
- **Surface File Detection**: Finds `molecule_vdw_surface.txt` automatically
- **Individual Scripts**: `etm-inputs.py` and `etm-extract.py` can auto-detect molecule from XYZ files when run directly

#### Redox Calculation Architecture
Redox calculations use a sophisticated 4-step thermodynamic cycle:

1. **Directory Structure**:
   ```
   molecule_redox_inputs/
     ├── neutral_gfec/    # Neutral Gibbs free energy correction calculations
     ├── neutral_sscf/    # Neutral solvated SCF calculations
     ├── anion_gfec/      # Anion Gibbs free energy correction calculations
     └── anion_sscf/      # Anion solvated SCF calculations
   ```

2. **Calculation Types**:
   - **GFEC**: Gibbs free energy correction (gas-phase frequency for thermal corrections)
   - **SSCF**: Solvated single-point for electronic energies (with solvent)

3. **Property Calculation**:
   - `reduction_gfe = (E_anion_sscf + G_anion_gfec) - (E_neutral_sscf + G_neutral_gfec)`
   - `redox_potential = -reduction_gfe / F - 4.43` (vs NHE in Volts)

### Failure Handling System

#### Automatic Retry Mechanism
- Each job gets **2 automatic retry attempts**
- Failed jobs logged to `failed_jobs.csv`
- Detailed error output for debugging

#### Graceful Degradation
```bash
# If >10 successful outputs: Continue with partial data
# If <10 successful outputs: Stop with error

"Found 195 successful outputs - continuing with extraction..."
# vs
"ERROR: Too few successful calculations (3 outputs)"
```

#### Missing Data Handling
- Individual surface points with missing data are skipped
- Properties with missing dependencies show as skipped
- Reference calculations still work with partial datasets

### Output Validation
1. **File Structure**: Verify all expected directories and files are created
2. **CSV Format**: Check column names, data types, and delta calculations
3. **MOL2 Format**: Validate molecular viewer compatibility  
4. **Property Values**: Ensure reasonable physical values and units
5. **Delta Calculations**: Verify reference subtraction is correct
6. **Redox Validation**: Check 4-step thermodynamic cycle completeness

## ⚠️ Important Notes

1. **Excited States**: The `-exs N` parameter sets how many excited states to calculate. Properties `exe` and `osc` automatically expand to generate properties for all states from S1 to SN.

2. **Solvent Independence**: `-sol` controls non-redox calculations. Redox calculations use default solvent settings from the JSON configuration file.

3. **Delta Properties**: All properties automatically get corresponding `delta_property` columns showing variation from the reference calculation.

4. **Redox Requirements**: Redox calculations require both gas-phase frequency and solvated SCF calculations for thermodynamic accuracy.

5. **File Cleanup**: By default, input files are cleaned after completion unless `-ki` is specified.

## 📚 Example Workflows

### Basic Property Mapping
```bash
python main.py benzene.xyz -prop homo lumo gap
```

### Comprehensive Analysis with Solvent
```bash
python main.py molecule.xyz -prop all -sol Water -exs 2
```

### Redox Potential Analysis
```bash
python main.py chromophore.xyz -prop redox
```

### High-Performance Calculation
```bash
python main.py large_molecule.xyz -prop all -ttt 16 -tpj 4 -det
```

## 🏷️ All Command-Line Flags & Options

Below is a complete list of all flags, options, and job types for the ETM Suite using short-form syntax:

### Core Arguments
* `-vds` : van der Waals scale
* `-dns` : surface density  
* `-bss` : basis set
* `-mtd` : QM method
* `-chg` : molecular charge
* `-mult` : spin multiplicity

### Ion States
* `-chga` : anion charge
* `-multa` : anion multiplicity
* `-chgc` : cation charge
* `-multc` : cation multiplicity

### Solvent
* `-sol` : solvent name (e.g., Water, Acetonitrile, DMSO)

### Redox Solvent
* `-rsol` : redox solvent name (e.g., Acetonitrile, DMSO, Water)
* `-rsolt` : redox solvent type (CPCM, IEFPCM)
* `-rsolr` : redox solvent radii set (UFF, BONDI)

### Excited States
* `-exs` : number of excited states

### Properties (space-separated)
* `-prop` : property list (e.g. `-prop homo lumo gap redox`)
  * `ground_state_energy` or `gse` : ground-state-energy
  * `homo` or `hom` : homo
  * `lumo` or `lum` : lumo
  * `gap` : gap
  * `dipole_moment` or `dpm` : dipole-moment
  * `electron_affinity` or `eaf` : electron-affinity
  * `ionization_energy` or `ioe` : ionization-energy
  * `chemical_potential` or `cpt` : chemical-potential
  * `electronegativity` or `eng` : electronegativity
  * `hardness` or `hrd` : hardness
  * `electrophilicity` or `efl` : electrophilicity
  * `nucleophilicity` or `nfl` : nucleophilicity
  * `excitation_energy` or `exe` : excitation-energy (expands to all states up to -exs)
  * `oscillator_strength` or `osc` : oscillator-strength (expands to all states up to -exs)
  * `redox` or `rdx` : redox

### Workflow Control
* `-ss1` : skip step 1 (inputs)
* `-ss2` : skip step 2 (jobs)
* `-ss3` : skip step 3 (extraction)
* `-ss4` : skip step 4 (MOL2)
* `-frs` : force rerun
* `-det` : detached/background

### Threading & Performance
* `-ttt` : total threads
* `-tpj` : threads per job

### File Management
* `-ki` : keep inputs
* `-co` : clean outputs
* `-csvp` : properties CSV path
* `-xcsv` : ETM CSV path

### Extraction/Debug
* `-si` : single index (e.g. `-si 0000` or `-si Ref`)

---

