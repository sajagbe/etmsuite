# ETM Suite Argument Mapping

## Comprehensive Command-Line Argument Reference

This document provides the complete mapping of long-form and short-form command-line arguments for the ETM Suite quantum chemistry workflow.

### Workflow Control

| Long Form | Short Form | Description |
|-----------|------------|-------------|
| `--skip-step1` | `-ss1` | Skip step 1 (input generation) |
| `--skip-step2` | `-ss2` | Skip step 2 (job submission) |
| `--skip-step3` | `-ss3` | Skip step 3 (property extraction) |
| `--skip-step4` | `-ss4` | Skip step 4 (MOL2 generation) |
| `--threads` | `-ttt` | Total number of threads (alias: --total-threads) |
| `--threads-per-job` | `-tpj` | Threads per individual job |
| `--detached` | `-det` | Run in detached mode |

### Core Parameters

| Long Form | Short Form | Description |
|-----------|------------|-------------|
| `--vdw-scale` | `-vds` | Van der Waals scaling factor |
| `--density` | `-dns` | Electron density value |
| `--basis` | `-bss` | Basis set |
| `--method` | `-mtd` | Quantum chemistry method |
| `--charge` | `-chg` | Molecular charge |
| `--multiplicity` | `-mult` | Spin multiplicity |

### Ion States

| Long Form | Short Form | Description |
|-----------|------------|-------------|
| `--charge-anion` | `-chga` | Anion charge state |
| `--multiplicity-anion` | `-multa` | Anion multiplicity |
| `--charge-cation` | `-chgc` | Cation charge state |
| `--multiplicity-cation` | `-multc` | Cation multiplicity |

### Solvent Configuration

| Long Form | Short Form | Description |
|-----------|------------|-------------|
| `--solvent` | `-sol` | Solvent name |
| `--solvent-type` | `-solt` | Solvent model type |
| `--solvent-radii` | `-solr` | Solvent radii set |
| `--redox-solvent` | `-rsol` | Redox calculation solvent |
| `--redox-solvent-type` | `-rsolt` | Redox solvent model type |
| `--redox-solvent-radii` | `-rsolr` | Redox solvent radii set |

### File Management

| Long Form | Short Form | Description |
|-----------|------------|-------------|
| `--force` | `-frs` | Force overwrite existing files |
| `--keep-inputs` | `-ki` | Keep input files after processing |
| `--clean-outputs` | `-co` | Clean output files |
| `--single-index` | `-si` | Process single charge index |
| `--csv-path` | `-csvp` | Properties CSV output path (default: properties.csv) |
| `--extract-csv` | `-xcsv` | Final ETM properties CSV (default: etm_properties.csv) |

### Excited States

| Long Form | Short Form | Description |
|-----------|------------|-------------|
| `--excited-state` | `-exs` | Excited state number |

### Properties

| Long Form | Short Form | Description |
|-----------|------------|-------------|
| `--properties` | `-prop` | List of properties to calculate |

### ETM-Submit Specific

| Long Form | Short Form | Description |
|-----------|------------|-------------|
| `--input` | `-inp` | Input directory path |
| `--total-threads` | `-ttt` | Total threads for submission |

### ETM-Extract Specific

| Long Form | Short Form | Description |
|-----------|------------|-------------|
| `--index` | `-idx` | Charge index to extract |
| `--outfolder` | `-ofd` | Output folder path |
| `--config` | `-cfg` | Configuration file |
| `--outfile` | `-ofl` | Output file name |
| `--surface-file` | `-sfl` | Surface file path |

## Property Values (for --properties flag)

### Electronic Properties

| Long Form | Short Form | Description |
|-----------|------------|-------------|
| `ground-state-energy` | `gse` | Ground state energy |
| `homo` | `hom` | Highest occupied molecular orbital |
| `lumo` | `lum` | Lowest unoccupied molecular orbital |
| `gap` | `gap` | HOMO-LUMO gap |

### Redox Properties

| Long Form | Short Form | Description |
|-----------|------------|-------------|
| `electron-affinity` | `eaf` | Electron affinity |
| `ionization-energy` | `ioe` | Ionization energy |
| `redox` | `rdx` | Redox potential |

### Chemical Descriptors

| Long Form | Short Form | Description |
|-----------|------------|-------------|
| `electronegativity` | `eng` | Electronegativity |
| `hardness` | `hrd` | Chemical hardness |
| `electrophilicity` | `efl` | Electrophilicity index |
| `nucleophilicity` | `nfl` | Nucleophilicity index |
| `chemical-potential` | `cpt` | Chemical potential |

### Molecular Properties

| Long Form | Short Form | Description |
|-----------|------------|-------------|
| `dipole-moment` | `dpm` | Dipole moment |
| `excitation-energy` | `exe` | Excitation energy |
| `oscillator-strength` | `osc` | Oscillator strength |

## Solvent Property Values

| Long Form | Short Form | Description |
|-----------|------------|-------------|
| `solver-type` | `slt` | PCM solver type |
| `radii-set` | `rds` | Atomic radii set |
| `cavity-type` | `cvt` | Cavity type |

## Job Type Values

| Long Form | Short Form | Description |
|-----------|------------|-------------|
| `neutral` | `neu` | Neutral state calculation |
| `anion` | `ani` | Anion state calculation |
| `cation` | `cat` | Cation state calculation |
| `neutral_gfec` | `neu-gfc` | Neutral gibbs free energy correction|
| `neutral_sscf` | `neu-ssc` | Neutral solvated SCF |
| `anion_gfec` | `ani-gfc` | Anion gibbs free energy correction|
| `anion_sscf` | `ani-ssc` | Anion solvated SCF |

---

**Total Arguments**: 62 command-line arguments + property/value aliases

**Usage Examples**:
- `python main.py lumiflavin.xyz -bss cc-pVDZ -mtd B3LYP -prop homo,lumo,gap`
- `python etm-inputs.py lumiflavin.xyz -vds 1.2 -dns 0.001 -exs 1`
- `python etm-submit.py -inp inputs/ -tpj 4 -ttt 16 -frs`

**Note**: All short forms use single dash prefix (`-`) while long forms use double dash prefix (`--`).
