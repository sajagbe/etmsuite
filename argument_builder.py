#!/usr/bin/env python3
"""
ETM Suite Argument Builder
=========================
Provides high-performance, maintainable argument construction for all ETM Suite workflow components.
Supports dual long/short argument system and proper None handling.
"""

from typing import Any, List

class ArgumentBuilder:
    def __init__(self):
        self._args: List[str] = []

    def add_if_exists(self, args_obj: Any, attr_name: str, cli_flag: str, transform=None, condition=None) -> 'ArgumentBuilder':
        if hasattr(args_obj, attr_name):
            value = getattr(args_obj, attr_name)
            if condition is None or condition(value):
                if value is not None:
                    if transform:
                        value = transform(value)
                    self._args.extend([cli_flag, value])
        return self

    def add_flag_if_true(self, args_obj: Any, attr_name: str, cli_flag: str) -> 'ArgumentBuilder':
        if hasattr(args_obj, attr_name) and getattr(args_obj, attr_name):
            self._args.append(cli_flag)
        return self

    def add_list_if_exists(self, args_obj: Any, attr_name: str, cli_flag: str) -> 'ArgumentBuilder':
        if hasattr(args_obj, attr_name):
            value = getattr(args_obj, attr_name)
            if value:
                self._args.append(cli_flag)
                if isinstance(value, list):
                    self._args.extend([str(v) for v in value])
                else:
                    self._args.append(str(value))
        return self

    def build(self) -> List[str]:
        return self._args

# Inputs builder
def build_etm_inputs_args(args):
    builder = ArgumentBuilder()
    return (builder
        .add_if_exists(args, 'vdw_scale', '--vdw-scale', str, lambda x: x is not None)
        .add_if_exists(args, 'density', '--density', str, lambda x: x is not None)
        .add_if_exists(args, 'basis', '--basis', condition=lambda x: x is not None)
        .add_if_exists(args, 'method', '--method', condition=lambda x: x is not None)
        .add_if_exists(args, 'charge', '--charge', str, lambda x: x is not None)
        .add_if_exists(args, 'multiplicity', '--multiplicity', str, lambda x: x is not None)
        .add_if_exists(args, 'charge_anion', '--charge-anion', str, lambda x: x is not None)
        .add_if_exists(args, 'multiplicity_anion', '--multiplicity-anion', str, lambda x: x is not None)
        .add_if_exists(args, 'charge_cation', '--charge-cation', str, lambda x: x is not None)
        .add_if_exists(args, 'multiplicity_cation', '--multiplicity-cation', str, lambda x: x is not None)
        .add_if_exists(args, 'excited_state', '--excited-state', str, lambda x: x is not None)
        .add_list_if_exists(args, 'properties', '--properties')
        .add_if_exists(args, 'solvent', '--solvent', str, condition=lambda x: x is not None)
        .add_if_exists(args, 'solvent_type', '--solvent-type', str, condition=lambda x: hasattr(args, 'solvent') and args.solvent is not None)
        .add_if_exists(args, 'solvent_radii', '--solvent-radii', str, condition=lambda x: hasattr(args, 'solvent') and args.solvent is not None)
        .add_if_exists(args, 'redox_solvent', '--redox-solvent', str, condition=lambda x: x is not None)
        .add_if_exists(args, 'redox_solvent_type', '--redox-solvent-type', str, condition=lambda x: hasattr(args, 'redox_solvent') and args.redox_solvent is not None)
        .add_if_exists(args, 'redox_solvent_radii', '--redox-solvent-radii', str, condition=lambda x: hasattr(args, 'redox_solvent') and args.redox_solvent is not None)
        .build())

# Submit builder
def build_etm_submit_args(args):
    builder = ArgumentBuilder()
    return (builder
        .add_if_exists(args, 'input', '--input', str, lambda x: x is not None)
        .add_if_exists(args, 'threads_per_job', '--threads-per-job', str, lambda x: x is not None)
        .add_if_exists(args, 'threads', '--total-threads', str, lambda x: x is not None)
        .add_flag_if_true(args, 'clean_outputs', '--clean-outputs')
        .add_flag_if_true(args, 'force', '--force')
        .add_if_exists(args, 'excited_state', '--excited-state', str, lambda x: x is not None)
        .build())

# Extract builder
def build_etm_extract_args(args):
    builder = ArgumentBuilder()
    return (builder
        .add_if_exists(args, 'index', '--index', str, lambda x: x is not None)
        .add_if_exists(args, 'outfolder', '--outfolder', str, lambda x: x is not None)
        .add_if_exists(args, 'config', '--config', str, lambda x: x is not None)
        .add_if_exists(args, 'outfile', '--outfile', str, lambda x: x is not None)
        .add_if_exists(args, 'surface_file', '--surface-file', str, lambda x: x is not None)
        .add_if_exists(args, 'excited_state', '--excited-state', str, lambda x: x is not None)
        .add_list_if_exists(args, 'properties', '--requested-properties')
        .build())

# Construct builder
def build_etm_construct_args(args):
    builder = ArgumentBuilder()
    return (builder
        .add_flag_if_true(args, 'keep_inputs', '--keep-inputs')
        .add_if_exists(args, 'excited_state', '--excited-state', str, lambda x: x is not None)
        .build())
