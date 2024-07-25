#!/usr/bin/env python
'''
Module to run CatEnergy to write mkm input files and json files
'''

from pathlib import Path
import yaml
from cathub.cat_energy.core import EnergeticAnalysis, write_energies


class ReturnValues:
    '''
    dummy class to return objects from methods defined inside other classes
    '''
    def __init__(self, input_dict):
        for key, value in input_dict.items():
            setattr(self, key, value)


def run_cat_energy(src_path, desired_surface, desired_facet, write_species,
                   write_mkm_input_files, verbose, latex):
    '''
    Function to run CatEnergy module
    '''
    config_file_path = f'{desired_surface}_{desired_facet}.yaml'
    with open(config_file_path, 'r', encoding='utf-8') as stream:
        try:
            params = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    config_params = ReturnValues(params)

    db_filepath = src_path / config_params.db_filename
    gas_jsondata_filepath = src_path / config_params.gas_jsondata_filename
    ads_jsondata_filepath = src_path / config_params.ads_jsondata_filename
    ts_jsondata_filepath = src_path / config_params.ts_jsondata_filename
    rxn_expressions_filepath = (src_path
                                / config_params.rxn_expressions_filename)

    system_parameters = {}
    system_parameters['desired_surface'] = desired_surface
    system_parameters['desired_facet'] = desired_facet
    system_parameters['temp'] = config_params.temp
    system_parameters['pH'] = config_params.pH
    for she_voltage in config_params.she_voltage_range:
        system_parameters['u_she'] = she_voltage
        analysis = EnergeticAnalysis(
            db_filepath, system_parameters, config_params.reference_gases, config_params.dummy_gases,
            config_params.dft_corrections_gases, config_params.fake_ads, config_params.beef_dft_helmholtz_offset,
            config_params.external_effects, config_params.facet_conditional, gas_jsondata_filepath,
            ads_jsondata_filepath, ts_jsondata_filepath, rxn_expressions_filepath, config_params.ts_data,
            write_mkm_input_files, verbose, latex)

        df_out = write_energies(analysis, write_species)
        df_out = df_out[df_out.site_name != 'gas']
    return None


def main():
    '''
    Main function
    '''
    cwd = Path.cwd()
    desired_surface = 'Cu'
    desired_facet = '100'
    write_species = {'gases': 1,
                     'adsorbates': 1,
                     'transition_states': 1,
                     }
    write_mkm_input_files = 1
    verbose = 1
    latex = 0
    run_cat_energy(cwd, desired_surface, desired_facet, write_species,
                   write_mkm_input_files, verbose, latex)
    return None

if __name__ == "__main__":
    main()
