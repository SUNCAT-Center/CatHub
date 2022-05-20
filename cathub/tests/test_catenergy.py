#!/usr/bin/env python
'''
Test script for CatEnergy package
'''

from pathlib import Path
import unittest

import yaml

from cathub.cat_energy.conversion import formula_to_chemical_symbols, \
    get_catmap_style_species


DESIRED_SURFACE = 'Cu'
DESIRED_FACET = '100'
SRC_PATH = Path.cwd()

class ReturnValues:
    '''
    dummy class to return objects from methods defined inside other classes
    '''
    def __init__(self, input_dict):
        for key, value in input_dict.items():
            setattr(self, key, value)

class CatEnergyTestCase(unittest.TestCase):
    '''
    Test case definitions for CatEnergy package
    '''
    def test1_read_params(self):
        '''
        Test function to read parameter file for a desired configuration
        '''
        config_file_path = f'cat_energy/{DESIRED_SURFACE}_{DESIRED_FACET}.yaml'
        with open(config_file_path, 'r', encoding='utf-8') as stream:
            try:
                params = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        config_params = ReturnValues(params)
        return config_params

    def test2_chemical_formula_to_symbols(self):
        '''
        Test function to convert formula to dictionary of chemical symbols
        '''
        species_list = ['H2', 'CO', '2CCO', 'OCCOH']
        result_check_list = [{'H': 2},
                             {'C': 1, 'O': 1},
                             {'C': 4, 'O': 2},
                             {'O': 2, 'C': 2, 'H': 1}]
        chemical_symbols_dict_list = []
        for index, species_value in enumerate(species_list):
            chemical_symbols_dict_list[index] = formula_to_chemical_symbols(
                                                                species_value)

        for index, result_dict in enumerate(result_check_list):
            assert result_dict == chemical_symbols_dict_list[index]

    def test3_catmap_style_species(self):
        '''
        Test function to convert formula to CatMAP style species definition
        '''
        species_list = ['H2_g', '2CO*_t', 'CCO_t', 'H_g']
        result_check_list = [('H2gas', 1),
                             ('COstar', 2.0),
                             ('CCOstar', 1),
                             ('H2gas', 0.5)]
        catmap_species_list = []
        for index, species_value in enumerate(species_list):
            catmap_species_list[index] = get_catmap_style_species(
                                                                species_value)

        for index, result_tuple in enumerate(result_check_list):
            assert result_tuple == catmap_species_list[index]

if __name__ == '__main__':
    unittest.main()
