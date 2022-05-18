#!/usr/bin/env python
'''
Test script for CatEnergy package
'''

from pathlib import Path
import unittest

import yaml


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
        config_file_path = f'{DESIRED_SURFACE}_{DESIRED_FACET}.yaml'
        with open(config_file_path, 'r', encoding='utf-8') as stream:
            try:
                params = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)
        config_params = ReturnValues(params)
        return config_params

if __name__ == '__main__':
    unittest.main()
