#!/usr/bin/env python
'''
Test script for input output related functions
'''

import json
from pathlib import Path
import unittest
from collections import Counter

import pandas as pd

from cathub.cathubsql import CathubSQL


PUB_ID = 'PengRole2020'
DB_FILEPATH = Path.cwd() / 'io' / PUB_ID / '.db'

class IOTestCase(unittest.TestCase):
    '''
    Test case definitions for IO functions
    '''
    def test1_read_db_file(self):
        '''
        Test function to read ASE DB file from Catalysis Hub vs. local dir
        '''

        # Data from local cathub .db file
        db_local = CathubSQL(filename=DB_FILEPATH)
        df_local = db_local.get_dataframe()

        # To get data on catalysis-hub.org
        db_remote = CathubSQL()
        df_remote = db_remote.get_dataframe(pub_id=PUB_ID)

        select_columns = [
            'chemical_composition', 'surface_composition', 'facet',
            'reaction_energy', 'activation_energy', 'dft_code',
            'dft_functional', 'username', 'pub_id', 'equation']
        assert df_local[select_columns].equals(df_remote[select_columns])

        select_columns = ['reactants', 'products']
        if isinstance(df_local.loc[0, select_columns[0]], str):
            list1 = [json.loads(d) for d in df_local[select_columns[0]]]
            list2 = [json.loads(d) for d in df_local[select_columns[1]]]
            df = pd.DataFrame(list(zip(list1, list2)), columns=select_columns)
            assert df.equals(df_remote[select_columns])
        else:
            assert df_local[select_columns].equals(df_remote[select_columns])

        select_column = 'atoms_name'
        dict_local = {
                select_column : [Counter(d) for d in df_local[select_column]]}
        dict_remote = {
                select_column : [Counter(d) for d in df_remote[select_column]]}
        df_temp_local = pd.DataFrame(dict_local)
        df_temp_remote = pd.DataFrame(dict_remote)
        assert df_temp_local.equals(df_temp_remote)
