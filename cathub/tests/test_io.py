#!/usr/bin/env python
'''
Test script for input/output related functions
'''
import unittest

from cathub.query import CathubQuery

PUB_ID = 'PengRole2020'

EXPECTED_COLUMNS = [
    'chemical_composition', 'surface_composition', 'facet',
    'reaction_energy', 'activation_energy', 'dft_code', 'dft_functional',
    'username', 'pub_id', 'equation', 'atoms_name', 'atoms_id', 'doi',
]


class IOTestCase(unittest.TestCase):
    def setUp(self):
        self.client = CathubQuery()

    def test_get_dataframe(self):
        df = self.client.get_dataframe(pub_id=PUB_ID)
        assert len(df) > 0, 'DataFrame should not be empty'
        for col in EXPECTED_COLUMNS:
            assert col in df.columns, 'Missing column: {}'.format(col)

    def test_get_dataframe_equation(self):
        df = self.client.get_dataframe(pub_id=PUB_ID)
        assert df['equation'].notna().all(), 'All rows should have an equation'
        assert df['equation'].str.contains('->').all(), \
            'All equations should contain ->'

    def test_get_dataframe_atoms_lists(self):
        df = self.client.get_dataframe(pub_id=PUB_ID)
        assert df['atoms_name'].apply(lambda x: isinstance(x, list)).all()
        assert df['atoms_id'].apply(lambda x: isinstance(x, list)).all()
