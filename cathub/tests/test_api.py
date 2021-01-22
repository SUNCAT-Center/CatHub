import os
import sys
import unittest
import tempfile
import pprint
import sqlite3
import json
import click
from cathub.cathubsql import CathubSQL

path = os.path.abspath(os.path.join(os.path.dirname(__file__)))


class ApiTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_remote_atoms(self):
        db = CathubSQL(user='postgres')

        dataframe = db.get_dataframe(pub_id='BajdichWO32018',
                                     include_atoms=True,
                                     elements=['Ti', 'W', 'O', 'H'])
        columns = ['surface_composition', 'reaction_energy', 'atoms']
        for c in columns:
            assert c in dataframe

    def test_remote_noatoms(self):
        db = CathubSQL(user='postgres')

        dataframe = db.get_dataframe(include_atoms=False,
                                     elements=['Ti', 'W', 'O', 'H'])
        assert dataframe[dataframe['reaction_id'] ==
                         34594]['reaction_energy'].values[0] == 4.59044414

        columns = ['surface_composition', 'reaction_energy']
        for c in columns:
            assert c in dataframe


if __name__ == '__main__':
    unittest.main()
