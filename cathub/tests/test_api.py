import os
import unittest

from cathub.query import CathubQuery

path = os.path.abspath(os.path.join(os.path.dirname(__file__)))


class ApiTestCase(unittest.TestCase):
    def setUp(self):
        self.client = CathubQuery()

    def test_remote_dataframe(self):
        df = self.client.get_reactions(pub_id='BajdichWO32018',
                                       n_results=1,
                                       include_structures=True)
        columns = ['surface_composition', 'reaction_energy',
                   'atoms_name', 'atoms_id']
        for c in columns:
            assert c in df.columns, 'Missing column: {}'.format(c)

    def test_remote_filter_elements(self):
        df = self.client.get_reactions(chemicalComposition='~Ti', n_results=2)
        assert len(df) > 0
        columns = ['surface_composition', 'reaction_energy']
        for c in columns:
            assert c in df.columns, 'Missing column: {}'.format(c)


if __name__ == '__main__':
    unittest.main()
