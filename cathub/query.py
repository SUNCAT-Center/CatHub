#!/usr/bin/env python
import re
import os
import io
import json
import requests
import pprint
import ase.db
import pandas as pd

from cathub.cathubsqlite import CathubSQLite
from cathub.postgresql import CathubPostgreSQL

all_columns = {'reactions': ['chemicalComposition', 'surfaceComposition',
                             'facet', 'sites', 'coverages', 'reactants',
                             'products', 'Equation',
                             'reactionEnergy', 'activationEnergy',
                             'dftCode', 'dftFunctional',
                             'username', 'pubId'],
               'publication': ['pubId', 'title', 'authors', 'journal',
                               'number', 'volume',
                               'pages', 'year', 'publisher', 'doi', 'tags'],
               'publications': ['pubId', 'title', 'authors', 'journal',
                                'volume', 'number',
                                'pages', 'year', 'publisher', 'doi', 'tags'],
               'reactionSystems': ['name', 'energyCorrection', 'aseId'],
               'publicationSystems': ['pubId', 'aseId']}


class CathubQuery:
    """
    Authenticated client for the Catalysis Hub GraphQL API.

    Parameters
    ----------
    api_key : str, optional
        API key for api.catalysis-hub.org. Falls back to the
        CATHUB_API_KEY environment variable if not provided.
    """

    _root = 'http://api.catalysis-hub.org/graphql'

    def __init__(self, api_key=None):
        self.api_key = api_key or os.environ.get('CATHUB_API_KEY')

    def _execute_graphql(self, query_string):
        print('Connecting to database at {}'.format(self._root))
        print('')
        print('Executing query:')
        print('')
        print(query_string)
        print('')
        print('Getting data from server...')
        print('')
        headers = {'X-API-Key': self.api_key} if self.api_key else {}
        data = requests.post(self._root, json={'query': query_string},
                             headers=headers)
        try:
            data = data.json()['data']
            print('Data fetched!')
        except BaseException:
            print(data.text)

        #for i, node in enumerate(data['reactions']['edges']):
        #    node = node['node']
        #    for key, value in list(node.items()):
        #        try:
        #            node[key] = json.loads(value)
        #        except (ValueError, TypeError):
        #            pass

        return data

    def query(self, table, columns, subtables, n_results, queries):
        query_string = graphql_query(table=table,
                                     subtables=subtables,
                                     columns=columns,
                                     n_results=n_results,
                                     queries=queries)
        return self._execute_graphql(query_string)
    """
    def get_reactions(self, columns='all', n_results=20, write_db=False,
                      **kwargs):
        # Get reactions from the server. Pass filter criteria as kwargs.
        if write_db or columns == 'all':
            columns = all_columns['reactions']
        queries = {}
        for key, value in kwargs.items():
            key = map_column_names(key)
            if key == 'distinct':
                if value in [True, 'True', 'true']:
                    queries.update({key: True})
                    continue
            if isinstance(value, (int, float)):
                queries.update({key: value})
            else:
                queries.update({key: '{0}'.format(value)})

        subtables = ['reactionSystems', 'publication'] if write_db else []
        data = self._query(table='reactions', subtables=subtables,
                           columns=columns, n_results=n_results,
                           queries=queries)


    """
    def get_publications(self, **kwargs):
        """Get publications from the server."""
        queries = {}
        for key, value in kwargs.items():
            key = map_column_names(key)
            try:
                value = int(value)
                queries.update({key: value})
            except (ValueError, TypeError):
                queries.update({key: '{0}'.format(value)})

        return self.query(table='publications',
                           columns=all_columns['publications'],
                           subtables=[],
                           n_results='all',
                           queries=queries)

    def get_structure(self, ase_id):
        queries = {'uniqueId': ase_id}

        data = self.query(table='systems',
                           subtables=[],
                           n_results=1,
                           columns=['Trajdata'],
                           queries=queries)
        edge = data.get('systems', {}).get('edges', [])
        if not data:
            return none

        import io
        trajdata_string = edge[0]['node']['Trajdata']
        trajfile = io.StringIO(trajdata_string)
        atoms = ase.io.read(trajfile, format='json')

        return atoms

    def get_reactions(self, pub_id=None, surface_composition=None,
                      facet=None, n_results='all', include_structures=False,
                      **kwargs):
        """
        Get reactions from the Catalysis Hub GraphQL API as a pandas DataFrame.

        Parameters
        ----------
        pub_id : str, optional
            Publication ID, e.g. 'PengRole2020'
        surface_composition : str, optional
            Match a specific surface composition
        facet : str, optional
            Match a specific surface facet
        n_results : int or 'all'
            Maximum number of results to fetch
        include_structures: Bool
            Fetch ASE Atoms objects
        **kwargs
            Additional GraphQL query parameters (camelCase or mapped names)
        """
        queries = {}
        if pub_id is not None:
            queries['pubId'] = pub_id
        if surface_composition is not None:
            queries['surfaceComposition'] = surface_composition
        if facet is not None:
            queries['facet'] = facet
        for key, value in kwargs.items():
            queries[map_column_names(key)] = value

        columns = [c for c in all_columns['reactions'] if c != 'Equation']
        data = self.query(table='reactions',
                          columns=columns,
                          subtables=['reactionSystems', 'publication'],
                          n_results=n_results,
                          queries=queries)

        edges = data.get('reactions', {}).get('edges', [])

        if not edges:
            print('No reactions found.')
            return pd.DataFrame()

        rows = []
        for edge in edges:
            node = edge['node']
            row = {}
            for key, value in node.items():
                if key in ('reactionSystems', 'publication'):
                    continue
                row[convert(key)] = value

            pub = node.get('publication') or {}
            row['doi'] = pub.get('doi')

            reaction_systems = node.get('reactionSystems') or []
            row['atoms_name'] = [rs['name'] for rs in reaction_systems]
            row['atoms_id'] = [rs['aseId'] for rs in reaction_systems]

            row['equation'] = get_equation(
                node.get('reactants', {}), node.get('products', {}))

            rows.append(row)
            if include_structures:
                row['atoms'] = [self.get_structure(ase_id=rs['aseId']) for rs
                                in reaction_systems]

        return pd.DataFrame(rows)

def graphql_query(table='reactions',
                  subtables=[],
                  columns=['chemicalComposition',
                           'reactants',
                           'products'],
                  n_results=10,
                  queries={}):

    statement = '{'
    statement += '{}('.format(table)
    if n_results != 'all':
        statement += 'first: {}'.format(n_results)
    for key, value in queries.items():
        if isinstance(value, str):
            statement += ', {}: "{}"'.format(key, value)
        elif isinstance(value, bool):
            if value:
                statement += ', {}: true'.format(key)
            else:
                statement += ', {}: false'.format(key)
        else:
            statement += ', {}: {}'.format(key, value)

    statement += ') {\n'
    statement += ' totalCount\n  edges {\n    node { \n'
    for column in columns:
        column = map_column_names(column)
        statement += '      {}\n'.format(column)
    for subtable in subtables:
        statement += '      {}'.format(subtable)
        statement += '{\n'
        for column in all_columns[subtable]:
            statement += '        {}\n'.format(column)
        statement += '      }\n'
    statement += '    }\n'
    statement += '  }\n'
    statement += '}}'

    return statement


def map_column_names(column):
    mapping = {'surface': 'chemicalComposition'}

    if column in mapping:
        return mapping[column]
    else:
        return column


def convert(name):
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


def get_equation(reactants, products):
    """Format a reaction equation string from reactants/products dicts."""
    r_str = ''
    for j, side in enumerate([reactants, products]):
        if isinstance(side, str):
            side = json.loads(side)
        i = 0
        for name in sorted(side.keys()):
            pf = side[name]
            name = name.replace('gas', '(g)').replace('star', '*')
            if i > 0 and not pf < 0:
                r_str += ' +'
            if pf == 1:
                r_str += ' {}'.format(name)
            elif pf == -1:
                r_str += ' - {}'.format(name)
            elif pf < 0:
                r_str += ' - {}{}'.format(abs(pf), name)
            else:
                r_str += ' {}{}'.format(pf, name)
            i += 1
        if j == 0:
            r_str += ' ->'
    return r_str.lstrip(' ')


if __name__ == '__main__':
    query = query(table='reactions',
                  columns=['chemicalComposition',
                           'reactants',
                           'products'],
                  n_results=10,
                  queries={'chemicalComposition': "~Pt"})
