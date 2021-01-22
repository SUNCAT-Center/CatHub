from sqlalchemy import create_engine
import sqlite3
from pandas import read_sql
import json
import ase.db
import ase.visualize

from cathub.postgresql import CathubPostgreSQL


class CathubSQL:
    """
    Generallized interface to CatHub local and server SQL databases
    """

    def __init__(self, user='apiuser', filename=None):

        if filename is not None:
            sql_url = 'sqlite:///' + str(filename)
            self.backend = 'sqlite'
        else:
            sql_url = CathubPostgreSQL(user).server_name
            self.backend = 'postgres'

        self.sql_url = sql_url

        self.connection = None

    def _connect(self):
        engine = create_engine(self.sql_url)
        return engine.connect()

    def __enter__(self):
        """Set connection upon entry using with statement"""
        assert self.connection is None
        self.connection = self._connect()
        return self

    def __exit__(self, exc_type, exc_value, tb):
        """Commit changes upon exit"""
        if exc_type is None:
            self.connection.commit()
        else:
            self.connection.rollback()
        self.connection.close()
        self.connection = None

    def get_dataframe(self,
                      include_atoms=False,
                      pub_id=None,
                      reactants=None,
                      products=None,
                      elements=None,
                      surface_composition=None,
                      facet=None):
        """
        Get pandas dataframe containing reactions for dataset

        Parameters:

        include_atoms: bool or filename
           True: Atoms objects are included in dataframe.
           filename.db: Atoms objects are writen to a local ASE db file.
        pub_id: str
           Catalysis-hub dataset pubid, such as 'PengRole2020',
           found at catalysis-hub.or/publications
        reactants: list of str
           list of species on the left-hand side of equation
           for example: ['CH4gas' + 'H2gas']
        products: list of str
           list of species on the right-hand side of equation
           for example: ['CH3*'] or [CH3star]
        elements: list of str
           List of atomic elements in the surface, for example: ['Cu', 'Zn'].
           Use '-' in front, as in ['Cu', '-Zn'], to exclude elements
        surface_composition: str
           Match a specific surface composition        
        facet: str
           Match a specific surface facet
        """

        # Query SQL table to get reaction, publication, and structure info.
        query = \
            """SELECT r.*, rs.name, rs.ase_id, p.doi FROM reaction as r 
LEFT JOIN
reaction_system as rs on r.id = rs.id
LEFT JOIN
publication as p on r.pub_id=p.pub_id"""

        query += get_sql_query(backend=self.backend,
                               pub_id=pub_id,
                               reactants=reactants,
                               products=products,
                               elements=elements,
                               surface_composition=surface_composition,
                               facet=facet)

        con = self.connection or self._connect()
        print('Querying database\n')
        dataframe = read_sql(query, con)
        if self.connection is None:
            con.close()

        if len(dataframe) == 0:
            print('No reactions in database for {} -> {} and elements={}'
                  .format(reactants, products, elements))
            print(query)
            return dataframe

        # load ase atoms objects to add to dataframe
        if include_atoms:
            atoms_list = []
            id_to_atoms = {}
            with ase.db.connect(self.sql_url.replace('sqlite:///', '')) as ase_db:
                if isinstance(include_atoms, str):
                    with ase.db.connect(include_atoms) as local_db:
                        for id in set(dataframe['ase_id'].values):
                            row = ase_db.get(unique_id=id)
                            try:
                                local_db.write(row, data=row.data)
                            except sqlite3.IntegrityError:
                                # pass if structure allready downloaded
                                pass
                else:
                    for id in set(dataframe['ase_id'].values):
                        id_to_atoms[id] = \
                            ase_db.get(unique_id=id).toatoms()

                    for id in dataframe['ase_id'].values:
                        atoms_list += [id_to_atoms[id]]

                    dataframe['atoms'] = atoms_list

        # group by reaction id and aggregate structure columns to list
        if 'textsearch' in dataframe.columns:
            dataframe = dataframe.drop(columns=['textsearch'])

        dataframe = dataframe.rename(columns={'id': 'reaction_id',
                                              'name': 'atoms_name',
                                              'ase_id': 'atoms_id'}, index={'id': 'id'})

        columns_group = {}

        for c in dataframe.columns.values:
            columns_group[c] = 'first'

        if 'atoms_name' in dataframe.columns.values:
            columns_group['atoms_name'] = list
            columns_group['atoms_id'] = list
        if include_atoms == True:
            columns_group['atoms'] = list

        dataframe = dataframe.groupby(['reaction_id'], as_index=False)\
                             .agg(columns_group)

        equations = []
        for reactants, products in dataframe[['reactants', 'products']].values:
            equations += [get_equation(reactants, products)]

        dataframe['equation'] = equations

        return dataframe

    def get_atoms_for_reaction(self, reaction_id):
        con = self.connection or self._connect()
        print('Querying database')
        query = "select ase_id from reaction_system where id='{}'".format(
            reaction_id)

        ids = con.execute(query)
        if self.connection is None:
            con.close()

        atoms_list = []
        with ase.db.connect(self.sql_url) as ase_db:
            for unique_id in ids:
                for row in ase_db.select(unique_id=unique_id[0]):
                    atoms_list += [row.toatoms()]

        return atoms_list

    def get_atoms_for_publication(self, pub_id):
        atoms_list = []
        with ase.db.connect(self.sql_url) as ase_db:
            total = ase_db.count('pub_id={}'.format(pub_id))
            print('Fetching {} atomic structures'.format(total))
            for i, row in enumerate(ase_db.select('pub_id={}'.format(pub_id))):
                if (i+1) % 10 == 0:
                    print('  {}/{}'.format(i+1, total))
                    atoms_list += [row.toatoms()]

        return atoms_list

    def get_atoms_for_id(self, atoms_id=None):
        """Get atoms for atoms_id"""

        if not isinstance(atoms_id, list):
            atoms_id = [atoms_id]

        with ase.db.connect(self.sql_url) as ase_db:
            for unique_id in atoms_id:
                for row in ase_db.select(unique_id=unique_id):
                    atoms_list += [row.toatoms()]

        return atoms_list


def get_sql_query(backend='postgres',
                  pub_id=None,
                  reactants=None,
                  products=None,
                  elements=None,
                  surface_composition=None,
                  facet=None):
    query = ''

    if pub_id is not None:
        query += " \nWHERE r.pub_id='{}'".format(pub_id)

    reaction_side = ['reactants', 'products']
    for i, reactant_list in enumerate([reactants, products]):
        if reactant_list is not None:
            if not 'WHERE' in query:
                query += ' \nWHERE '
            else:
                query += ' \nAND '
            query_list = []
            for species in reactant_list:
                species = species.replace(
                    '*', 'star').replace('(g)', 'gas')
                if backend == 'postgres':
                    query_list += ["r.{} ? '{}'".format(
                        reaction_side[i], species)]
                else:
                    query_list += ["r.{} like '%{}%'".format(
                        reaction_side[i], species)]
            query += ' AND '.join(query_list)
    if elements is not None:
        if not 'WHERE' in query:
            query += ' \nWHERE '
        else:
            query += ' \nAND '
        query_list = []
        for e in elements:
            if e[0] == '-':
                e = e[1:]
                query_list += [
                    "r.chemical_composition not like '%{}%'".format(e)]
            else:
                query_list += [
                    "r.chemical_composition like '%{}%'".format(e)]
        query += ' \nAND '.join(query_list)
    if surface_composition is not None:
        if not 'WHERE' in query:
            query += ' \nWHERE '
        else:
            query += ' \nAND '
        query += "r.surface_composition = '{}' or surface_composition like '{}-%'"\
            .format(surface_composition, surface_composition)
    if facet is not None:
        if not 'WHERE' in query:
            query += ' \nWHERE '
        else:
            query += ' \nAND '
        query += "r.facet ilike '{}%'".format(facet)

    return query


def get_equation(reactants, products):
    r_str = ''
    for j, side in enumerate([reactants, products]):
        if isinstance(side, str):
            side = json.loads(side)
        i = 0
        for name in sorted(side.keys()):
            pf = side[name]
            name = name.replace('gas', '(g)').replace('star', '*')
            if i > 0 and not pf < 0:
                r_str += ' + '
            if pf == 1:
                r_str += '{}'.format(name)
            elif pf == -1:
                r_str += ' -{}'.format(name)
            else:
                r_str += '{}{}'.format(pf, name)

            i += 1
        if j == 0:
            r_str += ' -> '

    return r_str
