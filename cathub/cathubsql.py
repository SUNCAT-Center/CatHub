from sqlalchemy import create_engine
from pandas import read_sql
import ase.db
import ase.visualize


class CathubSQL:
    """
    Generallized interface to CatHub local and server SQL databases
    """

    def __init__(self, filename=None):

        if filename is not None:
            sql_url = 'sqlite:///' + str(filename)
        else:
            sql_url = 'postgresql://catvisitor:eFjohbnD57WLYAJX@catalysishub.c8gwuc8jwb7l.us-west-2.rds.amazonaws.com:5432/catalysishub'

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

    def get_dataframe(self, pub_id=None, include_atoms=False):
        """
        Get pandas dataframe containing reactions for dataset

        Parameters:

        pub_id: str
           catalysis-hub dataset pub id, such as 'PengRole2020'
        include_atoms: bool
           To include atoms object or not. If set to false. 
        """

        # Query SQL table to get reaction, publication, and structure info.
        query = \
            """select r.*, rs.name, rs.ase_id, p.doi from reaction as r 
        left join
        reaction_system as rs on r.id = rs.id
        left join
        publication as p on r.pub_id=p.pub_id
        """

        if pub_id is not None:
            query += " where r.pub_id='{}'".format(pub_id)

        con = self.connection or self._connect()
        print('Querying database')
        dataframe = read_sql(query, con)
        print('Done')
        if self.connection is None:
            con.close()

        # load ase atoms objects to add to dataframe
        if include_atoms:
            atoms_list = []
            id_to_atoms = {}
            with ase.db.connect(self.sql_url) as ase_db:
                total = ase_db.count('pub_id={}'.format(pub_id))
                print('Fetching {} atomic structures'.format(total))
                for i, row in enumerate(ase_db.select('pub_id={}'.format(pub_id))):
                    if (i+1) % 10 == 0:
                        print('  {}/{}'.format(i+1, total))
                    id_to_atoms[row.unique_id] = row.toatoms()

            for id in dataframe['ase_id'].values:
                atoms_list += [id_to_atoms[id]]

            dataframe['atoms'] = atoms_list

        # group by reaction id and aggregate structure columns to list
        dataframe = dataframe.drop(columns=['textsearch'])
        dataframe = dataframe.rename(columns={'id': 'reaction_id',
                                              'name': 'atoms_name',
                                              'ase_id': 'atoms_id'}, index={'id': 'id'})

        columns_group = {}

        for c in dataframe.columns.values:
            columns_group[c] = 'first'
        columns_group['atoms_name'] = list
        columns_group['atoms_id'] = list
        if include_atoms:
            columns_group['atoms'] = list

        dataframe = dataframe.groupby(['reaction_id'], as_index=False)\
                             .agg(columns_group)

        equations = []
        print(dataframe[['reactants', 'products']].values)
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


def get_equation(reactants, products):
    r_str = ''
    print(reactants, products)
    for j, side in enumerate([reactants, products]):
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
