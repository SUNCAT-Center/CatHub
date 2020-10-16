import pandas as pd
from sqlalchemy import create_engine 
from ase.db import connect

from tabulate import tabulate

class CatmapInterface():
    def __init__(self, table_name, filename):
        self.df = db_to_dataframe(table_name, filename)
        

    def write_input(self):
        # function to write input file
        a = 1
        
        
def db_to_dataframe(table_name, filename):
    "Read cathub .db file into pandas dataframe"

    # define sql url
    sql_url = 'sqlite:///' + str(filename)

    # SQLAlchemy connectable
    cnx = create_engine(sql_url).connect()

    # table will be returned as a dataframe
    df = pd.read_sql_table(table_name, cnx)
    return df

def write_energies(db_filepath, critical_density):
    "Write formation energies to energies.txt after applying free energy corrections"

    table_name = 'systems'
    df = db_to_dataframe(table_name, str(db_filepath))
    gas_ids = list(df.id[df.mass / df.volume < critical_density])

    db = connect(str(db_filepath))
    gas_select_rows = [list(db.select(id=gas_id))[0] for gas_id in gas_ids]
    surface, site, species, formation_energies = [], [], [], []
    for row in gas_select_rows:
        mass = row.mass
        volume = row.volume
        if mass / volume < critical_density:
            surface.append('None')
            site.append('gas')
        species.append(row.formula)
        formation_energies.append(f'{row.energy:.3f}')

    df = pd.DataFrame(list(zip(surface, site, species, formation_energies)),
                      columns=['Surface Name', 'Site Name', 'Species Name', 'Formation Energy'])

    energies_filepath = db_filepath.parent / 'energies.txt'
    with open(energies_filepath, 'w') as energies_file:
        df.to_string(energies_file, index=False)
    return None

