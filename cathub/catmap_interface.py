from pandas import read_sql_table
from sqlalchemy import create_engine 


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
    df = read_sql_table(table_name, cnx)
    return df

