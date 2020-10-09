from pandas import read_sql_table

class CatmapInterface():
    def __init__(self, filename):
        self.df = db_to_dataframe(filename)
        

    def write_input(self):
        # function to write input file
        a = 1
        
        
def db_to_dataframe(filename):
    "Read cathub .db file into pandas dataframe"
    pd = read_sql_table(filename)

    return pd

    
