import os
import json
import pandas
import numpy as np
from cathub.postgresql import CathubPostgreSQL, get_value_str


init_commands = [

    """CREATE TABLE publication (
    id SERIAL PRIMARY KEY,
    pub_id text UNIQUE,
    title text,
    authors jsonb,
    journal text,
    volume text,
    number text,
    pages text,
    year smallint,
    publisher text,
    doi text,
    tags jsonb,
    stime numeric
    );""",

    """CREATE TABLE material (
    mat_id SERIAL PRIMARY KEY,
    pub_id text REFERENCES publication (pub_id),
    composition text,
    arrangement text,
    icsd_ids integer[],
    space_group text,
    lattice_parameter text,
    morphology text,
    notes jsonb
    );""",

    """CREATE TABLE sample (
    sample_id SERIAL PRIMARY KEY,
    mat_id integer REFERENCES material (mat_id),
    pub_id text REFERENCES publication (pub_id),
    data jsonb 
    );""",

    """CREATE TABLE xps (
    mat_id integer PRIMARY KEY REFERENCES material (mat_id),
    sample_id integer REFERENCES sample (sample_id),
    type text, 
    binding_energy DOUBLE PRECISION[],
    intensity DOUBLE PRECISION[]
    );""",  # type: pre, post

    """CREATE TABLE xrd (
    mat_id integer PRIMARY KEY REFERENCES material (mat_id),
    type text,
    degree DOUBLE PRECISION[],
    intensity DOUBLE PRECISION[]
    );""",

    """CREATE TABLE echemical (
    id SERIAL PRIMARY KEY,
    type text,
    total_time numeric, 
    time DOUBLE PRECISION[],
    potential DOUBLE PRECISION[],
    current DOUBLE PRECISION[],
    sample_id integer REFERENCES sample (sample_id)
    );""",  # type: CV_initial, CV_end, CP, CA
]

mat_columns = ['composition', 'arrangement', 'ICSD_ID', 'space_group',
               'lattice_parameter_a(nm)', 'morphology']  # , 'XPS_ID', 'XRD(CuKa)ID']

sample_columns = ['mat_id', 'pub_id', 'data']
xps_columns = ['mat_id', 'sample_id', 'type', 'binding_energy', 'intensity']
xrd_columns = ['mat_id', 'type', 'degree', 'intensity']
echemical_columns = ['sample_id', 'type',
                     'total_time', 'time', 'potential', 'current']


class ExpSQL(CathubPostgreSQL):

    def __init__(self, user='experimental', schema='experimental', password=None):
        super().__init__(schema=schema, user=user, password=password)
        #self.create_user('experimental', row_limit = None)

    def _initialize(self, con):
        if self.initialized:
            return
        cur = con.cursor()

        set_schema = 'SET search_path TO {0};'\
                     .format(self.schema)
        cur.execute(set_schema)

        self.stdout.write(
            "_initialize set schema to {self.schema}\n".format(**locals()))

        cur.execute("""SELECT to_regclass('sample');""")

        if cur.fetchone()[0] is not None:
            self.initialized = True
            return

        for init_command in init_commands:
            self.stdout.write(init_command + '\n')
            cur.execute(init_command)
        self.initialized = True

    def write(self, dataframes, doi=None):
        con = self.connection or self._connect()
        cur = con.cursor()
        self._initialize(con)

        main_table = dataframes['Tabulated data']
        main_table = main_table[main_table['DOI'] == doi]
        print(main_table)
        XPS_table = dataframes['XPS']
        XRD_table = dataframes['XRD']
        CV_table = dataframes['CV']
        CVend_table = dataframes['CVend']
        Stability_table = dataframes['Stability Test']

        import sys

        XPS_ids = []
        for index, row in main_table.iterrows():
            XPS_id = row['XPS_ID']
            XRD_id = row['XRD(CuKa)ID']
            CV_init_id = row['CV_intial_ID']
            CV_end_id = row['CV_end_ID']
            Stab_id = row['Stability_Test_ID']
            XPS_post_id = row['XPS_post_test_ID']

            if pandas.isnull(XPS_id) or pandas.isnull(XRD_id):
                continue

            elif not XPS_id in XPS_table:
                continue

            elif not XRD_id in XRD_table:
                continue

            print(XRD_id, XRD_id, CV_init_id)
            if not XPS_id in XPS_ids:
                XPS_ids += [XPS_id]

                # Material
                mat_value_list = [row[c] for c in mat_columns]
                mat_value_list[2] = [int(m) for m in str(
                    mat_value_list[2]).split('/')]
                key_str = ', '.join(mat_columns)
                key_str = key_str.replace(
                    '_a(nm)', '').replace('ICSD_ID', 'icsd_ids')

                value_str = get_value_str(mat_value_list)

                query_mat = \
                    """INSERT INTO material ({}) VALUES ({}) RETURNING mat_id
                    """.format(key_str, value_str)
                print(query_mat)
                cur.execute(query_mat)
                mat_id = cur.fetchone()[0]

                # XPS
                E = XPS_table[XPS_id].values[1:]
                idx = [not pandas.isnull(e) for e in E]
                E = np.array(E)[idx]
                I = np.array(XPS_table[XPS_id + '.1'].values[1:])[idx]

                key_str = ', '.join(xps_columns)
                xps_value_list = [mat_id, None, 'pre', list(E), list(I)]
                value_str = get_value_str(xps_value_list)

                query_xps = \
                    """INSERT INTO xps ({}) VALUES ({})
                    """.format(key_str, value_str)
                cur.execute(query_xps)

                # XRD
                Deg = XRD_table[XRD_id].values[1:]
                idx = [not (pandas.isnull(d) or d == '--') for d in Deg]
                Deg = np.array(Deg)[idx]
                I = np.array(XRD_table[XRD_id + '.1'].values[1:])[idx]

                key_str = ', '.join(xrd_columns)
                xrd_value_list = [mat_id, 'pre', list(Deg), list(I)]
                value_str = get_value_str(xrd_value_list)

                query_xrd = \
                    """INSERT INTO xrd ({}) VALUES ({})
                    """.format(key_str, value_str)
                cur.execute(query_xrd)

            # Sample table
            sample_dict = row.to_dict()
            clean_dict = {}
            for k, v in sample_dict.items():
                if pandas.isnull(v):
                    v = 'NULL'
                if k == 'ICSD ID':
                    v = [int(m) for m in v.split('/')]
                k = k.replace('(', '_').replace(')', '_').replace('=', '_').replace(',', '_').replace(
                    '/', '_o_').replace('%', 'perc').replace('+', 'plus').replace('-', 'minus').replace('.', 'd')
                clean_dict[k] = v

            key_str = ', '.join(sample_columns)
            value_str = get_value_str([mat_id, None, json.dumps(clean_dict)])

            query_sample = \
                """INSERT INTO sample ({}) VALUES ({}) RETURNING sample_id
                 """.format(key_str, value_str)

            cur.execute(query_sample)
            sample_id = cur.fetchone()[0]
            #print(mat_id, sample_id)

            """ echemical Table"""
            key_str = ', '.join(echemical_columns)
            if not pandas.isnull(CV_init_id):
                if not CV_init_id in CV_table:
                    print('ID not found cvinit!', CV_init_id)
                else:
                    V = CV_table[CV_init_id].values[1:]
                    idx = [not (pandas.isnull(v) or v == '--') for v in V]
                    V = np.array(V)[idx]
                    I = np.array(CV_table[CV_init_id + '.1'].values[1:])[idx]

                    cv_value_list = [sample_id, 'CV_initial',
                                     row['time_tested(h)'], None, list(V), list(I)]
                    value_str = get_value_str(cv_value_list)

                    query_cv = \
                        """INSERT INTO echemical ({}) VALUES ({})
                        """.format(key_str, value_str)
                    # print(query_cv)
                    cur.execute(query_cv)

            if not pandas.isnull(CV_end_id):
                if not CV_end_id in CVend_table:
                    print('Id not found cvend!', CV_end_id)
                else:
                    V = CVend_table[CV_end_id].values[1:]
                    idx = [not (pandas.isnull(v) or v == '--') for v in V]
                    V = np.array(V)[idx]
                    I = np.array(CVend_table[CV_end_id + '.1'].values[1:])[idx]

                    cv_value_list = [sample_id, 'CV_end',
                                     row['time_tested(h)'], None, list(V), list(I)]
                    value_str = get_value_str(cv_value_list)

                    query_cv = \
                        """INSERT INTO echemical ({}) VALUES ({})
                        """.format(key_str, value_str)
                    # print(query_cv)
                    cur.execute(query_cv)

            # if not pandas.isnull(Stab_id):
                #V = Stab_table[Sta_id].values[1:]
                #idx = [not (pandas.isnull(v) or v=='--') for v in V]
                #V = np.array(V)[idx]
                #I = np.array(CVend_table[CV_end_id + '.1'].values[1:])[idx]
               #
               # cv_value_list = [sample_id, 'CV_end', row['time_tested(h)'], None, list(V), list(I)]
               # value_str = get_value_str(cv_value_list)
               #
               # query_cv = \
               #     """INSERT INTO echemical ({}) VALUES ({})
               #     """.format(key_str, value_str)
               # print(query_cv)
               # cur.execute(query_cv)

            # break

        con.commit()
        con.close()

        #print(mat_id, sample_id)
        # sys.exit()

        #xps_id = row['XPS_ID']

        # print(xps_id)
        # print(mat_id)
        # return

        #


def load_table(filename, tablename):
    skiprows = {'Tabulated data': [0, 1]}
    data = {}

    data = pandas.read_excel(filename, sheet_name=tablename,
                             skiprows=skiprows.get(tablename),
                             mangle_dupe_cols=True)

    return data


def clean_column_names(dataframe):
    # remove empty spaces from column names
    columns = dataframe.columns

    column_rename = {}
    for c in columns:
        c_old = c
        c_new = c_old.replace(' (', '(').replace(') ', ')')
        c_new = '_'.join(c_new.split(' '))
        c_new = c_new.strip('_')
        column_rename[c_old] = c_new

    dataframe = dataframe.rename(columns=column_rename)
    return dataframe


def clear_duplicate_columns(data):

    values = np.array(data.values[1:, :], dtype=float)
    print(values)
    columns = data.columns

    duplicates_dict = {}
    for i, c in enumerate(columns):
        if data.values[0, i] == 'Binding Energy':
            continue

        duplicate_columns = []
        for v in duplicates_dict.values():
            duplicate_columns += v

        if c in duplicate_columns:
            continue

        duplicates_idx = \
            np.where([np.all(np.isclose(values[:, i], values[:, j], equal_nan=True))
                      for j in range(0, len(columns))])[0]

        duplicates_idx = duplicates_idx[1:]

        duplicates = [columns[k].rstrip('.1') for k in duplicates_idx]
        duplicates_dict[c.rstrip('.1')] = duplicates

    all_duplicates = []
    for v in duplicates_dict.values():
        all_duplicates += v

    all_duplicates += [d + '.1' for d in all_duplicates]

    all_duplicates = set(all_duplicates)

    data = data.drop(columns=all_duplicates)

    return data, duplicates_dict


def clear_duplicate_rows(dataframe):

    data = dataframe.groupby('XPS_ID', as_index=False).first()
    return data
    # print(data.values)
    # print('hello')
    #import sys
    # sys.exit()

    values = np.array(dataframe.values[1:, :])
    print(values)
    duplicates_dict = {}

    print(dataframe)

    for i in range(len(values)):
        duplicate_rows = []
        for v in duplicates_dict.values():
            # print(v)
            duplicate_rows += list(v)
        if i in duplicate_rows:
            continue

        duplicates_idx = \
            np.where([np.all(values[i] == values[j])
                      for j in range(0, len(values))])[0]

        duplicates_idx = duplicates_idx[1:]

        duplicates_dict[i] = duplicates_idx

    all_duplicates = []
    for v in duplicates_dict.values():
        all_duplicates += list(v)

    print(all_duplicates)
    dataframe = dataframe.drop(all_duplicates)

    import sys
    sys.exit()
    return dataframe, duplicates_dict


if __name__ == '__main__':

    dataframes = {}
    filename = '/Volumes/GoogleDrive/Shared drives/SUNCAT_experimental_data/OxR _HxR Jaramillo Database.xlsx'

    for table in ['Tabulated data', 'XPS', 'XRD', 'CV', 'CVend', 'Stability Test']:
        data = load_table(filename=filename, tablename=table)
        data = clean_column_names(data)
        dataframes[table] = data

    print(dataframes['Tabulated data']['DOI'].values)

    DB = ExpSQL(user='experimental', password='hGaPjHfo')  # 'catroot') #

    DB.write(dataframes, doi='10.1021/acscatal.0c02252')

    # print(data.columns)
    # print(len(data_material.to_dict('records')))

    # CREATED USER experimental WITH PASSWORD hGaPjHfo
