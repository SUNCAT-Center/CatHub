import warnings
from cathub.tools import get_pub_id, doi_request
from cathub.postgresql import CathubPostgreSQL, get_value_str, get_key_list, get_key_str
from sqlalchemy import create_engine
import pylab as p
import numpy as np
from pandas import read_sql
import sys
import os
import json
import pandas



warnings.filterwarnings('ignore', '.*SQL.*', )

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
    icdd_ids text[],
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
    id SERIAL PRIMARY KEY,
    mat_id integer REFERENCES material (mat_id),
    sample_id integer REFERENCES sample (sample_id),
    type text,
    binding_energy DOUBLE PRECISION[],
    intensity DOUBLE PRECISION[]
    );""",  # type: pre, post

    """CREATE TABLE xrd (
    id SERIAL PRIMARY KEY,
    mat_id integer PRIMARY KEY REFERENCES material (mat_id),
    type text,
    degree DOUBLE PRECISION[],
    intensity DOUBLE PRECISION[]
    );""",

    """CREATE TABLE xas_ledge (
    id SERIAL PRIMARY KEY,
    mat_id integer REFERENCES material (mat_id),
    type text,
    energy DOUBLE PRECISION[],
    intensity DOUBLE PRECISION[]
    );""",

    """CREATE TABLE xas_kedge (
    id SERIAL PRIMARY KEY,
    mat_id integer REFERENCES material (mat_id),
    type text,
    energy DOUBLE PRECISION[],
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

mat_columns = ['pub_id', 'composition', 'arrangement', 'ICSD_ID', 'ICDD_ID', 'space_group',
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

    def get_dataframe(self, table='sample', sample_ids=None, mat_ids=None, pub_id=None,
                      include_replicates=False):
        """
        Get pandas dataframe containing experimental data

        Parameters:

        """

        # Query SQL table to get reaction, publication, and structure info.
        query = \
            """SELECT * from experimental.{}""".format(table)
        if sample_ids is not None:
            query += " \nWHERE sample_id in ({})".format(
                ','.join([str(i) for i in sample_ids]))
        elif mat_ids is not None:
            query += " \nWHERE mat_id in ({})".format(
                ','.join([str(i) for i in mat_ids]))
        if pub_id is not None:
            if table in ['publication', 'sample', 'material']:
                query += " \nWHERE pub_id='{}'".format(pub_id)
                if table == 'sample' and not include_replicates:
                    query += "and (data->>'replicate' is null or data->>'replicate' in ('1', '-'))"
            elif table in ['xps', 'xrd']:
                query += " \nWHERE mat_id in (select mat_id from material where pub_id='{}')".format(
                    pub_id)
            elif table == 'echemical':
                if not include_replicates:
                    query += " \nWHERE sample_id in (select sample_id from sample where pub_id='{}' and (data->>'replicate' is null or  data->>'replicate' in ('1', '-')))".format(pub_id)
                else:
                    query += " \nWHERE sample_id in (select sample_id from sample where pub_id='{}')".format(
                        pub_id)

        con = self.connection or self._connect()
        #print('Querying database\n')
        dataframe = read_sql(query, con)
        if self.connection is None:
            con.close()

        if len(dataframe) == 0:
            print('No entries found')
            print(query)
            return dataframe
        if 'data' in dataframe:
            data_extr = dataframe['data'].values

            for key in data_extr[0].keys():
                dataframe[key] = [row.get(key, None) for row in data_extr]

            dataframe = dataframe.drop(axis=1, labels='data')

        return dataframe

    def delete(self, pub_id):
        con = self.connection or self._connect()
        cur = con.cursor()
        self._initialize(con)

        cur.execute(
            """SELECT id from publication where pub_id='{pub_id}'"""
            .format(pub_id=pub_id))
        id_p = cur.fetchone()

        if id_p is None:
            print('Dataset not found')
            return

        """Delete echemical"""
        cur.execute(
            """DELETE from echemical where sample_id in
        (select sample_id from sample where pub_id='{pub_id}');"""
            .format(pub_id=pub_id))

        """Delete sample"""
        cur.execute(
            """DELETE from sample where pub_id='{pub_id}';"""
            .format(pub_id=pub_id))
        """Delete xps"""
        cur.execute(
            """DELETE from xps where mat_id in
            (select mat_id from material where pub_id='{pub_id}');"""
            .format(pub_id=pub_id))
        """Delete xrd"""
        cur.execute(
            """DELETE from xrd where mat_id in
            (select mat_id from material where pub_id='{pub_id}');"""
            .format(pub_id=pub_id))
        """Delete material"""
        cur.execute(
            """DELETE from material where pub_id='{pub_id}';"""
            .format(pub_id=pub_id))
        """Delete publication"""
        cur.execute(
            """DELETE from publication where pub_id='{pub_id}';"""
            .format(pub_id=pub_id))
        print('Delete complete')
        con.commit()
        con.close()

    def write_from_datasheet(self, directory, exclude_dois=[]):
        dataframes = get_dataframes_from_sheet(directory)

        server_dois = self.get_dataframe(table='publication')['doi'].values
        dois = set(dataframes['Tabulated data']['DOI'].values)
        for doi in dois:
            if not '10.' in str(doi):
                continue
            if doi in server_dois or doi in exclude_dois:
                continue
            print('Adding publication with doi={} to server'.format(doi))
            self.write(dataframes, doi=doi)

    def write(self, dataframes, doi=None, dataset=None, pub_info=None):
        con = self.connection or self._connect()
        cur = con.cursor()
        self._initialize(con)

        main_table = dataframes['Tabulated data']
        if doi:
            main_table = main_table[main_table['DOI'] == doi]
        elif dataset:
            main_table = main_table[main_table['dataset_name'] == dataset]

        XPS_table = dataframes['XPS']
        XPS_HR_table = dataframes['XPS_HR']
        XRD_table = dataframes['XRD']
        CV_table = dataframes['CV']
        CVend_table = dataframes['CVend']
        ###Stability_table = dataframes['Stability Test']

        if not pub_info and doi is not None:
            pub_info = doi_request(doi)
        elif doi is None:
            assert pub_info is not None


        # np.unique(main_table['inputer_name'].values).tolist()}}
        tags = {'experimental': {'inputer_names': list(
            set([n for n in main_table['inputer_name'].values if not pandas.isnull(n)]))}}

        dataset_names = list(
            set([n for n in main_table['dataset_name'].values if not pandas.isnull(n)]))
        inputer_names = list(
            set([n for n in main_table['inputer_name'].values if not pandas.isnull(n)]))
        if dataset_names:
            tags['dataset_names'] = dataset_names

        if inputer_names:
            tags['inputer_names'] = inputer_names

        pub_info['tags'] = json.dumps(tags)

        """Write publication """
        pub_id = pub_info['pub_id']
        cur.execute(
            """SELECT id from publication where pub_id='{pub_id}'"""
            .format(pub_id=pub_id))
        id_p = cur.fetchone()

        #if id_p is not None:
        #    print('Dataset already exists: skipping!')
        #    return

        key_list = get_key_list('publication', start_index=1)
        pub_values = [pub_info[k] for k in key_list]
        key_str = get_key_str('publication', start_index=1)
        value_str = get_value_str(pub_values, start_index=0)

        insert_command = """INSERT INTO publication ({0}) VALUES
        ({1});""".format(key_str, value_str)
        print(insert_command)
        cur.execute(insert_command)

        main_table['ICDD_ID'] = main_table['ICDD_ID'].apply(
            lambda x: str(x).split('/'))

        main_table['ICSD_ID'] = main_table['ICSD_ID'].apply(lambda x: [int(float(i)) for i in str(
            x).split('/') if not (pandas.isnull(i) or i == 'nan' or i == '-')])

        XRD_ids = []
        mat_ids = {}
        for index, row in main_table.iterrows():
            row['pub_id'] = pub_id
            XPS_id = row['XPS_ID']
            XRD_id = row['XRD(CuKa)ID']
            CV_init_id = row['CV_intial_ID']
            CV_end_id = row['CV_end_ID']
            print(CV_init_id, CV_end_id)
            Stab_id = row['stability_test_ID']
            XPS_post_id = row['XPS_post_test_ID']

            if pandas.isnull(CV_init_id) and pandas.isnull(CV_end_id):
                print('No CV for row {}. Skipping!'.format(index))
                continue

            if XRD_id in XRD_ids: # already written
                mat_id = mat_ids[XRD_id]
            else:
                if not pandas.isnull(XRD_id):
                    XRD_ids += [XRD_id]
                # Material
                mat_value_list = [row[c] for c in mat_columns]
                # mat_value_list[4] =
                # mat_value_list[4] = [m.replace('-', '') for m in str(
                #    mat_value_list[4]).split('/')]
                key_str = ', '.join(mat_columns)
                key_str = key_str.replace(
                    '_a(nm)', '').replace('ICSD_ID', 'icsd_ids').replace('ICDD_ID', 'icdd_ids')

                value_str = get_value_str(mat_value_list)

                query_mat = \
                    """INSERT INTO material ({}) VALUES ({}) RETURNING mat_id
                    """.format(key_str, value_str)
                print('Wrote material')
                cur.execute(query_mat)
                mat_id = cur.fetchone()[0]
                if not pandas.isnull(XRD_id):
                    mat_ids[XRD_id] = mat_id
                if XPS_id in XPS_table:
                    # XPS
                    E = XPS_table[XPS_id].values[1:]
                    idx = [not pandas.isnull(e) for e in E]
                    E = np.array(E)[idx]
                    I = np.array(XPS_table[XPS_id + '.1'].values[1:])[idx]

                    key_str = ', '.join(xps_columns)
                    xps_value_list = [mat_id, None, 'survey', list(E), list(I)]
                    value_str = get_value_str(xps_value_list)

                    query_xps = \
                        """INSERT INTO xps ({}) VALUES ({})
                        """.format(key_str, value_str)
                    cur.execute(query_xps)
                    print('Wrote survey XPS')
                    xps_hr_elements = row['XPS_spectra']
                    if not pandas.isnull(xps_hr_elements):

                        xps_hr_elements = xps_hr_elements.lstrip('[').rstrip(']').replace(' ', '').split(';')
                        for e in xps_hr_elements:
                            XPS_ID_hr = '{}_{}'.format(XPS_id, e)

                            if XPS_ID_hr in XPS_HR_table:
                                E = XPS_HR_table[XPS_ID_hr].values[1:]
                                idx = [not pandas.isnull(e) for e in E]
                                E = np.array(E)[idx]
                                I = np.array(XPS_HR_table[XPS_ID_hr + '.1'].values[1:])[idx]

                                key_str = ', '.join(xps_columns)
                                xps_value_list = [mat_id, None, 'highres', list(E), list(I)]
                                value_str = get_value_str(xps_value_list)

                                query_xps = \
                                    """INSERT INTO xps ({}) VALUES ({})
                                    """.format(key_str, value_str)
                                cur.execute(query_xps)
                                print('Wrote XPS')

                if XRD_id in XRD_table:
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
                    print('wrote XRD')

            # Sample table
            pub_id = row.pop('pub_id')
            sample_dict = row.to_dict()
            clean_dict = {}
            for k, v in sample_dict.items():
                if isinstance(v, list):
                    continue
                if pandas.isnull(v) or v == '' or v == 'nan':
                    v = None
                if k == 'ICSD ID':
                    v = [int(m) for m in v.split('/')]
                if k == 'ICDD ID':
                    v = [m for m in v.split('/')]

                # k = k.replace('(', '\(')#.replace(')', '').replace('=', '-').replace(',', '_')#.replace(
                # '/', '_o_').replace('%', 'perc').replace('+', 'plus').replace('-', 'minus').replace('.', 'd')
                clean_dict[k] = v


            key_str = ', '.join(sample_columns)
            value_str = get_value_str([mat_id, pub_id, json.dumps(clean_dict)])

            query_sample = \
                """INSERT INTO sample ({}) VALUES ({}) RETURNING sample_id
                """.format(key_str, value_str)
            print('wrote sample')
            cur.execute(query_sample)
            sample_id = cur.fetchone()[0]

            """ Echemical table"""
            key_str = ', '.join(echemical_columns)
            if not pandas.isnull(CV_init_id):
                if not CV_init_id in CV_table:
                    print('ID not found cvinit!', CV_init_id)
                else:
                    V = CV_table[CV_init_id + '.1'].values[1:]
                    #print('shape', V.shape)
                    #if len(V.shape > 1):
                #        V = V[:,0]
                    idx = [not (pandas.isnull(v) or v == '--') for v in V]
                    V = np.array(V)[idx]
                    I = np.array(CV_table[CV_init_id + '.2'].values[1:])[idx]

                    cv_value_list = [sample_id, 'CV_initial',
                                     row['time_tested(h)'], None, list(V), list(I)]
                    value_str = get_value_str(cv_value_list)

                    query_cv = \
                        """INSERT INTO echemical ({}) VALUES ({})
                        """.format(key_str, value_str)
                    print('wrote CV')
                    cur.execute(query_cv)

            if not pandas.isnull(CV_end_id):
                if not CV_end_id in CVend_table:
                    print('Id not found cvend!', CV_end_id)
                else:
                    V = CVend_table[CV_end_id + '.1'].values[1:]
                    idx = [not (pandas.isnull(v) or v == '--') for v in V]
                    V = np.array(V)[idx]
                    I = np.array(CVend_table[CV_end_id + '.2'].values[1:])[idx]

                    cv_value_list = [sample_id, 'CV_end',
                                     row['time_tested(h)'], None, list(V), list(I)]
                    value_str = get_value_str(cv_value_list)

                    query_cv = \
                        """INSERT INTO echemical ({}) VALUES ({})
                        """.format(key_str, value_str)
                    print('wrote CV end')
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

        return

def load_table(filename, tablename):
    skiprows = {'Tabulated data': [0, 1], 'Tabulated Data': [0, 1]}
    data = {}

    data = pandas.read_excel(filename, sheet_name=tablename,
                             skiprows=skiprows.get(tablename),
                             mangle_dupe_cols=True)

    return data


def get_dataframes_from_sheet(directory, use_cache=True):
    filenames = {
        '1-OxR _HxR ActivityDatabase - 1.xlsx':
        [
        'Tabulated data',
        'pretest CV (t-i-V)',
        'posttest CV (t-i-V)',
        'stability (t-i-V)'
        ],
        '2- OxR _HxR ActivityDatabase - 2.xlsx':
        [
         'Tabulated data',
         'pretest CV (t-i-V)',
         'posttest CV (t-i-V)',
         'stability (t-i-V)'
         ],
        'XPS- OxR _HxR ActivityDatabase Survey.xlsx': ['XPS', 'XPSPosttest'],
        'HR XPS- OxR _HxR ActivityDatabase_carbon shifted .xlsx': ['HR data and background'],
        'XRD- OxR _HxR ActivityDatabase .xlsx': ['Sheet1']
    }
    name_map = {'pretest CV (t-i-V)': 'CV',
                'posttest CV (t-i-V)': 'CVend',
                'stability (t-i-V)': 'Stability Test',
                'Sheet1': 'XRD',
                'HR data and background': 'XPS_HR',
                #'Tabulated Data': 'Tabulated data'
                }
    dataframes = {}
    for filename, tables in filenames.items():
        for tableo in tables:
            print(tableo)
            table = name_map.get(tableo, tableo)
            pklname = '{}/{}.pkl'.format(directory, table)
            if use_cache and os.path.isfile(pklname):
                data = pandas.read_pickle(pklname)
                dataframes[table] = data
                continue

            data = load_table(filename=directory+filename, tablename=tableo)
            data = clean_column_names(data)

            if table in list(dataframes.keys()): ## append
                if table == 'Tabulated data':
                    axis = 0
                else:
                    axis = 1
                data = pandas.concat([dataframes[table], data], sort=False, axis=axis)
            data = data.loc[:,~data.columns.duplicated()]
            data = data.replace('-', np.nan)
            dataframes[table] = data
            if not filename == '1-OxR _HxR ActivityDatabase - 1.xlsx':
                data.to_pickle(pklname)
    return dataframes


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

    values = np.array(dataframe.values[1:, :])
    duplicates_dict = {}

    for i in range(len(values)):
        duplicate_rows = []
        for v in duplicates_dict.values():
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

    dataframe = dataframe.drop(all_duplicates)

    return dataframe, duplicates_dict
