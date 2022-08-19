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

import bokeh
from bokeh.plotting import figure, show
from bokeh.layouts import column, row, layout
from bokeh.models import Div, ColumnDataSource, DataTable, DateFormatter, TableColumn, HTMLTemplateFormatter
from bokeh.palettes import Bokeh, Spectral, RdYlBu
palette = Bokeh[8] + Spectral[11][::-1] + RdYlBu[11]

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

    """CREATE TABLE xas_ledge (
    mat_id integer PRIMARY KEY REFERENCES material (mat_id),
    type text,
    energy DOUBLE PRECISION[],
    intensity DOUBLE PRECISION[]
    );""",

    """CREATE TABLE xas_kedge (
    mat_id integer PRIMARY KEY REFERENCES material (mat_id),
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

    def write_from_datasheet(self, directory):
        dataframes = get_dataframes_from_sheet(directory)

        server_dois = self.get_dataframe(table='publication')['doi'].values
        dois = set(dataframes['Tabulated data']['DOI'].values)
        for doi in dois:
            if not '10.' in str(doi):
                continue
            if doi in server_dois:
                continue
            print('Adding publication with doi={} to server'.format(doi))
            self.write(dataframes, doi=doi)

    def write(self, dataframes, doi=None):
        con = self.connection or self._connect()
        cur = con.cursor()
        self._initialize(con)

        main_table = dataframes['Tabulated data']
        main_table = main_table[main_table['DOI'] == doi]
        XPS_table = dataframes['XPS']
        XRD_table = dataframes['XRD']
        CV_table = dataframes['CV']
        CVend_table = dataframes['CVend']
        #Stability_table = dataframes['Stability Test']

        pub_info = doi_request(doi)
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

        if id_p is not None:
            print('Dataset already exists: skipping!')
            return

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
        main_table['ICSD_ID'] = main_table['ICSD_ID'].apply(lambda x: [int(i) for i in str(
            x).split('/') if not (pandas.isnull(i) or i == 'nan' or i == '-')])

        XPS_ids = []
        for index, row in main_table.iterrows():
            row['pub_id'] = pub_id
            XPS_id = row['XPS_ID']
            XRD_id = row['XRD(CuKa)ID']
            CV_init_id = row['CV_intial_ID']
            CV_end_id = row['CV_end_ID']
            Stab_id = row['stability_test_ID']
            XPS_post_id = row['XPS_post_test_ID']
            if not XPS_id in XPS_ids:
                if XPS_id is not None:
                    XPS_ids += [XPS_id]

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
                print(query_mat)
                cur.execute(query_mat)
                mat_id = cur.fetchone()[0]
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

            cur.execute(query_sample)
            sample_id = cur.fetchone()[0]

            """ Echemical table"""
            key_str = ', '.join(echemical_columns)
            if not pandas.isnull(CV_init_id):
                if not CV_init_id in CV_table:
                    print('ID not found cvinit!', CV_init_id)
                else:
                    V = CV_table[CV_init_id + '.1'].values[1:]
                    idx = [not (pandas.isnull(v) or v == '--') for v in V]
                    V = np.array(V)[idx]
                    I = np.array(CV_table[CV_init_id + '.2'].values[1:])[idx]

                    cv_value_list = [sample_id, 'CV_initial',
                                     row['time_tested(h)'], None, list(V), list(I)]
                    value_str = get_value_str(cv_value_list)

                    query_cv = \
                        """INSERT INTO echemical ({}) VALUES ({})
                        """.format(key_str, value_str)
                    print(query_cv)
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

        return

    def show_publications(self):

        all_publications = self.get_dataframe(table='publication')[
            ['pub_id', 'doi', 'title', 'authors', 'journal', 'year']]

        columns = all_publications.columns.values
        prettycolumns = [''.join([s.capitalize()
                                 for s in c.split('_')]) for c in columns]
        source = ColumnDataSource(all_publications)
        Columns = [
            TableColumn(field=c, title='DOI', formatter=HTMLTemplateFormatter(
                template='<a href=" https://doi.org/<%= doi %>" target="_blank"><%= value %></a>'))
            if c == 'doi' else TableColumn(field=c, title=prettycolumns[i], formatter=HTMLTemplateFormatter(template='<%= value %>')) for i, c in enumerate(columns)]

        data_table = DataTable(
            source=source, columns=Columns, width=1600, selectable=True, editable=True)

        logo = Div(text="""<img src="https://www.catalysis-hub.org/329d0ae97b6a9b2266d9cc3d41fffcc1.png" alt="div_image" width=500>""",
                   height=180, width=500)

        div = Div(
            text="""
            Welcome to the Catalysis-Hub Experimental Database!
            """,
            style={'font-size': '200%', 'color': '#8C1515', 'padding': '0px'},
            width=500)

        div2 = Div(
            text="""
            Available datasets are listed in the table below. Query a dataset from the command line interface with <br> <code>$ cathub exp PubId </code>
            """,
            style={'font-size': '120%',
                   'padding': '0px'},
            width=400)

        l = layout([
            row(column(div, div2), logo),  # [div],logo],
            # [div2],
            [data_table],
        ]
        )
        show(l)

    def show_dataset(self, pub_id):

        dataframe_publication = self.get_dataframe(
            'publication', pub_id=pub_id)

        dataframe_sample = self.get_dataframe('sample', pub_id=pub_id)
        if len(dataframe_sample) == 0:
            print('No data')
            return

        dataframe_material = self.get_dataframe(
            'material', pub_id=pub_id).replace([None], [np.nan])
        dataframe_material = dataframe_material.dropna(
            axis=1, thresh=len(dataframe_material))

        dataframe_xps = self.get_dataframe('xps', pub_id=pub_id)
        dataframe_xps = dataframe_xps.set_index('mat_id').join(
            dataframe_material.set_index('mat_id'))
        dataframe_xps = dataframe_xps.sort_values(by=['composition'])

        dataframe_xrd = self.get_dataframe('xrd', pub_id=pub_id)
        dataframe_xrd = dataframe_xrd.set_index('mat_id').join(
            dataframe_material.set_index('mat_id'))
        dataframe_xrd = dataframe_xrd.sort_values(by=['composition'])

        dataframe_cv = self.get_dataframe('echemical', pub_id=pub_id)
        dataframe_cv = dataframe_cv.set_index('sample_id').join(
            dataframe_sample.set_index('sample_id'))
        dataframe_cv = dataframe_cv.sort_values(by=['composition'])

        all_plots = {}

        columns = list(dataframe_material.columns.values)
        columns.remove('mat_id')
        columns.remove('pub_id')
        columns.remove('icdd_ids')
        prettycolumns = [' '.join([s for s in c.split('_')]) for c in columns]
        source = ColumnDataSource(dataframe_material)
        Columns = [
            TableColumn(field=c, title=prettycolumns[i]) for i, c in enumerate(columns)]  # , formatter=HTMLTemplateFormatter(template='<%= value %>')) for i,c in enumerate(columns)]

        data_table_mat = DataTable(
            source=source, columns=Columns, selectable=True, editable=True)

        doi = dataframe_sample['DOI'].values[0]
        title = dataframe_publication['title'].values[0]
        year = dataframe_publication['year'].values[0]
        journal = dataframe_publication['journal'].values[0]
        authors = ','.join([' '.join(d.split(',')[::-1])
                           for d in dataframe_publication['authors'].values[0]])

        div = Div(
            text="""
        <p><b>Viewing dataset: {} </b></p>
        <i>"{}"</i>. {}. {} ({}). <a href=" https://doi.org/{}" target="_blank">{}</a>
        """.format(pub_id, title, authors, journal, year, doi, doi),
            style={'font-size': '120%'}, width=600)

        logo = Div(text="""<img src="https://www.catalysis-hub.org/329d0ae97b6a9b2266d9cc3d41fffcc1.png" alt="div_image" width=500>""",
                   height=180, width=500)

        material_header = Div(
            text=""" Materials """, style={'font-size': '200%'})
        char_header = Div(
            text="""Characterization """, style={'font-size': '200%'})
        echem_header = Div(
            text="""Electrochemical testing """, style={'font-size': '200%'})
        cv_header = Div(
            text="""CV curves """, style={'font-size': '200%'})

        p1 = plot_overpotential(dataframe_sample)
        p2 = plot_cvs(dataframe_cv, cv_type='initial')
        p3 = plot_cvs(dataframe_cv, cv_type='end')

        p4 = plot_xps(dataframe_xps)
        p5 = plot_xrd(dataframe_xrd)

        cv_plots = []
        for cv in [p2, p3]:
            if cv is None:
                continue
            cv_plots += [cv]

        char_plots = []

        for c in [p4, p5]:
            if c is None:
                continue
            char_plots += [c]

        l = layout([
            [div, logo],
            [column(echem_header, p1), column(
                material_header, data_table_mat)],
            [char_header],
            char_plots,
            [cv_header],
            cv_plots
        ])
        show(l)


def load_table(filename, tablename):
    skiprows = {'Tabulated data': [0, 1]}
    data = {}

    data = pandas.read_excel(filename, sheet_name=tablename,
                             skiprows=skiprows.get(tablename),
                             mangle_dupe_cols=True)

    return data


def get_dataframes_from_sheet(directory):
    filenames = {
        '1-OxR _HxR ActivityDatabase - 1.xlsx':
        ['Tabulated data',
         'pretest CV (t-i-V)',
         'posttest CV (t-i-V)',
         'stability (t-i-V)'],
        'XPS- OxR _HxR ActivityDatabase_.xlsx': ['XPS'],  # , 'XPSPosttest'],
        'XRD- OxR _HxR ActivityDatabase_.xlsx': ['Sheet1']
    }
    name_map = {'pretest CV (t-i-V)': 'CV',
                'posttest CV (t-i-V)': 'CVend',
                'stability (t-i-V)': 'Stability Test',
                'Sheet1': 'XRD'}
    dataframes = {}
    for filename, tables in filenames.items():
        for table in tables:
            data = load_table(filename=directory+filename, tablename=table)
            data = clean_column_names(data)
            table = name_map.get(table, table)
            dataframes[table] = data
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


def plot_overpotential(dataframe, currents=[0.01, 0.05, 0.1, 1, 10]):

    source = ColumnDataSource(dataframe)
    p = figure(title='{} performance'.format('+'.join(set(dataframe['reaction'].values))),
               x_axis_label='Material',
               y_axis_label='Onset potential(V)')

    V = {}
    idx = None
    NV = []

    for I in currents:
        V[I] = np.array([v if not v == '-' else None for v in
                         dataframe['onset_potential(+/-{}_mA/cm2)'.format(I)].values])
        NV += [len([v for v in V[I] if v is not None])]
    if len(set(NV)) > 1:
        I_sort = np.argsort(NV)[::-1]
    else:
        I_sort = range(len(currents))
    for i, I in enumerate([currents[i] for i in I_sort]):
        if np.all([i is None for i in V[I]]):
            continue
        if idx is None:
            sort_nonone = [v if v is not None else 100 for v in V[I]]
            idx = np.argsort(sort_nonone)
        p.line(range(len(idx)), V[I][idx], color=palette[i],
               legend_label='{} mA/cm2'.format(I))
        p.circle(range(len(idx)), V[I][idx], color=palette[i])
        #p.plot(range(len(idx)), V[I][idx], '.-', markersize=12, linewidth=2, label = '{} mA/cm2'.format(I))

    p.xaxis.ticker = list(range(len(V[I])))
    materials_list = [mat if not mat is None else 'Support [{}]'
                      .format(dataframe['conductive_support_ID'].values[0])
                      for mat in dataframe['composition'].values[idx]]
    tick_dct = {}
    for i, v in enumerate(list(range(len(V[I])))):
        tick_dct[v] = materials_list[i]

    p.xaxis.major_label_overrides = tick_dct
    p.xaxis.major_label_orientation = 3.1415/3
    p.legend.location = 'top_left'

    return p


def plot_cvs(dataframe, cv_type='initial', cv_ids=None):
    dataframe = dataframe[dataframe['type'] == 'CV_{}'.format(cv_type)]
    #source = ColumnDataSource(dataframe)
    p = figure(title='CV ({})'.format(cv_type),
               x_axis_label='Potential (V)',
               y_axis_label='Current (mA/cm2)')

    if len(dataframe) == 0:
        return None

    if cv_ids is None:
        cv_ids = range(len(dataframe))
    for ii, i in enumerate(cv_ids):
        material = '{}({})'.format(dataframe['composition'].values[i],
                                   dataframe['id'].values[i])
        p.circle(dataframe['potential'].values[i], dataframe['current'].values[i],
                 color=palette[ii], legend_label=material)
        # p.plot(dataframe['potential'].values[i], dataframe['current'].values[i],
        #    linewidth=2, label=material)
    p.legend.location = 'top_left'
    return p


def plot_xps(dataframe, type='survey', xps_ids=None):

    pub_id = set(dataframe['pub_id'].values)
    dataframe = dataframe[dataframe['type'] == '{}'.format(type)]

    p = figure(title='XPS ({})'.format(type),
               x_axis_label='Binding Energy (eV)',
               y_axis_label='Intensity')

    if len(dataframe) == 0:
        return None

    if xps_ids is None:
        xps_ids = range(len(dataframe))
    for ii, i in enumerate(xps_ids):
        p.circle(dataframe['binding_energy'].values[i], dataframe['intensity'].values[i],
                 color=palette[ii], legend_label=dataframe['composition'].values[i])

    p.legend.location = 'top_left'

    return p


def plot_xrd(dataframe, xrd_ids=None):

    pub_id = set(dataframe['pub_id'].values)
    #dataframe = dataframe[dataframe['type']== '{}'.format('pre')]

    p = figure(title='XRD',
               x_axis_label='2 theta (degree)',
               y_axis_label='Intensity')

    if len(dataframe) == 0:
        return None

    if xrd_ids is None:
        xrd_ids = range(len(dataframe))
    for ii, i in enumerate(xrd_ids):
        p.circle(dataframe['degree'].values[i], dataframe['intensity'].values[i],
                 color=palette[ii], legend_label=dataframe['composition'].values[i])
    p.legend.location = 'top_left'

    return p


def get_publication_label(dataframe):
    title = dataframe['title'].values[0]
    authors = dataframe['authors'].values[0]


if __name__ == '__main__':
    import sys
    DB = ExpSQL(user='expvisitor',
                password=os.environ.get('DB_PASSWORD_expvis'))
    # DB.show_publications()
    DB.show_dataset('HubertAcidic2020')
