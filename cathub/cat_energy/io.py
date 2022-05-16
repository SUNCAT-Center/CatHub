'''
Module to perform input/output operations
'''

from pathlib import Path

import pandas as pd
from sqlalchemy import create_engine


NUM_DECIMAL_PLACES = 4
write_columns = ['surface_name', 'site_name', 'species_name', 'raw_energy',
                 'elec_energy', 'dft_corr', 'zpe', 'enthalpy', 'entropy',
                 'rhe_corr', 'solv_corr', 'formation_energy', 'energy_vector',
                 'frequencies', 'reference']


# NOTE: use cathub function to convert db to dataframe
def db_to_dataframe(table_name, filename):
    '''
    Read cathub .db file into pandas dataframe
    '''

    # define sql url
    sql_url = 'sqlite:///' + str(filename)

    # SQLAlchemy connectable
    cnx = create_engine(sql_url).connect()

    # table will be returned as a dataframe
    df = pd.read_sql_table(table_name, cnx)
    return df

def read_qe_log(file_path, wf_dipole_index):
    '''
    read quantum espresso log file and return converged final energy
    and workfunction
    '''
    wf_line_index = -1
    energy_line_index = -1
    with open(file_path, encoding='utf8') as file:
        for line_index, line in enumerate(file.readlines()):
            if 'wf' in line:
                wf_line_index = line_index + 1
                energy_line_index = line_index + 3
            if line_index == wf_line_index:
                if wf_dipole_index == 0:
                    workfunction = float(line.split('[')[1].split(',')[0])
                elif wf_dipole_index == 1:
                    workfunction = float(line.split(']')[0].split(' ')[1])
            if line_index == energy_line_index:
                final_energy = float(line.split()[0])
    return (final_energy, workfunction)

def make_mkm_input_files(db_filepath, adsorbate_parameters, field_effects,
                         df_out):
    '''
    Function to write input files for mkm simulations using CatMAP
    '''
    system_dir_path = (db_filepath.parent
                       / f'{adsorbate_parameters["desired_surface"]}'
                       f'_{adsorbate_parameters["desired_facet"]}')
    Path.mkdir(system_dir_path, parents=True, exist_ok=True)
    energies_filepath = (system_dir_path
                         / f'energies_f{field_effects["epsilon"]:.2e}.txt')

    header = '\t'.join(['surface_name', 'site_name', 'species_name',
                        'formation_energy', 'frequencies', 'reference'])
    lines = [] # list of lines in the output
    for _, row in df_out.iterrows():
        line = '\t'.join([row['surface_name'], row['site_name'],
                          row['species_name'], f'{row["formation_energy"]:.4f}',
                          str(row['frequencies']), row['reference']])
        lines.append(line)

    lines = [header] + lines #add header to top
    input_file = '\n'.join(lines) #Join the lines with a line break

    with open(energies_filepath, 'w', encoding='utf8') as energies_file:
        energies_file.write(input_file)
    return None
