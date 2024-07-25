'''
Module to perform input/output operations
'''

from pathlib import Path


NUM_DECIMAL_PLACES = 4
write_columns = ['surface_name', 'site_name', 'species_name', 'raw_energy',
                 'elec_energy', 'dft_corr', 'zpe', 'enthalpy', 'entropy',
                 'rhe_corr', 'formation_energy', 'energy_vector',
                 'frequencies', 'reference']


def read_qe_log(file_path, wf_dipole_index):
    '''
    Read Quantum Espresso log file and return converged final energy and workfunction.

    This function parses a Quantum Espresso log file to extract the converged 
    final energy and workfunction based on specified indices.

    Parameters:
    -----------
    file_path : pathlib.Path
        Path to the Quantum Espresso log file.
    wf_dipole_index : int
        Index to determine how to extract the workfunction from the log file 
        (0 or 1 based on the specific format in the log file).

    Returns:
    --------
    tuple
        A tuple containing:
        - final_energy (float): The converged final energy.
        - workfunction (float): The workfunction.
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


def make_mkm_input_files(db_filepath, system_parameters, df_out):
    '''
    Write input files for MKM simulations using CatMAP.

    This function generates the necessary input files for microkinetic modeling
    (MKM) simulations using CatMAP, based on the provided system parameters and 
    reaction data. The input file is a tab-separated text file with specific 
    column headers required by CatMAP.

    Parameters:
    -----------
    db_filepath : pathlib.Path
        Path to the ASE database file containing the reaction data.
    system_parameters : dict
        Dictionary containing system parameters such as:
            'desired_surface' (str): The surface material being analyzed.
            'desired_facet' (str): The facet of the surface material being analyzed.
            'temp' (float): Temperature in Kelvin.
            'pH' (float): pH value of the system.
            'u_she' (float): Standard Hydrogen Electrode potential.
    df_out : pandas.DataFrame
        DataFrame containing the computed energetics of the species.

    Returns:
    --------
    None

    Notes:
    ------
    The generated input file will be a tab-separated text file with the following
    required column headers:
        - surface_name
        - site_name
        - species_name
        - formation_energy
        - frequencies
        - reference
    '''
    system_dir_path = (db_filepath.parent
                       / f'{system_parameters["desired_surface"]}'
                       f'_{system_parameters["desired_facet"]}')
    Path.mkdir(system_dir_path, parents=True, exist_ok=True)
    energies_filepath = (system_dir_path
                         / f'energies_she_{system_parameters["u_she"]:.2f}V.txt')

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
