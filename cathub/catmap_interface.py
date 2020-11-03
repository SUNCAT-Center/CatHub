import re

import pandas as pd
from sqlalchemy import create_engine 
from ase.db import connect
from ase.data import chemical_symbols, atomic_numbers
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

def write_energies(db_filepath, critical_density, reference_gases, dummy_gases, dft_corrections, offset, field_effects, adsorbate_parameters, write_gases, write_adsorbates):

    df_out = pd.DataFrame(columns=['Surface Name', 'Site Name', 'Species Name', 'Formation Energy'])

    if write_gases:
        df_out = write_gas_energies(db_filepath, df_out, critical_density, reference_gases, dummy_gases, dft_corrections, offset, field_effects)

    if write_adsorbates:
        df_out = write_adsorbate_energies(db_filepath, df_out, adsorbate_parameters)

    # write corrected energy data to file
    energies_filepath = db_filepath.parent / f'energies_f{field_effects["epsilon"]:.2e}.txt'
    with open(energies_filepath, 'w') as energies_file:
        df_out.to_string(energies_file, index=False)
    return None

def write_gas_energies(db_filepath, df_out, critical_density, reference_gases, dummy_gases, dft_corrections, offset, field_effects):
    "Write formation energies to energies.txt after applying free energy corrections"

    # identify system ids for gaseous species
    table_name = 'systems'
    df1 = db_to_dataframe(table_name, str(db_filepath))
    gas_ids = list(df1.id[df1.mass / df1.volume < critical_density])
    num_decimal_places = 5

    # record energies for reference gases
    db = connect(str(db_filepath))
    gas_select_rows = [list(db.select(id=gas_id))[0] for gas_id in gas_ids]
    surface, site, species, formation_energies = [], [], [], []
    reference_gas_energies = {}
    for row in gas_select_rows:
        if row.formula in reference_gases:
            reference_gas_energies[row.formula] = row.energy + dft_corrections[row.formula]

    # build dataframe data for dummy gases
    dummy_gas_energy = 0.0
    for dummy_gas in dummy_gases:
        surface.append('None')
        site.append('gas')
        species.append(dummy_gas)
        formation_energies.append(f'{dummy_gas_energy:.{num_decimal_places}f}')

    # build dataframe data for gaseous species
    for row in gas_select_rows:
        mass = row.mass
        volume = row.volume
        surface.append('None')
        site.append('gas')
        if row.formula in reference_gases:
            relative_energy = 0.0
        else:
            chemical_symbols_dict = formula_to_chemical_symbols(row.formula)
            for chemical_symbol in chemical_symbols_dict.keys():
                count = chemical_symbols_dict[chemical_symbol]

            # xCO + (x-z+y/2)H2 --> CxHyOz + (x-z)H2O
            if 'C' in chemical_symbols_dict:
                x = chemical_symbols_dict['C']
            else:
                x = 0
            if 'H' in chemical_symbols_dict:
                y = chemical_symbols_dict['H']
            else:
                y = 0
            if 'O' in chemical_symbols_dict:
                z = chemical_symbols_dict['O']
            else:
                z = 0
            relative_energy = (row.energy
                               + (x - z) * reference_gas_energies['H2O']
                               - x * reference_gas_energies['CO']
                               - (x - z + y / 2) * reference_gas_energies['H2'])

        # Apply offset
        if row.formula in offset:
            relative_energy += offset[row.formula]

        # Apply field effects
        epsilon = field_effects['epsilon']
        pH = field_effects['pH']
        d = field_effects['d']
        n = field_effects['n']
        UM_PZC = field_effects['UM_PZC']
        mu = field_effects['mu']
        alpha = field_effects['alpha']

        # U_RHE-scale field effect
        U_SHE = epsilon * d + UM_PZC
        U_RHE = U_SHE + 0.059 * pH
        if row.formula in n:
            relative_energy += n[row.formula] * U_RHE

        # U_SHE-scale field effects
        if row.formula in mu:
            relative_energy += (mu[row.formula] * epsilon
                                - 0.5 * alpha[row.formula] * epsilon**2)

        species.append(row.formula)
        formation_energies.append(f'{relative_energy:.{num_decimal_places}f}')

    df2 = pd.DataFrame(list(zip(surface, site, species, formation_energies)),
                       columns=['Surface Name', 'Site Name', 'Species Name', 'Formation Energy'])
    df_out = df_out.append(df2)
    return df_out

def write_adsorbate_energies(db_filepath, df_out, adsorbate_parameters):
    "Write formation energies to energies.txt after applying free energy corrections"

    # identify system ids for adsorbate species
    table_name = 'reaction'
    df1 = db_to_dataframe(table_name, str(db_filepath))

    desired_surface = adsorbate_parameters['desired_surface']
    desired_facet = adsorbate_parameters['desired_facet']
    df2 = df1.loc[df1['surface_composition'] == desired_surface]
    df1 = df1.loc[df1['facet'].str.contains(desired_facet)]

    df2 = pd.DataFrame(list(zip(surface, site, species, formation_energies)),
                       columns=['Surface Name', 'Site Name', 'Species Name', 'Formation Energy'])
    df_out = df_out.append(df2)
    return df_out

def formula_to_chemical_symbols(formula):
    "Return dictionary mapping chemical symbols to number of atoms"

    chemical_symbols_dict = {}

    # split chemical formula string into alpha and numeric characters
    regex = re.compile('(\d+|\s+)')
    split_formula = regex.split(formula)
    split_formula_list = []

    # count number of formula units if any
    start_index = 0
    formula_unit_count = 1
    if str.isdigit(split_formula[0]):
        formula_unit_count = int(split_formula[0])
        start_index = 1

    # identify chemical symbol and map to its count
    for string in split_formula[start_index:]:
        if str.isdigit(string):
            chemical_symbols_dict[last_chemical_symbol] = int(string)
        else:
            if len(string) == 0:
                pass
            elif len(string) == 1:
                last_chemical_symbol = string
                chemical_symbols_dict[last_chemical_symbol] = 1
            elif len(string) == 2:
                if string in chemical_symbols:
                    last_chemical_symbol = string
                    chemical_symbols_dict[last_chemical_symbol] = 1
                else:
                    chemical_symbols_dict[string[0]] = 1
                    last_chemical_symbol = string[1]
                    chemical_symbols_dict[last_chemical_symbol] = 1
            elif len(string) == 3:
                if string[0] in chemical_symbols:
                    chemical_symbols_dict[string[0]] = 1
                    last_chemical_symbol = string[1:]
                    chemical_symbols_dict[last_chemical_symbol] = 1
                else:
                    chemical_symbols_dict[string[:2]] = 1
                    last_chemical_symbol = string[2]
                    chemical_symbols_dict[last_chemical_symbol] = 1

    # multiply number of atoms for each chemical symbol with number of formula units
    for key in chemical_symbols_dict.keys():
        chemical_symbols_dict[key] = formula_unit_count * chemical_symbols_dict[key]
    return chemical_symbols_dict

