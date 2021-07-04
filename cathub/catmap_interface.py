import re
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sqlalchemy import create_engine 
from ase.db import connect
from ase.data import chemical_symbols, atomic_numbers
from ase.thermochemistry import IdealGasThermo, HarmonicThermo
from ase.units import _e, _hplanck, _c, pi
from tabulate import tabulate


write_columns = ['surface_name', 'site_name', 'species_name', 'raw_energy','elec_energy', 'dft_corr','zpe', 'enthalpy', 'entropy', 'rhe_corr', 'solv_corr', 'formation_energy', 'energy_vector', 'frequencies', 'reference']
mkm_write_columns = ['surface_name', 'site_name', 'species_name', 'formation_energy', 'frequencies', 'reference']
num_decimal_places = 4
CM2M = 1E-02
cm2eV = _hplanck / _e * _c / CM2M
kB = 8.617333262145e-5


def db_to_dataframe(table_name, filename):
    "Read cathub .db file into pandas dataframe"

    # define sql url
    sql_url = 'sqlite:///' + str(filename)

    # SQLAlchemy connectable
    cnx = create_engine(sql_url).connect()

    # table will be returned as a dataframe
    df = pd.read_sql_table(table_name, cnx)
    return df

def write_energies(db_filepath, reference_gases, dummy_gases,
                   dft_corrections_gases, beef_dft_helmholtz_offset,
                   field_effects, adsorbate_parameters, write_gases,
                   write_adsorbates, gas_jsondata_filepath,
                   ads_jsondata_filepath, temp, write_mkm_input_files,
                   verbose=True):

    df_out = pd.DataFrame(columns=write_columns)

    if verbose:
        system_header = f'###### {adsorbate_parameters["desired_surface"]}_{adsorbate_parameters["desired_facet"]}: electric field strength = {field_effects["epsilon"]:.2f} V/A ######'
        print('-' * len(system_header))
        print(system_header)
        print('-' * len(system_header))

        print()
        print('Term1 = Calculated Electronic Energy + DFT Correction')
        print('Term2 = Enthalpy Contribution + Entropy Contribution')
        print('Term3 = RHE-scale Dependency')
        print('Term4 = Solvation Correction + Electric Field Correction')
        print('Chemical Potential, Mu = Term1 + Term2 + Term3 + Term4')
        print('Free Energy Change, G = Mu_species - Mu_ref. For example, G_CH4 = Mu_CH4 - (Mu_CO + 3 * Mu_H2(ref) - Mu_H2O)')
        print()
        
    if write_gases:
        df_out = write_gas_energies(db_filepath, df_out, gas_jsondata_filepath,
                                    reference_gases, dummy_gases,
                                    dft_corrections_gases,
                                    beef_dft_helmholtz_offset, field_effects,
                                    temp, verbose)

    if write_adsorbates:
        df_out = write_adsorbate_energies(db_filepath, df_out,
                                          ads_jsondata_filepath,
                                          adsorbate_parameters, reference_gases,
                                          dft_corrections_gases, field_effects,
                                          temp, verbose)

    # write corrected energy data to mkm input file
    if write_mkm_input_files:
        make_mkm_input_files(db_filepath, adsorbate_parameters, field_effects,
                             df_out)
    return None

def write_gas_energies(db_filepath, df_out, gas_jsondata_filepath,
                       reference_gases, dummy_gases, dft_corrections_gases,
                       beef_dft_helmholtz_offset, field_effects, temp, verbose):

    db = connect(str(db_filepath))
    gas_atoms_rows = list(db.select(state='gas'))

    surface, site, species, raw_energy, elec_energy_calc = [], [], [], [], []
    dft_corr, helm_offset, zpe, enthalpy, entropy = [], [], [], [], []
    rhe_corr, solv_corr, efield_corr, formation_energy = [], [], [], []
    energy_vector, xyz, frequencies, references = [], [], [], []

    # Load vibrational data
    with open(gas_jsondata_filepath) as f:
        gas_data = json.load(f)

    vibrational_energies = {}
    species_list = [species_data['species'] for species_data in gas_data]
    for species_data in gas_data:
        gas_species = species_data['species']
        vibrational_energies[gas_species] = []
        for vibration in species_data['frequencies']:
            vibrational_energies[gas_species].append(vibration * cm2eV)

    reference_gas_energies = {}
    for row in gas_atoms_rows:
        if row.formula in dft_corrections_gases:
            reference_gas_energies[row.formula] = row.energy + dft_corrections_gases[row.formula]
        else:
            reference_gas_energies[row.formula] = row.energy

    # build dataframe data for dummy gases
    dummy_gas_energy = 0.0
    for dummy_gas in dummy_gases:
        surface.append('None')
        site.append('gas')
        species.append(dummy_gas)
        xyz.append([])
        raw_energy.append(0.0)
        elec_energy_calc.append(0.0)
        dft_corr.append(0.0)
        helm_offset.append(0.0)
        zpe.append(0.0)
        enthalpy.append(0.0)
        entropy.append(0.0)
        rhe_corr.append(0.0)
        solv_corr.append(0.0)
        formation_energy.append(0.0)
        energy_vector.append([0.0, 0.0, 0.0, 0.0, 0.0])
        frequencies.append([])
        references.append('')

    # build dataframe data for gaseous species
    for row in gas_atoms_rows:
        species_name = row.formula
        if species_name == 'H2' and 'H2_ref' in species_list:
            species_names = ['H2', 'H2_ref']
        else:
            species_names = [species_name]
        for species_name in species_names:
            surface.append('None')
            site.append('gas')
            species.append(species_name)
    
            chemical_symbols_dict = formula_to_chemical_symbols(species_name.replace('_ref', ''))
            
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
            xyz.append([x, y, z])
            raw_energy.append(row.energy)
            elec_energy_calc.append(row.energy
                                    + (x - z) * reference_gas_energies['H2O']
                                    - x * reference_gas_energies['CO']
                                    - (x - z + y / 2) * reference_gas_energies['H2'])
            dft_corr.append(dft_corrections_gases[species_name] if species_name in dft_corrections_gases else 0.0)
            helm_offset.append(beef_dft_helmholtz_offset[species_name] if species_name in beef_dft_helmholtz_offset else 0.0)
            
            if species_name in species_list:
                species_index = species_list.index(species_name)
                thermo = IdealGasThermo(vib_energies=vibrational_energies[species_name],
                                        geometry=gas_data[species_index]['geometry'],
                                        atoms=row.toatoms(),
                                        symmetrynumber=gas_data[species_index]['symmetry'],
                                        spin=gas_data[species_index]['spin'])
                # zero point energy correction
                zpe.append(thermo.get_ZPE_correction())
                # enthalpy contribution
                enthalpy.append(thermo.get_enthalpy(temp, verbose=False))
                S = thermo.get_entropy(temp, gas_data[species_index]['fugacity'],verbose=False)
                # entropy contribution
                entropy.append(- temp * S)
            else:
                zpe.append(0.0)
                enthalpy.append(0.0)
                entropy.append(0.0)
            # RHE-scale dependency. Zero for gaseous species
            rhe_corr.append(0.0)
            # Solvation correction. Zero for gaseous species.
            solv_corr.append(0.0)
    
            # Apply field effects
            if field_effects:
                (U_RHE_energy_contribution, U_SHE_energy_contribution) = get_electric_field_contribution(field_effects, species_name)
                electric_field_contribution = U_RHE_energy_contribution + U_SHE_energy_contribution
                efield_corr.append(electric_field_contribution)
    
            # compute energy vector
            term1 = elec_energy_calc[-1] + dft_corr[-1]
            term2 = enthalpy[-1] + entropy[-1]
            term3 = rhe_corr[-1]
            term4 = solv_corr[-1] + efield_corr[-1]
            mu = term1 + term2 + term3 + term4
            energy_vector.append([term1, term2, term3, term4, mu])

            # formation energy
            formation_energy.append(term1 + term4 + helm_offset[-1])
            
            if species_name in species_list:
                frequencies.append(gas_data[species_list.index(species_name)]['frequencies'])
                references.append(gas_data[species_list.index(species_name)]['reference'])
            else:
                frequencies.append([])
                references.append('')
    
    reference_mu = {}
    for species_name in reference_gases:
        species_index = species.index(species_name)
        reference_mu[species_name] = energy_vector[species_index][-1]
    
    for species_index, species_name in enumerate(species):
        if species_name in dummy_gases or species_name in reference_gases:
            G = 0.0
        else:
            [x, y, z] = xyz[species_index]
            G = (energy_vector[species_index][-1]
                 + (x - z) * reference_mu['H2O']
                 - x * reference_mu['CO']
                 - (x - z + y / 2) * reference_mu['H2_ref'])
        energy_vector[species_index].append(G)
        
    df = pd.DataFrame(list(zip(surface, site, species, raw_energy,
                               elec_energy_calc, dft_corr, zpe, enthalpy,
                               entropy, rhe_corr, solv_corr, formation_energy,
                               energy_vector, frequencies, references)),
                       columns=write_columns)
    df_out = df_out.append(df, ignore_index=True, sort=False)

    if verbose:
        gas_phase_header = 'Gas Phase Free Energy Correction:'
        print(gas_phase_header)
        print('-' * len(gas_phase_header))
        
        table = []
        table_headers = ["Species", "Term1 (eV)", "Term2 (eV)", "Term3 (eV)", "Term4 (eV)", "mu (eV)", "G (eV)"]
        for index, species_name in enumerate(df['species_name']):
            sub_table = []
            sub_table.extend([species_name,
                              f'{df["energy_vector"][index][0]:.{num_decimal_places}f}',
                              f'{df["energy_vector"][index][1]:.{num_decimal_places}f}',
                              f'{df["energy_vector"][index][2]:.{num_decimal_places}f}',
                              f'{df["energy_vector"][index][3]:.{num_decimal_places}f}',
                              f'{df["energy_vector"][index][4]:.{num_decimal_places}f}',
                              f'{df["energy_vector"][index][5]:.{num_decimal_places}f}'])
            table.append(sub_table)
        print(tabulate(table, headers=table_headers, tablefmt='psql', colalign=("right", ) * len(table_headers), disable_numparse=True))
        print('\n')
    return df_out

def get_electric_field_contribution(field_effects, species_value, reactants=None):
    epsilon = field_effects['epsilon']
    pH = field_effects['pH']
    U_RHE = field_effects['U_RHE']
    mu = field_effects['mu']
    alpha = field_effects['alpha']

    ## U_RHE-scale dependency
    # x CO + y (H++e-) = CxH(y-2x+2z)Oz + (x-z) H2O
    # x CO + y/2 H2 = CxH(y-2x+2z)Oz + (x-z) H2O
    # Based on computational hydrogen electrode, n should be twice the number of H2 gas molecules that are required for the reduction reaction
    if reactants:
        if 'H2gas' in reactants:
            n = 2 * reactants['H2gas']
        else:
            n = 0
    else:
        n = 0

    # # U_RHE derived from epsilon
    # UM_PZC = -0.54  # zero-charge potential of the metal electrode in V vs. SHE for Cu(100)
    # d = 1.2  # thickness in A (angstrom)
    # U_SHE = UM_PZC + d * epsilon
    # U_RHE = U_SHE + 0.059 * pH

    U_RHE_energy_contribution = n * U_RHE

    ## U_SHE-scale dependency
    if reactants == None:
        if species_value + '_g' in mu:
            species_value = species_value + '_g'
    if species_value in mu:
        U_SHE_energy_contribution = mu[species_value] * epsilon - 0.5 * alpha[species_value] * epsilon**2
    else:
        U_SHE_energy_contribution = 0.0

    return (U_RHE_energy_contribution, U_SHE_energy_contribution)

def write_adsorbate_energies(db_filepath, df_out, ads_jsondata_filepath,
                              adsorbate_parameters, reference_gases,
                              dft_corrections_gases, field_effects, temp,
                              verbose):
    "Write formation energies to energies.txt after applying free energy corrections"

    # identify system ids for adsorbate species
    table_name = 'reaction'
    df1 = db_to_dataframe(table_name, str(db_filepath))

    # Load vibrational data
    with open(ads_jsondata_filepath) as f:
        ads_data = json.load(f)

    vibrational_energies = {}
    json_species_list = [species_data['species'] for species_data in ads_data]
    for species_data in ads_data:
        ads_species = species_data['species']
        vibrational_energies[ads_species] = []
        for vibrational_frequency in species_data['frequencies']:
            vibrational_energies[ads_species].append(vibrational_frequency * cm2eV)

    desired_surface = adsorbate_parameters['desired_surface']
    desired_facet = adsorbate_parameters['desired_facet']
    df1 = df1.loc[df1['surface_composition'] == desired_surface]
    df1 = df1.loc[df1['facet'].str.contains(desired_facet)]
    
    ## build dataframe data for adsorbate species
    db = connect(str(db_filepath))
    surface, site, species, raw_energy, facet, elec_energy_calc = [], [], [], [], [], []
    dft_corr, zpe, enthalpy, entropy, rhe_corr = [], [], [], [], []
    solv_corr, formation_energy, efield_corr, energy_vector = [], [], [], []
    frequencies, references = [], []

    # simple reaction species: only one active product and filter out reactions without any adsorbed species
    index_list = []
    for index, product in enumerate(df1['products']):
        if product.count('star') == 1 and 'star' not in json.loads(product):
            index_list.append(index)
    df2 = df1.iloc[index_list]

    products_list = []
    species_list = []
    for index, products_string in enumerate(df2.products):
        products_list.append(json.loads(products_string))
        for product in products_list[-1]:
            if 'star' in product:
                species_list.append(product.replace('star', ''))
    unique_species = sorted(list(set(species_list)), key=len)
    for species_name in unique_species:
        if '-' in desired_surface:
            surface.append(desired_surface.split('-')[0])
        else:
            surface.append(desired_surface)
        site.append(desired_facet)
        species.append(species_name)

        # [adsorption_energy_RHE0, U_RHE_energy_contribution, U_SHE_energy_contribution, solvation_correction]
        (site_wise_energy_contributions, facet_list) = get_adsorption_energies(
                    df2, df_out, species_list, species_name, products_list,
                    reference_gases, dft_corrections_gases, adsorbate_parameters,
                    field_effects)
        site_wise_energy_contributions = np.asarray(site_wise_energy_contributions)
        
        site_wise_adsorption_energies = np.sum(site_wise_energy_contributions, axis=1)
        min_adsorption_energy = min(site_wise_adsorption_energies)
        min_index = np.where(site_wise_adsorption_energies == min_adsorption_energy)[0][0]
        facet.append(facet_list[min_index])
        raw_energy.append(float("nan"))
        elec_energy_calc.append(site_wise_energy_contributions[min_index][0])
        # Zero DFT Correction for Adsorbates
        dft_corr.append(0.0)

        if species_name in json_species_list:
            species_index = json_species_list.index(species_name)
            thermo = HarmonicThermo(vib_energies=vibrational_energies[species_name])
            
            # zero point energy correction
            zpe.append(np.sum(vibrational_energies[species_name]) / 2.0)
            
            # enthalpy contribution
            enthalpy.append(thermo.get_internal_energy(temp,verbose=False))
            
            S = thermo.get_entropy(temp, verbose=False)
            # entropy contribution
            entropy.append(- temp * S)
        else:
            zpe.append(0.0)
            enthalpy.append(0.0)
            entropy.append(0.0)

        rhe_corr.append(site_wise_energy_contributions[min_index][1])
        solv_corr.append(site_wise_energy_contributions[min_index][3])
        
        if field_effects:
            efield_corr.append(site_wise_energy_contributions[min_index][2])

        # compute energy vector
        term1 = elec_energy_calc[-1] + dft_corr[-1]
        term2 = enthalpy[-1] + entropy[-1]
        term3 = rhe_corr[-1]
        term4 = solv_corr[-1] + efield_corr[-1]
        mu = term1 + term2 + term3 + term4
        energy_vector.append([term1, term2, term3, term4, mu])
        formation_energy.append(term1 + term4)
        
        if species_name in json_species_list:
            frequencies.append(ads_data[json_species_list.index(species_name)]['frequencies'])
            references.append(ads_data[json_species_list.index(species_name)]['reference'])
        else:
            frequencies.append([])
            references.append('')

    reference_mu = {}
    for species_name in reference_gases:
        species_index = list(df_out["species_name"]).index(species_name)
        reference_mu[species_name] = df_out["energy_vector"][species_index][-2]
    
    for species_index, species_name in enumerate(species):
        chemical_symbols_dict = formula_to_chemical_symbols(species_name)
        
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

        G = (energy_vector[species_index][-1]
             + (x - z) * reference_mu['H2O']
             - x * reference_mu['CO']
             - (x - z + y / 2) * reference_mu['H2_ref'])
        energy_vector[species_index].append(G)

    df3 = pd.DataFrame(list(zip(surface, site, species, raw_energy,
                                elec_energy_calc, dft_corr, zpe, enthalpy,
                                entropy, rhe_corr, solv_corr, formation_energy,
                                energy_vector, frequencies, references)),
                       columns=write_columns)
    df_out = df_out.append(df3, ignore_index=True, sort=False)

    if verbose:
        adsorbate_phase_header = 'Adsorbate Free Energy Correction:'
        print(adsorbate_phase_header)
        print('-' * len(adsorbate_phase_header))
        
        table = []
        table_headers = ["Species", "Term1 (eV)", "Term2 (eV)", "Term3 (eV)", "Term4 (eV)", "mu (eV)", "G (eV)"]
        for index, species_name in enumerate(df3['species_name']):
            sub_table = []
            sub_table.extend([species_name,
                              f'{df3["energy_vector"][index][0]:.{num_decimal_places}f}',
                              f'{df3["energy_vector"][index][1]:.{num_decimal_places}f}',
                              f'{df3["energy_vector"][index][2]:.{num_decimal_places}f}',
                              f'{df3["energy_vector"][index][3]:.{num_decimal_places}f}',
                              f'{df3["energy_vector"][index][4]:.{num_decimal_places}f}',
                              f'{df3["energy_vector"][index][5]:.{num_decimal_places}f}'])
            table.append(sub_table)
        print(tabulate(table, headers=table_headers, tablefmt='psql', colalign=("right", ) * len(table_headers), disable_numparse=True))
        print('\n')
    return df_out

def get_adsorption_energies(df, df_out, species_list, species_value,
                            products_list, reference_gases, dft_corrections_gases,
                            adsorbate_parameters, field_effects):
    "Compute electronic adsorption energies for a given species at all suitable adsorption sites at a given U_SHE/RHE"
    
    indices = [index for index, value in enumerate(species_list) if value == species_value]
    facet_list = df.facet.iloc[indices].tolist()

    site_wise_energy_contributions = []
    for index, reaction_index in enumerate(indices):
        facet = facet_list[index]
        # NOTE: Reactions with unspecified adsorption site in the facet label are constant-charge NEB calculations and irrelevant for formation_energy calculations.
        # Thus, considering only reactions with specified adsorption site in this code.
        if '-' in facet:
            reactants = json.loads(df.reactants.iloc[reaction_index])
            products = products_list[reaction_index]
            reaction_energy = df.reaction_energy.iloc[reaction_index]
            (adsorption_energy_RHE0,
             U_RHE_energy_contribution,
             U_SHE_energy_contribution,
             solvation_correction) = get_adsorption_energy(
                 df_out, species_value, reactants, products, reaction_energy,
                 reference_gases, dft_corrections_gases, adsorbate_parameters,
                 field_effects)
            site_wise_energy_contributions.append([adsorption_energy_RHE0, U_RHE_energy_contribution, U_SHE_energy_contribution, solvation_correction])
    return (site_wise_energy_contributions, facet_list)

def get_adsorption_energy(df_out, species_value, reactants, products, reaction_energy, reference_gases, dft_corrections_gases, adsorbate_parameters, field_effects):
    "Compute adsorption energy for an adsorbate species in a given reaction"
    
    product_energy = 0
    for product, num_units in products.items():
        if 'star' not in product:
            if 'gas' in product:
                gas_product = product.replace('gas', '')
                if gas_product not in reference_gases:
                    row_index = df_out.index[df_out['species_name'] == gas_product][0]
                    product_energy += float(df_out['formation_energy'].iloc[row_index]) * num_units
                    
                if gas_product in dft_corrections_gases:
                    product_energy += dft_corrections_gases[gas_product] * num_units

    reactant_energy = 0
    for reactant, num_units in reactants.items():
        if 'star' not in reactant:
            if 'gas' in reactant:
                gas_product = reactant.replace('gas', '')
                if gas_product not in reference_gases and (gas_product + '_ref') not in reference_gases:
                    row_index =  df_out.index[df_out['species_name'] == gas_product][0]
                    reactant_energy += float(df_out['formation_energy'].iloc[row_index]) * num_units
        
                if gas_product in dft_corrections_gases:
                    reactant_energy += dft_corrections_gases[gas_product] * num_units

    # Compute Adsorption Energy at U_RHE = 0 V
    adsorption_energy_RHE0 = reaction_energy + product_energy - reactant_energy

    # Apply solvation energy corrections
    if species_value in adsorbate_parameters['solvation_corrections_adsorbates']:
        solvation_correction = adsorbate_parameters['solvation_corrections_adsorbates'][species_value]
    else:
        solvation_correction = 0.0

    # Apply field effects
    (U_RHE_energy_contribution, U_SHE_energy_contribution) = get_electric_field_contribution(field_effects, species_value, reactants)
    return (adsorption_energy_RHE0, U_RHE_energy_contribution, U_SHE_energy_contribution, solvation_correction)

def make_mkm_input_files(db_filepath, adsorbate_parameters, field_effects,
                         df_out):
    system_dir_path = db_filepath.parent / f'{adsorbate_parameters["desired_surface"]}_{adsorbate_parameters["desired_facet"]}'
    Path.mkdir(system_dir_path, parents=True, exist_ok=True)
    energies_filepath = system_dir_path / f'energies_f{field_effects["epsilon"]:.2e}.txt'
    
    # with open(energies_filepath, 'w') as energies_file:
    #     df_out[mkm_write_columns].to_string(energies_file, index=False, justify='left')
    
    header = '\t'.join(['surface_name', 'site_name', 'species_name',
                        'formation_energy', 'frequencies', 'reference'])
    lines = [] # list of lines in the output
    for index, row in df_out.iterrows():
        line = '\t'.join([row['surface_name'], row['site_name'],
                          row['species_name'], f'{row["formation_energy"]:.4f}',
                          str(row['frequencies']), row['reference']])
        lines.append(line)

    lines = [header] + lines #add header to top
    input_file = '\n'.join(lines) #Join the lines with a line break

    with open(energies_filepath, 'w') as energies_file:
        energies_file.write(input_file)
    return None

def populate_chemical_symbols_dict(chemical_symbols_dict, last_chemical_symbol):
    if last_chemical_symbol in chemical_symbols_dict:
        chemical_symbols_dict[last_chemical_symbol] += 1
    else:
        chemical_symbols_dict[last_chemical_symbol] = 1
    return chemical_symbols_dict

def formula_to_chemical_symbols(formula):
    "Return dictionary mapping chemical symbols to number of atoms"

    chemical_symbols_dict = {}

    # split chemical formula string into alpha and numeric characters
    regex = re.compile('(\d+|\s+)')
    split_formula = regex.split(formula)
    if '' in split_formula:
        split_formula.remove('')

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
            str_len = len(string)
            while str_len > 0:
                if string[:2] in chemical_symbols:
                    last_chemical_symbol = string[:2]
                    string = string[2:]
                    str_len -= 2
                elif string[:1] in chemical_symbols:
                    last_chemical_symbol = string[:1]
                    string = string[1:]
                    str_len -= 1
                chemical_symbols_dict = populate_chemical_symbols_dict(chemical_symbols_dict, last_chemical_symbol)
                    

    # multiply number of atoms for each chemical symbol with number of formula units
    for key in chemical_symbols_dict.keys():
        chemical_symbols_dict[key] = formula_unit_count * chemical_symbols_dict[key]
    return chemical_symbols_dict

def U_RHE_to_field(U_RHE, pH, U_M_PZC, d, temp):
    U_SHE = U_RHE - kB * temp * np.log(10) * pH
    field = (U_SHE - U_M_PZC) / d
    return (U_SHE, field)

def U_SHE_to_field(U_SHE, pH, U_M_PZC, d, temp):
    U_RHE = U_SHE + kB * temp * np.log(10) * pH
    field = (U_SHE - U_M_PZC) / d
    return (U_RHE, field)

def field_to_voltage(field, pH, U_M_PZC, d, temp):
    U_SHE = field * d + U_M_PZC
    U_RHE = U_SHE + kB * temp * np.log(10) * pH
    return (U_SHE, U_RHE)
