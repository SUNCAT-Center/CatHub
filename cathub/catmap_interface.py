import re
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sqlalchemy import create_engine 
from ase.db import connect
from ase.neb import fit0
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

# NOTE: use cathub function to convert db to dataframe
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
                   write_adsorbates, write_transition_states,
                   gas_jsondata_filepath, ads_jsondata_filepath,
                   ts_jsondata_filepath, rxn_expressions_filepath, ts_data,
                   temp, write_mkm_input_files, verbose=True):

    df_out = pd.DataFrame(columns=write_columns)

    if verbose:
        system_header = f'###### {adsorbate_parameters["desired_surface"]}_{adsorbate_parameters["desired_facet"]}: electric field strength = {field_effects["epsilon"]:.2f} V/A ######'
        print('-' * len(system_header))
        print(system_header)
        print('-' * len(system_header))

        print()
        print('Term1 = Calculated Electronic Energy + DFT Correction')
        print('Term2 = Enthalpic Temperature Correction + Entropy Contribution')
        print('Term3 = RHE-scale Dependency')
        print('Term4 = Solvation Correction + Electric Field Correction')
        print('Chemical Potential, µ = Term1 + Term2 + Term3 + Term4')
        print(f'Free Energy Change, ∆G = µ_species - µ_ref. For example, ∆G_CH4 = µ_CH4 - (µ_CO + 3 * µ_H2(ref) - µ_H2O)')
        print(f'∆G at U_RHE=0.0 V = ∆G - Term3')
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
    
    if write_transition_states:
        df_out = write_ts_energies(db_filepath, df_out, ts_jsondata_filepath,
                                   rxn_expressions_filepath, ts_data,
                                   adsorbate_parameters, reference_gases,
                                   dft_corrections_gases, field_effects,
                                   temp, verbose)
        
    # write corrected energy data to mkm input file
    if write_mkm_input_files:
        make_mkm_input_files(db_filepath, adsorbate_parameters, field_effects,
                             df_out)
    return df_out

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
        gas_vibration_data = json.load(f)

    vibrational_energies = {}
    species_list = [species_data['species'] for species_data in gas_vibration_data]
    for species_data in gas_vibration_data:
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
            if species_name == 'H':
                print(field_effects['epsilon'])
                print(raw_energy)
                print()
            elec_energy_calc.append(row.energy
                                    + (x - z) * reference_gas_energies['H2O']
                                    - x * reference_gas_energies['CO']
                                    - (x - z + y / 2) * reference_gas_energies['H2'])
            dft_corr.append(dft_corrections_gases[species_name] if species_name in dft_corrections_gases else 0.0)
            helm_offset.append(beef_dft_helmholtz_offset[species_name] if species_name in beef_dft_helmholtz_offset else 0.0)
            
            if species_name in species_list:
                species_index = species_list.index(species_name)
                thermo = IdealGasThermo(vib_energies=vibrational_energies[species_name],
                                        geometry=gas_vibration_data[species_index]['geometry'],
                                        atoms=row.toatoms(),
                                        symmetrynumber=gas_vibration_data[species_index]['symmetry'],
                                        spin=gas_vibration_data[species_index]['spin'])
                # zero point energy correction
                zpe.append(thermo.get_ZPE_correction())
                # enthalpic temperature correction
                enthalpy.append(thermo.get_enthalpy(temp, verbose=False))
                S = thermo.get_entropy(temp, gas_vibration_data[species_index]['fugacity'],verbose=False)
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
                frequencies.append(gas_vibration_data[species_list.index(species_name)]['frequencies'])
                references.append(gas_vibration_data[species_list.index(species_name)]['reference'])
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
        table_headers = ["Species", "Term1 (eV)", "Term2 (eV)", "Term3 (eV)", "Term4 (eV)", "µ (eV)", "∆G (eV)", "∆G at U_RHE=0 (eV)"]
        for index, species_name in enumerate(df['species_name']):
            sub_table = []
            sub_table.extend([species_name,
                              f'{df["energy_vector"][index][0]:.{num_decimal_places}f}',
                              f'{df["energy_vector"][index][1]:.{num_decimal_places}f}',
                              f'{df["energy_vector"][index][2]:.{num_decimal_places}f}',
                              f'{df["energy_vector"][index][3]:.{num_decimal_places}f}',
                              f'{df["energy_vector"][index][4]:.{num_decimal_places}f}',
                              f'{df["energy_vector"][index][5]:.{num_decimal_places}f}',
                              f'{df["energy_vector"][index][5] - df["energy_vector"][index][2]:.{num_decimal_places}f}'])
            table.append(sub_table)
        print(tabulate(table, headers=table_headers, tablefmt='psql', colalign=("right", ) * len(table_headers), disable_numparse=True))
        print('\n')
    return df_out

def get_electric_field_contribution(field_effects, species_value,
                                    reactants=None, beta=None):
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

    if beta:
        U_RHE_energy_contribution = beta * U_RHE
    else:
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
    desired_surface = adsorbate_parameters['desired_surface']
    desired_facet = adsorbate_parameters['desired_facet']
    df1 = df1[df1['surface_composition'] == desired_surface]
    df1 = df1[df1['facet'].str.contains(desired_facet)]

    # Load vibrational data
    with open(ads_jsondata_filepath) as f:
        ads_vibration_data = json.load(f)
    vibrational_energies = {}
    json_species_list = [species_data['species'] for species_data in ads_vibration_data]
    for species_data in ads_vibration_data:
        ads_species = species_data['species']
        vibrational_energies[ads_species] = []
        for vibrational_frequency in species_data['frequencies']:
            vibrational_energies[ads_species].append(vibrational_frequency * cm2eV)

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
            
            # enthalpic temperature correction
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
            frequencies.append(ads_vibration_data[json_species_list.index(species_name)]['frequencies'])
            references.append(ads_vibration_data[json_species_list.index(species_name)]['reference'])
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
        table_headers = ["Species", "Term1 (eV)", "Term2 (eV)", "Term3 (eV)", "Term4 (eV)", "µ (eV)", "∆G (eV)", "∆G at U_RHE=0 (eV)"]
        for index, species_name in enumerate(df3['species_name']):
            sub_table = []
            sub_table.extend([species_name,
                              f'{df3["energy_vector"][index][0]:.{num_decimal_places}f}',
                              f'{df3["energy_vector"][index][1]:.{num_decimal_places}f}',
                              f'{df3["energy_vector"][index][2]:.{num_decimal_places}f}',
                              f'{df3["energy_vector"][index][3]:.{num_decimal_places}f}',
                              f'{df3["energy_vector"][index][4]:.{num_decimal_places}f}',
                              f'{df3["energy_vector"][index][5]:.{num_decimal_places}f}',
                              f'{df3["energy_vector"][index][5] - df3["energy_vector"][index][2]:.{num_decimal_places}f}'])
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

def get_solvation_layer_charge(src_path, adsorbate, bond_distance_cutoff):

    bohr = 0.52917721092  # angstrom
    bader_charges_filename = 'bader_charges.txt'
    coordinates_filename = 'ACF.dat'

    chemical_symbols_dict = formula_to_chemical_symbols(adsorbate)

    element_list = []
    bader_charge_list = []
    bader_charges_filepath = src_path / bader_charges_filename
    coordinates_filepath = src_path / coordinates_filename
    with open(bader_charges_filepath, 'r') as bader_charges_file:
        for line_index, line in enumerate(bader_charges_file):
            line_elements = line.split()
            element_list.append(line_elements[3])
            bader_charge_list.append(float(line_elements[5]))

    bader_charge_array = np.asarray(bader_charge_list)
    num_atoms = len(element_list)
    coordinates = np.loadtxt(coordinates_filepath, skiprows=2, max_rows=num_atoms)[:, 1:4] * bohr
    z_coordinates = coordinates[:, 2]
    total_indices = np.arange(len(z_coordinates)).tolist()

    chemical_symbol_to_index_list = {}
    for chemical_symbol in chemical_symbols_dict:
        chemical_symbol_to_index_list[chemical_symbol] = [i for i, x in enumerate(element_list) if x == chemical_symbol]

    chemical_symbols_to_sorted_indices = {}
    anchor_first_run = 1  # first run
    for chemical_symbol, num_atoms in chemical_symbols_dict.items():
        # identify indices of the chemical symbol with lowest z-coordinate
        chemical_symbol_indices = np.asarray([i for i, x in enumerate(element_list) if x == chemical_symbol])
        chemical_symbol_z_coordinates = z_coordinates[chemical_symbol_indices]
        sort_indices = chemical_symbol_z_coordinates.argsort()
        chemical_symbols_to_sorted_indices[chemical_symbol] = chemical_symbol_indices[sort_indices]
        if anchor_first_run:
            anchor_chemical_symbol = chemical_symbol
            anchor_atom_index = chemical_symbols_to_sorted_indices[chemical_symbol][0]
            anchor_first_run = 0
        elif z_coordinates[chemical_symbols_to_sorted_indices[chemical_symbol][0]] < z_coordinates[anchor_atom_index]:
            anchor_chemical_symbol = chemical_symbol
            anchor_atom_index = chemical_symbols_to_sorted_indices[chemical_symbol][0]

    anchor_z_coordinate = z_coordinates[anchor_atom_index]
    substrate_indices = np.where(z_coordinates < anchor_z_coordinate)[0].tolist()
    non_substrate_indices = [index for index in total_indices if index not in substrate_indices]

    adsorbate_scrape = {}
    num_atoms_to_scrape = sum(chemical_symbols_dict.values())
    for chemical_symbol, num_atoms in chemical_symbols_dict.items():
        adsorbate_scrape[chemical_symbol] = num_atoms

    adsorbate_indices = [anchor_atom_index]
    adsorbate_scrape[anchor_chemical_symbol] -= 1
    num_atoms_to_scrape -= 1

    reference_atom_indices = [anchor_atom_index]
    while num_atoms_to_scrape:
        for reference_atom_index in reference_atom_indices:
            distance_to_ref = np.linalg.norm(coordinates[non_substrate_indices] - coordinates[reference_atom_index], axis=1)
            bonding_subindices_to_ref = np.where((distance_to_ref > 0) & (distance_to_ref < bond_distance_cutoff))[0]
            distance_to_subindices = distance_to_ref[bonding_subindices_to_ref]
            sorted_bonding_subindices_to_ref = bonding_subindices_to_ref[np.argsort(distance_to_subindices)]
            bonding_atom_indices_to_ref = [non_substrate_indices[index] for index in sorted_bonding_subindices_to_ref if non_substrate_indices[index] not in adsorbate_indices]
            reference_atom_indices = bonding_atom_indices_to_ref[:]
            for atom_index in bonding_atom_indices_to_ref:
                chemical_symbol = element_list[atom_index]
                if chemical_symbol in adsorbate_scrape:
                    if adsorbate_scrape[chemical_symbol]:
                        adsorbate_indices.append(atom_index)
                        adsorbate_scrape[chemical_symbol] -= 1
                        num_atoms_to_scrape -= 1
    
    solvation_layer_indices = [index for index in non_substrate_indices if index not in adsorbate_indices]
    solvation_layer_charges = bader_charge_array[solvation_layer_indices]
    solvation_layer_charge = solvation_layer_charges.sum()
    return solvation_layer_charge

def read_qe_log(file_path, wf_dipole_index):
    # read quantum espresso log file and return converged final energy and workfunction
    wf_line_index = -1
    energy_line_index = -1
    with open(file_path) as file:
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

def compute_barrier_extrapolation(ts_species, beta, phi_correction, v_extra,
                                  energy_data, workfunction_data, charge_data):

    # energy_data = [E_TS, E_FS]
    # workfunction_data = [phi_TS, phi_FS]
    # charge_data = [q_TS, q_FS]

    phi_TS_corr = workfunction_data[0] - phi_correction
    phi_FS_corr = workfunction_data[1] - phi_correction

    del_E = energy_data[0] - energy_data[1]
    del_q = charge_data[0] - charge_data[1]
    del_phi = phi_TS_corr - phi_FS_corr
    
    # backward barrier
    E_r = del_E + 0.5 * del_q * del_phi
    
    ## convert v_extra from SHE to RHE at given pH_out
    E_r_extrapolated = E_r + del_q * (phi_FS_corr - v_extra)
    return E_r_extrapolated

def write_ts_energies(db_filepath, df_out, ts_jsondata_filepath,
                      rxn_expressions_filepath, ts_data, adsorbate_parameters,
                      reference_gases, dft_corrections_gases, field_effects,
                      temp, verbose):

    # identify system ids for transition state species
    table_name = 'reaction'
    df1 = db_to_dataframe(table_name, str(db_filepath))
    desired_surface = adsorbate_parameters['desired_surface']
    desired_facet = adsorbate_parameters['desired_facet']
    df1 = df1[df1['surface_composition'] == desired_surface]
    df1 = df1[df1['facet'].str.contains(desired_facet)]

    # Load vibrational data
    with open(ts_jsondata_filepath) as f:
        ts_vibration_data = json.load(f)
    vibrational_energies = {}
    json_species_list = [species_data['species'] for species_data in ts_vibration_data]
    for species_data in ts_vibration_data:
        ts_species = species_data['species']
        vibrational_energies[ts_species] = []
        for vibrational_frequency in species_data['frequencies']:
            vibrational_energies[ts_species].append(vibrational_frequency * cm2eV)
    
    # Load reaction expression data
    (rxn_expressions,
     reactants_rxn_expressions,
     products_rxn_expressions,
     ts_states_rxn_expressions,
     beta_list_rxn_expressions) = read_reaction_expression_data(
                                                    rxn_expressions_filepath)
    df_activation = df1[df1['activation_energy'].notna()]
    ts_states_user_input = ts_data['ts_states']
    reaction_index_map = []
    beta_list_map = []
    for ts_state in ts_states_user_input:
        if ts_state in ts_states_rxn_expressions:
            reaction_index = ts_states_rxn_expressions.index(ts_state)
            reaction_index_map.append(reaction_index)
            beta_list_map.append(beta_list_rxn_expressions[reaction_index])

    ## build dataframe data for transition state species
    db = connect(str(db_filepath))
    surface, site, species, raw_energy, facet = [], [], [], [], []
    forward_barrier, backward_barrier = [], []
    dft_corr, zpe, enthalpy, entropy, rhe_corr = [], [], [], [], []
    solv_corr, formation_energy, efield_corr, alk_corr, extrapolation_corr = [], [], [], [], []
    energy_vector, frequencies, references = [], [], []

    # simple reaction species: only one active product and filter out reactions without any transition state species
    df_activation_copy = df_activation.copy()

    df_activation_copy_reactant_dict_col = df_activation_copy.reactants.apply(json.loads)
    for index, reactant in enumerate(df_activation_copy_reactant_dict_col):
        row_index = df_activation_copy.index[index]
        new_dict = {}
        if 0 in reactant.values():
            for key, value in reactant.items():
                if value:
                    new_dict[key] = value
            df_activation_copy.at[row_index, 'reactants'] = json.dumps(new_dict)
    df_activation_copy_reactant_dict_col = df_activation_copy.reactants.apply(json.loads)

    df_activation_copy_product_dict_col = df_activation_copy.products.apply(json.loads)
    for index, product in enumerate(df_activation_copy_product_dict_col):
        row_index = df_activation_copy.index[index]
        new_dict = {}
        if 0 in product.values():
            for key, value in product.items():
                if value:
                    new_dict[key] = value
            df_activation_copy.at[row_index, 'products'] = json.dumps(new_dict)
    df_activation_copy_product_dict_col = df_activation_copy.products.apply(json.loads)
    
    df_index_list = []
    df_index_map = []
    for reaction_index in reaction_index_map:
        df_indices_product = df_activation_copy_product_dict_col[df_activation_copy_product_dict_col == products_rxn_expressions[reaction_index]].index.values
        df_indices_reactant = df_activation_copy_reactant_dict_col[df_activation_copy_reactant_dict_col == reactants_rxn_expressions[reaction_index]].index.values
        df_indices = np.intersect1d(df_indices_reactant, df_indices_product)
        if len(df_indices):
            df_index_list.append(df_indices[0])
            df_index_map.append(df_indices[0])
        else:
            df_index_map.append('')
    df_activation_rxns = df_activation.loc[df_index_list]

    products_list = []
    species_list = []
    beta_list = []
    snapshot_range_list = []
    for index, df_index in enumerate(df_index_map):
        if df_index:
            species_list.append(ts_states_user_input[index])
            products_list.append(json.loads(df_activation.products.loc[df_index]))
            beta_list.append(beta_list_map[index])
            snapshot_range_list.append(ts_data['rxn_pathway_image_ids'][index])
    
    phi_correction = ts_data['phi_correction']
    v_extra = ts_data['v_extra']
    energy_data = ts_data['energy_data']
    workfunction_data = ts_data['workfunction_data']
    charge_data = ts_data['charge_data']
    for species_index, species_name in enumerate(species_list):
        if '-' in desired_surface:
            surface.append(desired_surface.split('-')[0])
        else:
            surface.append(desired_surface)
        site.append(desired_facet)
        species.append(species_name)

        # [adsorption_energy_RHE0, U_RHE_energy_contribution, U_SHE_energy_contribution, solvation_correction]
        beta = beta_list[species_index]
        snapshot_range = snapshot_range_list[species_index]
        (site_wise_energy_contributions, facet_list) = get_ts_energies(
                    df_activation_rxns, df_out, db_filepath, species_list,
                    species_name, products_list, snapshot_range, reference_gases,
                    dft_corrections_gases, adsorbate_parameters, field_effects,
                    beta)
        site_wise_energy_contributions = np.asarray(site_wise_energy_contributions)

        site_wise_adsorption_energies = np.sum(site_wise_energy_contributions, axis=1)
        min_adsorption_energy = min(site_wise_adsorption_energies)
        min_index = np.where(site_wise_adsorption_energies == min_adsorption_energy)[0][0]
        facet.append(facet_list[min_index])
        raw_energy.append(float("nan"))
        # forward barrier
        forward_barrier.append(site_wise_energy_contributions[min_index][0])
        # backward barrier
        backward_barrier.append(site_wise_energy_contributions[min_index][1])
        # Zero DFT Correction for transition states
        dft_corr.append(0.0)

        if species_name in json_species_list:
            thermo = HarmonicThermo(vib_energies=vibrational_energies[species_name])
            
            # zero point energy correction
            zpe.append(np.sum(vibrational_energies[species_name]) / 2.0)
            
            # enthalpic temperature correction
            enthalpy.append(thermo.get_internal_energy(temp,verbose=False))
            
            S = thermo.get_entropy(temp, verbose=False)
            # entropy contribution
            entropy.append(- temp * S)
        else:
            zpe.append(0.0)
            enthalpy.append(0.0)
            entropy.append(0.0)

        rhe_corr.append(site_wise_energy_contributions[min_index][2])
        solv_corr.append(site_wise_energy_contributions[min_index][4])
        
        if field_effects:
            efield_corr.append(site_wise_energy_contributions[min_index][3])

        # Apply alkaline correction
        alk_corr.append(ts_data['alk_corr'] if beta else 0.0)
        
        # Apply charge extrapolation scheme
        fin_ads_energy = 0
        for product, num_products in products_list[species_index].items():
            if 'gas' in product:
                noncatmap_style_species = product.replace('gas', '')
                idx1 = df_out.index[df_out['site_name'] == 'gas']
                idx2 = df_out.index[df_out['species_name'] == noncatmap_style_species]
                idx = idx1.intersection(idx2)
                if len(idx) == 1:
                    fin_ads_energy += num_products * df_out.formation_energy[idx[0]]
            elif 'star' in product:
                noncatmap_style_species = product.replace('star', '')
                if noncatmap_style_species:
                    idx1 = df_out.index[df_out['site_name'] != 'gas']
                    idx2 = df_out.index[df_out['species_name'] == noncatmap_style_species]
                    idx = idx1.intersection(idx2)
                    if len(idx) == 1:
                        fin_ads_energy += num_products * df_out.formation_energy[idx[0]]
        extrapolation_corr.append(compute_barrier_extrapolation(
                                    species_name, beta, phi_correction, v_extra,
                                    energy_data[species_index],
                                    workfunction_data[species_index],
                                    charge_data[species_index]))

        # compute energy vector
        term1_forward = forward_barrier[-1] + dft_corr[-1]
        term1_backward = backward_barrier[-1] + dft_corr[-1]
        term2 = enthalpy[-1] + entropy[-1]
        term3 = rhe_corr[-1]
        term4 = solv_corr[-1] + efield_corr[-1] + alk_corr[-1] + extrapolation_corr[-1] + fin_ads_energy
        G = mu = term1_backward + term2 + term3 + term4
        energy_vector.append([term1_backward, term2, term3, term4, mu, G])
        formation_energy.append(term1_backward + term4)
        
        if species_name in json_species_list:
            frequencies.append(
                ts_vibration_data[json_species_list.index(species_name)]['frequencies'])
            references.append(
                ts_vibration_data[json_species_list.index(species_name)]['reference'])
        else:
            frequencies.append([])
            references.append('')

    df3 = pd.DataFrame(list(zip(surface, site, species, raw_energy,
                                backward_barrier, dft_corr, zpe, enthalpy,
                                entropy, rhe_corr, solv_corr, formation_energy,
                                energy_vector, frequencies, references)),
                       columns=write_columns)
    df_out = df_out.append(df3, ignore_index=True, sort=False)

    if verbose:
        ts_phase_header = 'Transition State Free Energy Correction:'
        print(ts_phase_header)
        print('-' * len(ts_phase_header))
        
        table = []
        table_headers = ["Species", "Term1 (eV)", "Term2 (eV)", "Term3 (eV)", "Term4 (eV)", "∆G (eV)", "∆G at U_RHE=0 (eV)"]
        for index, species_name in enumerate(df3['species_name']):
            sub_table = []
            sub_table.extend([species_name,
                              f'{df3["energy_vector"][index][0]:.{num_decimal_places}f}',
                              f'{df3["energy_vector"][index][1]:.{num_decimal_places}f}',
                              f'{df3["energy_vector"][index][2]:.{num_decimal_places}f}',
                              f'{df3["energy_vector"][index][3]:.{num_decimal_places}f}',
                              f'{df3["energy_vector"][index][5]:.{num_decimal_places}f}',
                              f'{df3["energy_vector"][index][5] - df3["energy_vector"][index][2]:.{num_decimal_places}f}'])
            table.append(sub_table)
        print(tabulate(table, headers=table_headers, tablefmt='psql', colalign=("right", ) * len(table_headers), disable_numparse=True))
        print('\n')
    return df_out

def get_catmap_style_species(species):
    if 'H_g' in species:
        split_species = re.split('(\d+)', species)
        if '' in split_species:
            split_species.remove('')
        catmap_species = 'H2gas'
        if len(split_species) > 1:
            num_species = float(split_species[0]) * 0.5
        else:
            num_species = 0.5
    else:
        if species[0].isnumeric():
            split_species = re.split('(\d+)', species)
            if '' in split_species:
                split_species.remove('')
            num_species_digits = len(split_species[0])
            catmap_species = species[num_species_digits:]
            num_species = float(split_species[0])
        else:
            split_species = [species]
            num_species = 1
            catmap_species = species
        if '*_t' in catmap_species:
            catmap_species = catmap_species.replace('*_t', 'star')
        elif '_t' in catmap_species:
            catmap_species = catmap_species.replace('_t', 'star')
        elif '_g' in catmap_species:
            catmap_species = catmap_species.replace('_g', 'gas')
    return (catmap_species, num_species)

def read_reaction_expression_data(rxn_expressions_filepath):
    exec(compile(open(rxn_expressions_filepath, 'rb').read(), '<string>', 'exec'))
    discard_species_list = ['*_t', '_t', '_g']
    reactants_list, products_list, ts_states, beta_list = [], [], [], []
    if 'rxn_expressions' in locals():
        for rxn_expression in locals()['rxn_expressions']:
            if ';' in rxn_expression:
                (rxn, beta) = rxn_expression.split(';')
            else:
                rxn = rxn_expression
                beta = 0.0  # chemical reaction if beta not specified

            if isinstance(beta, str):
                if '=' in beta:
                    beta_list.append(float(beta.split('=')[1]))
            else:
                beta_list.append(beta)

            rxn_no_spaces = re.sub(' ', '', rxn)
            split_rxn = re.split('->|<->', rxn_no_spaces)
            
            reactant_term = split_rxn[0]
            product_term = split_rxn[-1]
            
            if '+' in reactant_term:
                reactant_species = reactant_term.split('+')
                for discard_species in discard_species_list:
                    if discard_species in reactant_species:
                        reactant_species.remove(discard_species)
                reactants_list.append(reactant_species)
            else:
                reactants_list.append([reactant_term])

            if '+' in product_term:
                product_species = product_term.split('+')
                for discard_species in discard_species_list:
                    if discard_species in reactant_species:
                        product_species.remove(discard_species)
                products_list.append(product_species)
            else:
                products_list.append([product_term])

            if rxn_no_spaces.count('->') == 2:
                ts_term = split_rxn[1]
                if '+' in ts_term:
                    ts_species = ts_term.split('+')
                    for discard_species in discard_species_list:
                        if discard_species in ts_species:
                            ts_species.remove(discard_species)
                    ts_states.append(ts_species[0])
                else:
                    ts_states.append(ts_term)
            else:
                ts_states.append('')

    new_reactants_list = []
    for reactant_list in reactants_list:
        new_reactant_dict = {}
        for reactant in reactant_list:
            if 'ele_g' not in reactant:
                (catmap_species, num_species) = get_catmap_style_species(reactant)
                if catmap_species in new_reactant_dict:
                    new_reactant_dict[catmap_species] += num_species
                else:
                    new_reactant_dict[catmap_species] = num_species
        new_reactants_list.append(new_reactant_dict)

    new_products_list = []
    for product_list in products_list:
        new_product_dict = {}
        for product in product_list:
            if 'ele_g' not in product:
                (catmap_species, num_species) = get_catmap_style_species(product)
                if catmap_species in new_product_dict:
                    new_product_dict[catmap_species] += num_species
                else:
                    new_product_dict[catmap_species] = num_species
        new_products_list.append(new_product_dict)

    new_ts_states = []
    for ts_state in ts_states:
        for discard_species in discard_species_list:
            ts_state = ts_state.replace(discard_species, '') 
        new_ts_states.append(ts_state)
    return (locals()['rxn_expressions'], new_reactants_list,
            new_products_list, new_ts_states, beta_list)

def get_ts_energies(df, df_out, db_filepath, species_list, species_value,
                    products_list, snapshot_range, reference_gases,
                    dft_corrections_gases, adsorbate_parameters, field_effects,
                    beta):
    "Compute electronic transition state energies for a given species at all suitable adsorption sites at a given U_SHE/RHE"
    
    indices = [index for index, value in enumerate(species_list) if value == species_value]
    facet_list = df.facet.iloc[indices].tolist()

    site_wise_energy_contributions = []
    for index, reaction_index in enumerate(indices):
        facet = facet_list[index]
        reactants = json.loads(df.reactants.iloc[reaction_index])
        products = products_list[reaction_index]
        reaction_energy = df.reaction_energy.iloc[reaction_index]
        (forward_barrier,
         backward_barrier,
         U_RHE_energy_contribution,
         U_SHE_energy_contribution,
         solvation_correction) = get_ts_energy(
             df_out, db_filepath, species_value, reactants, snapshot_range,
             reference_gases, dft_corrections_gases, adsorbate_parameters,
             field_effects, beta)
        site_wise_energy_contributions.append([forward_barrier, backward_barrier, U_RHE_energy_contribution, U_SHE_energy_contribution, solvation_correction])
    return (site_wise_energy_contributions, facet_list)

def get_ts_energy(df_out, db_filepath, species_value, reactants, snapshot_range,
                  reference_gases, dft_corrections_gases, adsorbate_parameters,
                  field_effects, beta):
    "Compute energy barrier for an transition state species"

    db = connect(str(db_filepath))
    snapshot_positions = []
    snapshot_energies = []
    snapshot_forces = []
    for index, snapshot_id in enumerate(range(snapshot_range[0], snapshot_range[1]+1)):
        snapshot_positions.append(db.get(id=snapshot_id).toatoms().positions)
        snapshot_energies.append(db.get(id=snapshot_id).toatoms().get_potential_energy())
        snapshot_forces.append(db.get(id=snapshot_id).toatoms().get_forces())
        if index == 0:
            lattice_vectors = db.get(id=snapshot_id).toatoms().cell
            pbc = db.get(id=snapshot_id).toatoms().pbc

    _, snapshot_energies, _, snapshot_energies_fit, _ = fit0(snapshot_energies,
                                                             snapshot_forces,
                                                             snapshot_positions,
                                                             lattice_vectors,
                                                             pbc)
    initial_energy = snapshot_energies[0]
    final_energy = snapshot_energies[-1]    
    ts_energy = max(snapshot_energies_fit)

    # backward barrier
    backward_barrier = ts_energy - final_energy
    # forward barrier
    forward_barrier = ts_energy - initial_energy

    # Apply solvation energy corrections
    if species_value in adsorbate_parameters['solvation_corrections_adsorbates']:
        solvation_correction = adsorbate_parameters['solvation_corrections_adsorbates'][species_value]
    else:
        solvation_correction = 0.0

    # Apply field effects
    (U_RHE_energy_contribution, U_SHE_energy_contribution) = get_electric_field_contribution(field_effects, species_value, reactants, beta)
    return (forward_barrier, backward_barrier, U_RHE_energy_contribution, U_SHE_energy_contribution, solvation_correction)

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

def get_free_energy_change_species(df, species_name):
    species_index = df['species_name'][df['species_name'] == species_name].index[0]
    free_energy_change = df['energy_vector'][species_index][-1]
    return free_energy_change

def get_free_energy_change_reaction(df, reactants, products):
    free_energy_change_reactants = 0
    for reactant in reactants:
        free_energy_change_reactants += get_free_energy_change_species(df, reactant)

    free_energy_change_products = 0
    for product in products:
        free_energy_change_products += get_free_energy_change_species(df, product)

    free_energy_change = free_energy_change_products - free_energy_change_reactants
    return free_energy_change

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
