'''
Module with function definitions relating to adsorption state species
'''
import json

import numpy as np
import pandas as pd
from ase.thermochemistry import HarmonicThermo
from tabulate import tabulate

from .io import NUM_DECIMAL_PLACES, db_to_dataframe, write_columns
from .conversion import formula_to_chemical_symbols, CM2EV
from .conversion import get_electric_field_contribution


def write_adsorbate_energies(
        db_filepath, df_out, ads_jsondata_filepath, adsorbate_parameters,
        facet_conditional, reference_gases, dft_corrections_gases,
        field_effects, temp, verbose, latex):
    '''
    Write formation energies to energies.txt after applying free energy
    corrections
    '''

    # identify system ids for adsorbate species
    table_name = 'reaction'
    df1 = db_to_dataframe(table_name, str(db_filepath))
    desired_surface = adsorbate_parameters['desired_surface']
    desired_facet = adsorbate_parameters['desired_facet']
    df1 = df1[df1['surface_composition'] == desired_surface]
    df1 = df1[df1['facet'].str.contains(desired_facet)]

    # Load vibrational data
    with open(ads_jsondata_filepath, encoding='utf8') as f:
        ads_vibration_data = json.load(f)
    vibrational_energies = {}
    json_species_list = [species_data['species']
                         for species_data in ads_vibration_data]
    for species_data in ads_vibration_data:
        ads_species = species_data['species']
        vibrational_energies[ads_species] = []
        for vibrational_frequency in species_data['frequencies']:
            vibrational_energies[ads_species].append(vibrational_frequency
                                                     * CM2EV)

    # build dataframe data for adsorbate species
    surface, site, species, raw_energy, facet = [], [], [], [], []
    elec_energy_calc, dft_corr, zpe, enthalpy, entropy = [], [], [], [], []
    rhe_corr, solv_corr, formation_energy, efield_corr = [], [], [], []
    energy_vector, frequencies, references = [], [], []

    # simple reaction species: only one active product and filter out reactions
    # without any adsorbed species
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

        # [adsorption_energy_rhe0, rhe_energy_contribution,
        # she_energy_contribution, solvation_correction]
        (site_wise_energy_contributions, facet_list) = get_adsorption_energies(
            df2, df_out, species_list, species_name, products_list,
            reference_gases, dft_corrections_gases, adsorbate_parameters,
            facet_conditional, field_effects)
        site_wise_energy_contributions = np.asarray(
                                                site_wise_energy_contributions)

        site_wise_adsorption_energies = np.sum(site_wise_energy_contributions,
                                               axis=1)
        min_adsorption_energy = min(site_wise_adsorption_energies)
        min_index = np.where(
                site_wise_adsorption_energies == min_adsorption_energy)[0][0]
        facet.append(facet_list[min_index])
        raw_energy.append(float("nan"))
        elec_energy_calc.append(site_wise_energy_contributions[min_index][0])
        # Zero DFT Correction for Adsorbates
        dft_corr.append(0.0)

        if species_name in json_species_list:
            species_index = json_species_list.index(species_name)
            thermo = HarmonicThermo(
                                vib_energies=vibrational_energies[species_name])

            # zero point energy correction
            zpe.append(np.sum(vibrational_energies[species_name]) / 2.0)

            # enthalpic temperature correction
            enthalpy.append(thermo.get_internal_energy(temp, verbose=False))

            S = thermo.get_entropy(temp, verbose=False)
            # entropy contribution
            entropy.append(- temp * S)
        else:
            zpe.append(0.0)
            enthalpy.append(0.0)
            entropy.append(0.0)

        rhe_corr.append(site_wise_energy_contributions[min_index][1])
        solv_corr.append(site_wise_energy_contributions[min_index][3])
        efield_corr.append(site_wise_energy_contributions[min_index][2]
                           if field_effects else 0.0)

        # compute energy vector
        term1 = elec_energy_calc[-1] + dft_corr[-1]
        term2 = enthalpy[-1] + entropy[-1]
        term3 = rhe_corr[-1]
        term4 = solv_corr[-1] + efield_corr[-1]
        mu = term1 + term2 + term3 + term4
        energy_vector.append([term1, term2, term3, term4, mu])
        formation_energy.append(term1 + term4)

        if species_name in json_species_list:
            frequencies.append(ads_vibration_data[json_species_list.index(
                                                species_name)]['frequencies'])
            references.append(ads_vibration_data[json_species_list.index(
                                                    species_name)]['reference'])
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

        # CO2 Reduction Reaction
        if set(reference_gases) == set(['CO2', 'H2_ref', 'H2O']):
            G = (energy_vector[species_index][-1]
                 + (2 * x - z) * reference_mu['H2O']
                 - x * reference_mu['CO2']
                 - (2 * x - z + y / 2) * reference_mu['H2_ref'])
        # CO Reduction Reaction
        elif set(reference_gases) == set(['CO', 'H2_ref', 'H2O']):
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
        table_headers = ["Species", "Term1 (eV)", "Term2 (eV)", "Term3 (eV)",
                         "Term4 (eV)", "µ (eV)", "∆G (eV)",
                         "∆G at U_RHE=0 (eV)"]
        for index, species_name in enumerate(df3['species_name']):
            sub_table = []
            delg_at_zero_u_rhe = (df3["energy_vector"][index][5]
                                  - df3["energy_vector"][index][2])
            sub_table.extend(
                [species_name,
                f'{df3["energy_vector"][index][0]:.{NUM_DECIMAL_PLACES}f}',
                f'{df3["energy_vector"][index][1]:.{NUM_DECIMAL_PLACES}f}',
                f'{df3["energy_vector"][index][2]:.{NUM_DECIMAL_PLACES}f}',
                f'{df3["energy_vector"][index][3]:.{NUM_DECIMAL_PLACES}f}',
                f'{df3["energy_vector"][index][4]:.{NUM_DECIMAL_PLACES}f}',
                f'{df3["energy_vector"][index][5]:.{NUM_DECIMAL_PLACES}f}',
                f'{delg_at_zero_u_rhe:.{NUM_DECIMAL_PLACES}f}'])
            table.append(sub_table)
        if latex:
            print(tabulate(
                table, headers=table_headers, tablefmt='latex',
                colalign=("right", ) * len(table_headers),
                disable_numparse=True))
        else:
            print(tabulate(
                table, headers=table_headers, tablefmt='psql',
                colalign=("right", ) * len(table_headers),
                disable_numparse=True))
        print('\n')
    return df_out


def get_adsorption_energies(
        df, df_out, species_list, species_value, products_list, reference_gases,
        dft_corrections_gases, adsorbate_parameters, facet_conditional,
        field_effects):
    '''
    Compute electronic adsorption energies for a given species at all suitable
    adsorption sites at a given u_she/RHE
    '''

    indices = [index for index, value in enumerate(species_list)
                                                    if value == species_value]
    facet_list = df.facet.iloc[indices].tolist()

    site_wise_energy_contributions = []
    for index, reaction_index in enumerate(indices):
        facet = facet_list[index]
        if facet_conditional in facet:
            reactants = json.loads(df.reactants.iloc[reaction_index])
            products = products_list[reaction_index]
            reaction_energy = df.reaction_energy.iloc[reaction_index]
            (adsorption_energy_rhe0,
             rhe_energy_contribution,
             she_energy_contribution,
             solvation_correction) = get_adsorption_energy(
                 df_out, species_value, reactants, products, reaction_energy,
                 reference_gases, dft_corrections_gases, adsorbate_parameters,
                 field_effects)
            site_wise_energy_contributions.append(
                            [adsorption_energy_rhe0, rhe_energy_contribution,
                            she_energy_contribution, solvation_correction])
    return (site_wise_energy_contributions, facet_list)


def get_adsorption_energy(
        df_out, species_value, reactants, products, reaction_energy,
        reference_gases, dft_corrections_gases, adsorbate_parameters,
        field_effects):
    '''
    Compute adsorption energy for an adsorbate species in a given reaction
    '''

    product_energy = 0
    for product, num_units in products.items():
        if 'star' not in product:
            if 'gas' in product:
                gas_product = product.replace('gas', '')
                if gas_product not in reference_gases:
                    row_index = df_out.index[
                                    df_out['species_name'] == gas_product][0]
                    product_energy += float(
                        df_out['formation_energy'].iloc[row_index]) * num_units

                if gas_product in dft_corrections_gases:
                    product_energy += (dft_corrections_gases[gas_product]
                                                                    * num_units)

    reactant_energy = 0
    for reactant, num_units in reactants.items():
        if 'star' not in reactant:
            if 'gas' in reactant:
                gas_product = reactant.replace('gas', '')
                if (gas_product not in reference_gases
                            and (gas_product + '_ref') not in reference_gases):
                    row_index = (df_out.index[df_out['species_name']
                                                            == gas_product][0])
                    reactant_energy += (
                            float(df_out['formation_energy'].iloc[row_index])
                            * num_units)

                if gas_product in dft_corrections_gases:
                    reactant_energy += (dft_corrections_gases[gas_product]
                                        * num_units)

    # Compute Adsorption Energy at u_rhe = 0 V
    adsorption_energy_rhe0 = reaction_energy + product_energy - reactant_energy

    # Apply solvation energy corrections
    if species_value in adsorbate_parameters['solvation_corrections']:
        solvation_correction = adsorbate_parameters['solvation_corrections'][
                                                                species_value]
    else:
        solvation_correction = 0.0

    # Apply field effects
    (rhe_energy_contribution,
    she_energy_contribution) = get_electric_field_contribution(
                field_effects, species_value, reference_gases, reactants)
    return (adsorption_energy_rhe0, rhe_energy_contribution,
            she_energy_contribution, solvation_correction)