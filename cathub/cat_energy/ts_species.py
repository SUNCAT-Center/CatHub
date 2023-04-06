'''
Module with function definitions relating to transition state species
'''
import json

import numpy as np
import pandas as pd
from ase.db import connect
from ase.utils.forcecurve import fit_raw
from ase.thermochemistry import HarmonicThermo
from tabulate import tabulate

from cathub.cathubsql import CathubSQL
from .io import NUM_DECIMAL_PLACES, write_columns
from .conversion import read_reaction_expression_data, \
    formula_to_chemical_symbols, KB, CM2EV, get_electric_field_contribution


def compute_barrier_extrapolation(workfunction_data, phi_correction, phi_ref,
                                  beta, temp, pH):
    '''
    Compute charge extrapolated transition state barrier
    '''

    # workfunction_data = [phi_TS, phi_FS]

    phi_ts_corr = workfunction_data[0] - phi_correction
    phi_fs_corr = workfunction_data[1] - phi_correction

    del_phi = phi_ts_corr - phi_fs_corr

    # size extrapolation
    size_extrapolation = 0.5 * beta * del_phi

    # extrapolation to vacuum
    u_she = phi_fs_corr - phi_ref
    # converting from SHE to RHE scale
    u_rhe = u_she + np.log(10) * KB * temp * pH
    vacuum_extrapolation = beta * u_rhe

    energy_extrapolation = size_extrapolation + vacuum_extrapolation
    return energy_extrapolation

def write_ts_energies(db_filepath, df_out, ts_jsondata_filepath,
                      rxn_expressions, ts_data, system_parameters,
                      external_effects, verbose, latex):
    '''
    Function to compute and return energetics of transition state species
    '''

    # Data from local cathub .db file
    db = CathubSQL(filename=db_filepath)
    df1 = db.get_dataframe()

    desired_surface = system_parameters['desired_surface']
    desired_facet = system_parameters['desired_facet']
    df1 = df1[df1['surface_composition'] == desired_surface]
    df1 = df1[df1['facet'].str.contains(desired_facet)]

    temp = system_parameters['temp']
    pH = system_parameters['pH']

    # Load vibrational data
    with open(ts_jsondata_filepath, encoding='utf8') as f:
        ts_vibration_data = json.load(f)
    vibrational_energies = {}
    json_species_list = [species_data['species']
                         for species_data in ts_vibration_data]
    for species_data in ts_vibration_data:
        ts_species = species_data['species']
        vibrational_energies[ts_species] = []
        for vibrational_frequency in species_data['frequencies']:
            vibrational_energies[ts_species].append(vibrational_frequency
                                                    * CM2EV)

    # Load reaction expression data
    reactants_rxn_expressions, products_rxn_expressions = [], []
    ts_states_rxn_expressions, beta_list_rxn_expressions = [], []
    for rxn_expression in rxn_expressions:
        (reactant_dict, product_dict,
         ts_states, beta) = read_reaction_expression_data(rxn_expression)
        reactants_rxn_expressions.append(reactant_dict)
        products_rxn_expressions.append(product_dict)
        ts_states_rxn_expressions.append(ts_states)
        beta_list_rxn_expressions.append(beta)

    df_activation = df1[df1['activation_energy'].notna()]
    ts_states_user_input = ts_data['ts_states']
    reaction_index_map = []
    beta_list_map = []
    for ts_state in ts_states_user_input:
        if ts_state in ts_states_rxn_expressions:
            reaction_index = ts_states_rxn_expressions.index(ts_state)
            reaction_index_map.append(reaction_index)
            beta_list_map.append(beta_list_rxn_expressions[reaction_index])

    # build dataframe data for transition state species
    surface, site, species, raw_energy, facet = [], [], [], [], []
    forward_barrier, backward_barrier = [], []
    dft_corr, zpe, enthalpy, entropy, rhe_corr = [], [], [], [], []
    solv_corr, formation_energy, efield_corr, alk_corr = [], [], [], []
    if ts_data['extrapolation']:
        extrapolation_corr = []
    energy_vector, frequencies, references = [], [], []

    # simple reaction species: only one active product and filter out
    # reactions without any transition state species
    df_activation_copy = df_activation.copy()

    df_activation_copy_reactant_dict_col = df_activation_copy.reactants.apply(
        json.loads)
    for index, reactant in enumerate(df_activation_copy_reactant_dict_col):
        row_index = df_activation_copy.index[index]
        new_dict = {}
        if 0 in reactant.values():
            for key, value in reactant.items():
                if value:
                    new_dict[key] = value
            df_activation_copy.at[row_index, 'reactants'] = json.dumps(new_dict)
    df_activation_copy_reactant_dict_col = df_activation_copy.reactants.apply(
        json.loads)

    df_activation_copy_product_dict_col = df_activation_copy.products.apply(
        json.loads)
    for index, product in enumerate(df_activation_copy_product_dict_col):
        row_index = df_activation_copy.index[index]
        new_dict = {}
        if 0 in product.values() or "star" in product.keys():
            for key, value in product.items():
                if value and key != "star":
                    new_dict[key] = value
            df_activation_copy.at[row_index, 'products'] = json.dumps(new_dict)
    df_activation_copy_product_dict_col = df_activation_copy.products.apply(
        json.loads)

    df_index_list = []
    df_index_map = []
    for reaction_index in reaction_index_map:
        df_indices_product = df_activation_copy_product_dict_col[
            df_activation_copy_product_dict_col
            == products_rxn_expressions[reaction_index]].index.values
        df_indices_reactant = df_activation_copy_reactant_dict_col[
            df_activation_copy_reactant_dict_col
            == reactants_rxn_expressions[reaction_index]].index.values
        df_indices = np.intersect1d(df_indices_reactant, df_indices_product)
        if df_indices:
            df_index_list.append(df_indices[0])
            df_index_map.append(df_indices[0])
        else:
            df_index_map.append('')
    df_activation_rxns = df_activation.loc[df_index_list]

    products_list, species_list, beta_list, snapshot_range_list = [], [], [], []
    for index, df_index in enumerate(df_index_map):
        if df_index:
            species_list.append(ts_states_user_input[index])
            products_list.append(json.loads(
                df_activation.products.loc[df_index]))
            beta_list.append(beta_list_map[index])
            snapshot_range_list.append(ts_data['rxn_pathway_image_ids'][index])

    phi_correction = ts_data['phi_correction']
    phi_ref = ts_data['phi_ref']
    workfunction_data = ts_data['workfunction_data']
    for species_index, species_name in enumerate(species_list):
        if '-' in desired_surface:
            surface.append(desired_surface.split('-')[0])
        else:
            surface.append(desired_surface)
        site.append(desired_facet)
        species.append(species_name)

        # [adsorption_energy_rhe0, rhe_energy_contribution,
        # she_energy_contribution, solvation_correction]
        beta = beta_list[species_index]
        snapshot_range = snapshot_range_list[species_index]
        (site_wise_energy_contributions, facet_list) = get_ts_energies(
                        df_activation_rxns, db_filepath, species_list,
                        species_name, snapshot_range,
                        system_parameters, external_effects, beta)

        site_wise_energy_contributions = np.asarray(
                                                site_wise_energy_contributions)

        site_wise_adsorption_energies = np.sum(site_wise_energy_contributions,
                                               axis=1)
        min_adsorption_energy = min(site_wise_adsorption_energies)
        min_index = np.where(
                site_wise_adsorption_energies == min_adsorption_energy)[0][0]
        facet.append(facet_list[min_index])
        raw_energy.append(float("nan"))
        # forward barrier
        forward_barrier.append(site_wise_energy_contributions[min_index][0])
        # backward barrier
        backward_barrier.append(site_wise_energy_contributions[min_index][1])
        # Zero DFT Correction for transition states
        dft_corr.append(0.0)

        if species_name in json_species_list:
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

        rhe_corr.append(site_wise_energy_contributions[min_index][2])
        solv_corr.append(site_wise_energy_contributions[min_index][4])
        efield_corr.append(site_wise_energy_contributions[min_index][3]
                           if external_effects else 0.0)

        # Apply alkaline correction
        alk_corr.append(ts_data['alk_corr'] if beta else 0.0)

        # Compute final state energy
        fin_ads_energy = 0
        for product, num_products in products_list[species_index].items():
            if 'gas' in product:
                noncatmap_style_species = product.replace('gas', '')
                idx1 = df_out.index[df_out['site_name'] == 'gas']
                idx2 = (df_out.index[df_out['species_name']
                        == noncatmap_style_species])
                idx = idx1.intersection(idx2)
                if len(idx) == 1:
                    fin_ads_energy += (
                        num_products * df_out.formation_energy[idx[0]])
            elif 'star' in product:
                noncatmap_style_species = product.replace('star', '')
                if noncatmap_style_species:
                    idx1 = df_out.index[df_out['site_name'] != 'gas']
                    idx2 = (df_out.index[df_out['species_name']
                            == noncatmap_style_species])
                    idx = idx1.intersection(idx2)
                    if len(idx) == 1:
                        fin_ads_energy += (
                            num_products * df_out.formation_energy[idx[0]])

        # Apply charge extrapolation scheme
        if ts_data['extrapolation']:
            extrapolation_corr.append(compute_barrier_extrapolation(
                                            workfunction_data[species_index],
                                            phi_correction, phi_ref,
                                            beta_list[species_index], temp, pH))

        # compute energy vector
        # term1_forward = forward_barrier[-1] + dft_corr[-1]
        term1_backward = backward_barrier[-1] + dft_corr[-1]
        term2 = enthalpy[-1] + entropy[-1]
        term3 = rhe_corr[-1]
        if ts_data['extrapolation']:
            term4 = (solv_corr[-1] + efield_corr[-1] + alk_corr[-1]
                     + extrapolation_corr[-1] + fin_ads_energy)
        else:
            term4 = (solv_corr[-1] + efield_corr[-1]
                     + alk_corr[-1] + fin_ads_energy)
        G = mu = term1_backward + term2 + term3 + term4
        energy_vector.append([term1_backward, term2, term3, term4, mu, G])
        formation_energy.append(term1_backward + term4)

        if species_name in json_species_list:
            frequencies.append(ts_vibration_data[json_species_list.index(
                species_name)]['frequencies'])
            references.append(ts_vibration_data[json_species_list.index(
                species_name)]['reference'])
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

        print('Term1 = Backward Electronic Activation Energy + DFT Correction')
        print('Term2 = Enthalpic Temperature Correction + Entropy Contribution')
        print('Term3 = RHE-scale Dependency')  
        if ts_data['extrapolation']:
            print('Term4 = External Effect Corrections + Alkaline Correction'
                  ' + Charge Extrapolation Correction + Final Adsorbate Energy')
        else:
            print('Term4 = External Effect Corrections + Alkaline Correction'
                  ' + Final Adsorbate Energy')
        print('Free Energy Change, ∆G = Term1 + Term2 + Term3 + Term4')
        print('∆G at U_RHE=0.0 V = ∆G - Term3')
        print()

        table = []
        table_headers = ["Species", "Term1 (eV)", "Term2 (eV)", "Term3 (eV)",
                         "Term4 (eV)", "∆G (eV)", "∆G at U_RHE=0 (eV)"]
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
                 f'{df3["energy_vector"][index][5]:.{NUM_DECIMAL_PLACES}f}',
                 f'{delg_at_zero_u_rhe:.{NUM_DECIMAL_PLACES}f}'])
            table.append(sub_table)
        if latex:
            print(tabulate(table, headers=table_headers, tablefmt='latex',
                           colalign=("right", ) * len(table_headers),
                           disable_numparse=True))
        else:
            print(tabulate(table, headers=table_headers, tablefmt='psql',
                           colalign=("right", ) * len(table_headers),
                           disable_numparse=True))
        print('\n')
    return df_out

def get_ts_energies(
        df, db_filepath, species_list, species_value,
        snapshot_range, adsorbate_parameters, field_effects, beta):
    '''
    Compute electronic transition state energies for a given species at all
    suitable adsorption sites at a given u_she/RHE
    '''

    indices = [index for index, value in enumerate(species_list)
               if value == species_value]
    facet_list = df.facet.iloc[indices].tolist()

    site_wise_energy_contributions = []
    for reaction_index in indices:
        # facet = facet_list[index]
        reactants = json.loads(df.reactants.iloc[reaction_index])
        # products = products_list[reaction_index]
        # reaction_energy = df.reaction_energy.iloc[reaction_index]
        (forward_barrier,
         backward_barrier,
         rhe_energy_contribution,
         she_energy_contribution,
         solvation_correction) = get_ts_energy(
             db_filepath, species_value, reactants, snapshot_range,
             adsorbate_parameters, field_effects, beta)
        site_wise_energy_contributions.append(
            [forward_barrier, backward_barrier, rhe_energy_contribution,
            she_energy_contribution, solvation_correction])
    return (site_wise_energy_contributions, facet_list)

def get_ts_energy(db_filepath, species_value, reactants, snapshot_range,
                  adsorbate_parameters, field_effects, beta):
    '''
    Compute energy barrier for an transition state species
    '''

    db = connect(str(db_filepath))
    snapshot_positions = []
    snapshot_energies = []
    snapshot_forces = []
    for index, snapshot_id in enumerate(range(snapshot_range[0],
                                              snapshot_range[1]+1)):
        snapshot_positions.append(db.get(id=snapshot_id).toatoms().positions)
        snapshot_energies.append(db.get(
                            id=snapshot_id).toatoms().get_potential_energy())
        snapshot_forces.append(db.get(id=snapshot_id).toatoms().get_forces())
        if index == 0:
            lattice_vectors = db.get(id=snapshot_id).toatoms().cell
            pbc = db.get(id=snapshot_id).toatoms().pbc

    _, snapshot_energies, _, snapshot_energies_fit, _ = fit_raw(snapshot_energies,
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
    if species_value in adsorbate_parameters['solvation_corrections']:
        solvation_correction = adsorbate_parameters['solvation_corrections'][
                                                                species_value]
    else:
        solvation_correction = 0.0

    # Apply field effects
    (rhe_energy_contribution, she_energy_contribution) = (
                get_electric_field_contribution(field_effects, species_value,
                                                reactants, beta))
    return (forward_barrier, backward_barrier, rhe_energy_contribution,
            she_energy_contribution, solvation_correction)

def get_solvation_layer_charge(src_path, adsorbate, bond_distance_cutoff):
    '''
    Compute the total charge of the solvation layer
    '''

    bohr = 0.52917721092  # angstrom
    bader_charges_filename = 'bader_charges.txt'
    coordinates_filename = 'ACF.dat'

    chemical_symbols_dict = formula_to_chemical_symbols(adsorbate)

    element_list = []
    bader_charge_list = []
    bader_charges_filepath = src_path / bader_charges_filename
    coordinates_filepath = src_path / coordinates_filename
    with open(bader_charges_filepath, 'r', encoding='utf8') as bader_charges_file:
        for line in bader_charges_file:
            line_elements = line.split()
            element_list.append(line_elements[3])
            bader_charge_list.append(float(line_elements[5]))

    bader_charge_array = np.asarray(bader_charge_list)
    num_atoms = len(element_list)
    coordinates = np.loadtxt(coordinates_filepath, skiprows=2,
                             max_rows=num_atoms)[:, 1:4] * bohr
    z_coordinates = coordinates[:, 2]
    total_indices = np.arange(len(z_coordinates)).tolist()

    chemical_symbol_to_index_list = {}
    for chemical_symbol in chemical_symbols_dict:
        chemical_symbol_to_index_list[chemical_symbol] = [
            i for i, x in enumerate(element_list) if x == chemical_symbol]

    chemical_symbols_to_sorted_indices = {}
    anchor_first_run = 1  # first run
    for chemical_symbol, num_atoms in chemical_symbols_dict.items():
        # identify indices of the chemical symbol with lowest z-coordinate
        chemical_symbol_indices = np.asarray(
            [i for i, x in enumerate(element_list) if x == chemical_symbol])
        chemical_symbol_z_coordinates = z_coordinates[chemical_symbol_indices]
        sort_indices = chemical_symbol_z_coordinates.argsort()
        chemical_symbols_to_sorted_indices[chemical_symbol] = (
            chemical_symbol_indices[sort_indices])
        if anchor_first_run:
            anchor_chemical_symbol = chemical_symbol
            anchor_atom_index = (
                chemical_symbols_to_sorted_indices[chemical_symbol][0])
            anchor_first_run = 0
        elif z_coordinates[chemical_symbols_to_sorted_indices[
                chemical_symbol][0]] < z_coordinates[anchor_atom_index]:
            anchor_chemical_symbol = chemical_symbol
            anchor_atom_index = chemical_symbols_to_sorted_indices[
                chemical_symbol][0]

    anchor_z_coordinate = z_coordinates[anchor_atom_index]
    substrate_indices = np.where(z_coordinates
                                 < anchor_z_coordinate)[0].tolist()
    non_substrate_indices = [
        index for index in total_indices if index not in substrate_indices]

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
            distance_to_ref = np.linalg.norm(
                coordinates[non_substrate_indices]
                - coordinates[reference_atom_index], axis=1)
            bonding_subindices_to_ref = np.where(
                (distance_to_ref > 0) & (distance_to_ref < bond_distance_cutoff))[0]
            distance_to_subindices = distance_to_ref[bonding_subindices_to_ref]
            sorted_bonding_subindices_to_ref = bonding_subindices_to_ref[
                np.argsort(distance_to_subindices)]
            bonding_atom_indices_to_ref = [
                non_substrate_indices[index]
                for index in sorted_bonding_subindices_to_ref
                if non_substrate_indices[index] not in adsorbate_indices]
            reference_atom_indices = bonding_atom_indices_to_ref[:]
            for atom_index in bonding_atom_indices_to_ref:
                chemical_symbol = element_list[atom_index]
                if chemical_symbol in adsorbate_scrape:
                    if adsorbate_scrape[chemical_symbol]:
                        adsorbate_indices.append(atom_index)
                        adsorbate_scrape[chemical_symbol] -= 1
                        num_atoms_to_scrape -= 1

    solvation_layer_indices = [index for index in non_substrate_indices
                               if index not in adsorbate_indices]
    solvation_layer_charges = bader_charge_array[solvation_layer_indices]
    solvation_layer_charge = solvation_layer_charges.sum()
    return solvation_layer_charge
