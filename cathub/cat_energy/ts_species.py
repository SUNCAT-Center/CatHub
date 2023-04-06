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
    formula_to_chemical_symbols, CM2EV, PHI_REF, get_rhe_contribution


def get_constant_charge_barriers(db_filepath, neb_image_id_range):
    '''
    Compute energy barrier for an transition state species
    '''

    db = connect(str(db_filepath))
    neb_image_id_positions = []
    neb_energies = []
    neb_forces = []
    for index, neb_image_id in enumerate(range(neb_image_id_range[0],
                                               neb_image_id_range[1] + 1)):
        neb_image_id_positions.append(db.get(id=neb_image_id).toatoms().positions)
        neb_energies.append(db.get(
            id=neb_image_id).toatoms().get_potential_energy())
        neb_forces.append(db.get(id=neb_image_id).toatoms().get_forces())
        if index == 0:
            lattice_vectors = db.get(id=neb_image_id).toatoms().cell
            pbc = db.get(id=neb_image_id).toatoms().pbc

    _, neb_energies, _, neb_energies_fit, _ = \
        fit_raw(neb_energies, neb_forces, neb_image_id_positions,
                lattice_vectors, pbc)
    initial_energy = neb_energies[0]
    final_energy = neb_energies[-1]
    ts_energy = max(neb_energies_fit)

    constant_charge_forward_barrier = ts_energy - initial_energy
    constant_charge_backward_barrier = ts_energy - final_energy
    return (constant_charge_forward_barrier, constant_charge_backward_barrier)


def get_constant_potential_barrier(constant_charge_barrier,
                                   delq, del_phi, barrier_workfunction_state='TS'):
    if barrier_workfunction_state == 'TS':
        constant_potential_barrier = constant_charge_barrier - delq * del_phi / 2
    else:
        constant_potential_barrier = constant_charge_barrier + delq * del_phi / 2
    return constant_potential_barrier


def get_extrapolated_barrier(constant_potential_barrier, barrier_workfunction,
                             extrapolated_workfunction, delq):
    extrapolated_barrier = (constant_potential_barrier
                            - delq * (extrapolated_workfunction - barrier_workfunction))
    return extrapolated_barrier


def get_charge_extrapolated_constant_potential_barriers(
        db_filepath, neb_image_id_range, species_workfunction_data, beta, u_she,
        extrapolate, alk_corr):
    '''Compute charge extrapolated constant potential energy barrier for a
    transition state species '''

    (constant_charge_forward_barrier, constant_charge_backward_barrier) = \
        get_constant_charge_barriers(db_filepath, neb_image_id_range)

    barrier_natures = ['forward', 'backward']
    barrier_workfunction_state = 'TS'  # TS is the default
    for barrier_nature in barrier_natures:
        if barrier_nature == 'forward':
            barrier = constant_charge_forward_barrier
            del_phi = species_workfunction_data[1] - species_workfunction_data[0]
            delq = - beta
            constant_potential_forward_barrier = get_constant_potential_barrier(
                barrier, delq, del_phi, barrier_workfunction_state)
        else:
            barrier = constant_charge_backward_barrier
            del_phi = species_workfunction_data[1] - species_workfunction_data[2]
            delq = 1 - beta
            constant_potential_backward_barrier = get_constant_potential_barrier(
                barrier, delq, del_phi, barrier_workfunction_state)

    # Apply charge extrapolation scheme
    if extrapolate:
        if barrier_workfunction_state == 'IS':
            barrier_workfunction = species_workfunction_data[0]
        elif barrier_workfunction_state == 'TS':
            barrier_workfunction = species_workfunction_data[1]
        elif barrier_workfunction_state == 'FS':
            barrier_workfunction = species_workfunction_data[2]

        workfunction_at_u_she = u_she + PHI_REF
        for barrier_nature in barrier_natures:
            if barrier_nature == 'forward':
                delq = - beta
                charge_extrapolated_constant_potential_forward_barrier = \
                    get_extrapolated_barrier(constant_potential_forward_barrier,
                                             barrier_workfunction,
                                             workfunction_at_u_she, delq)
            else:
                delq = 1 - beta
                charge_extrapolated_constant_potential_backward_barrier = \
                    get_extrapolated_barrier(constant_potential_backward_barrier,
                                             barrier_workfunction,
                                             workfunction_at_u_she, delq)
    else:
        charge_extrapolated_constant_potential_forward_barrier = \
            constant_potential_forward_barrier
        charge_extrapolated_constant_potential_backward_barrier = \
            constant_potential_backward_barrier

    # Apply alkaline correction
    charge_extrapolated_constant_potential_forward_barrier += alk_corr
    charge_extrapolated_constant_potential_backward_barrier += alk_corr

    return (charge_extrapolated_constant_potential_forward_barrier,
            charge_extrapolated_constant_potential_backward_barrier)


def write_ts_energies(
        db_filepath, df_out, ts_jsondata_filepath, rxn_expressions, ts_data,
        system_parameters, reference_gases, external_effects, verbose, latex):
    '''
    Write formation energies of transition state species to energies.txt after
    applying free energy corrections
    '''

    # Data from local cathub .db file
    db = CathubSQL(filename=db_filepath)
    df1 = db.get_dataframe()

    desired_surface = system_parameters['desired_surface']
    desired_facet = system_parameters['desired_facet']
    df1 = df1[df1['surface_composition'] == desired_surface]
    df1 = df1[df1['facet'].str.contains(desired_facet)]

    temp = system_parameters['temp']
    u_rhe = system_parameters['u_rhe']
    u_she = system_parameters['u_she']

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
    num_rxns = len(rxn_expressions)
    reactants_rxn_expressions, products_rxn_expressions = [], []
    ts_states_rxn_expressions, beta_list_rxn_expressions = [], []
    for rxn_expression in rxn_expressions:
        (reactant_dict, product_dict,
         ts_states, beta) = read_reaction_expression_data(rxn_expression)
        reactants_rxn_expressions.append(reactant_dict)
        products_rxn_expressions.append(product_dict)
        ts_states_rxn_expressions.append(ts_states)
        beta_list_rxn_expressions.append(beta)
    valid_ts_states_reaction_indices = [reaction_index for reaction_index, ts_state in enumerate(ts_states_rxn_expressions) if ts_state]

    df_activation = df1[df1['activation_energy'].notna()]

    # build dataframe data for transition state species
    surface, site, species, raw_energy = [], [], [], []
    charge_extrapolated_constant_potential_barriers = []
    dft_corr, zpe, enthalpy, entropy, rhe_corr = [], [], [], [], []
    term1_forward_list, term4_forward_list = [], []
    formation_energy, alk_corr, energy_vector, frequencies = [], [], [], []
    references = []

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

    df_index_rxn_expressions = []
    for reaction_index in range(num_rxns):
        df_indices_product = df_activation_copy_product_dict_col[
            df_activation_copy_product_dict_col
            == products_rxn_expressions[reaction_index]].index.values
        df_indices_reactant = df_activation_copy_reactant_dict_col[
            df_activation_copy_reactant_dict_col
            == reactants_rxn_expressions[reaction_index]].index.values
        df_indices = np.intersect1d(df_indices_reactant, df_indices_product)
        if df_indices:
            df_index_rxn_expressions.append(df_indices[0])
        else:
            df_index_rxn_expressions.append('')
    available_df_indices = [df_index for df_index in df_index_rxn_expressions if df_index != '']
    df_activation_rxns = df_activation.loc[available_df_indices]

    species_list, beta_list, reactants_list, products_list = [], [], [], []
    phi_correction = ts_data['phi_correction']
    for reaction_index in valid_ts_states_reaction_indices:
        ts_state = ts_states_rxn_expressions[reaction_index]
        species_list.append(ts_state)
        beta_list.append(beta_list_rxn_expressions[reaction_index])
        df_index = df_index_rxn_expressions[reaction_index]
        if df_index and ts_state in ts_data['ts_states']:
            reactants_list.append(json.loads(
                df_activation.reactants.loc[df_index]))
            products_list.append(json.loads(
                df_activation.products.loc[df_index]))
            neb_image_id_range = ts_data['ts_states'][ts_state]['neb_image_id_range']
            species_workfunction_data = ts_data['ts_states'][ts_state]['wf_data']
            corrected_species_workfunction_data = []
            for workfunction_value in species_workfunction_data:
                if workfunction_value != 'nan':
                    corrected_species_workfunction_data.append(workfunction_value - phi_correction)
                else:
                    corrected_species_workfunction_data.append(float('nan'))
            alk_corr.append(ts_data['alk_corr'] if beta_list[-1] else 0.0)
            charge_extrapolated_constant_potential_barriers.append(
                get_charge_extrapolated_constant_potential_barriers(
                    db_filepath, neb_image_id_range,
                    corrected_species_workfunction_data, beta_list[-1], u_she,
                    ts_data['extrapolation'], alk_corr[-1]))
        else:
            reactants_list.append(reactants_rxn_expressions[reaction_index])
            products_list.append(products_rxn_expressions[reaction_index])
            barrier_data = [float('nan'), float('nan')]
            species_barrier_fit = ts_data['barrier_fits'][ts_state]
            if 'forward' in species_barrier_fit:
                barrier_data[0] = np.poly1d(species_barrier_fit['forward'])(u_she)
            if 'backward' in species_barrier_fit:
                barrier_data[1] = np.poly1d(species_barrier_fit['backward'])(u_she)
            charge_extrapolated_constant_potential_barriers.append(tuple(barrier_data))

    for species_index, species_name in enumerate(species_list):
        if '-' in desired_surface:
            surface.append(desired_surface.split('-')[0])
        else:
            surface.append(desired_surface)
        site.append(desired_facet)
        species.append(species_name)

        raw_energy.append(float("nan"))
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

        # RHE correction
        reaction_index = species_list.index(species_name)
        df_index = df_index_rxn_expressions[reaction_index]
        if df_index:
            reactants = json.loads(df_activation_rxns.reactants.loc[df_index])
        else:
            reactants = reactants_rxn_expressions[reaction_index]
        rhe_corr.append(get_rhe_contribution(u_rhe, species_name,
                                             reference_gases, reactants,
                                             beta_list[species_index]))

        # Compute initial state energy
        ini_ads_energy_term1 = 0
        ini_ads_energy_term4 = 0
        for reactant, num_reactants in reactants_list[species_index].items():
            if 'gas' in reactant:
                noncatmap_style_species = reactant.replace('gas', '')
                idx1 = df_out.index[df_out['site_name'] == 'gas']
                idx2 = (df_out.index[df_out['species_name']
                        == noncatmap_style_species])
                idx = idx1.intersection(idx2)
                if len(idx) == 1:
                    ini_ads_energy_term1 += (
                        num_reactants * df_out.energy_vector[idx[0]][0])
                    ini_ads_energy_term4 += (
                        num_reactants * df_out.energy_vector[idx[0]][3])
            elif 'star' in reactant:
                noncatmap_style_species = reactant.replace('star', '')
                if noncatmap_style_species:
                    idx1 = df_out.index[df_out['site_name'] != 'gas']
                    idx2 = (df_out.index[df_out['species_name']
                            == noncatmap_style_species])
                    idx = idx1.intersection(idx2)
                    if len(idx) == 1:
                        ini_ads_energy_term1 += (
                            num_reactants * df_out.energy_vector[idx[0]][0])
                        ini_ads_energy_term4 += (
                            num_reactants * df_out.energy_vector[idx[0]][3])

        # Compute final state energy
        fin_ads_energy_term1 = 0
        fin_ads_energy_term4 = 0
        for product, num_products in products_list[species_index].items():
            if 'gas' in product:
                noncatmap_style_species = product.replace('gas', '')
                idx1 = df_out.index[df_out['site_name'] == 'gas']
                idx2 = (df_out.index[df_out['species_name']
                        == noncatmap_style_species])
                idx = idx1.intersection(idx2)
                if len(idx) == 1:
                    fin_ads_energy_term1 += (
                        num_products * df_out.energy_vector[idx[0]][0])
                    fin_ads_energy_term4 += (
                        num_products * df_out.energy_vector[idx[0]][3])
            elif 'star' in product:
                noncatmap_style_species = product.replace('star', '')
                if noncatmap_style_species:
                    idx1 = df_out.index[df_out['site_name'] != 'gas']
                    idx2 = (df_out.index[df_out['species_name']
                            == noncatmap_style_species])
                    idx = idx1.intersection(idx2)
                    if len(idx) == 1:
                        fin_ads_energy_term1 += (
                            num_products * df_out.energy_vector[idx[0]][0])
                        fin_ads_energy_term4 += (
                            num_products * df_out.energy_vector[idx[0]][3])

        # compute energy vector
        term1_forward_list.append(
            charge_extrapolated_constant_potential_barriers[species_index][0]
            + dft_corr[-1] + ini_ads_energy_term1)
        term1_backward = (charge_extrapolated_constant_potential_barriers[species_index][1]
                          + dft_corr[-1] + fin_ads_energy_term1)
        term2 = enthalpy[-1] + entropy[-1]
        term3 = rhe_corr[-1]
        if species_name in external_effects:
            external_effect_contribution = np.poly1d(external_effects[species_name])(u_she)
        else:
            external_effect_contribution = 0.0
        term4_forward_list.append(external_effect_contribution + ini_ads_energy_term4)
        term4_backward = external_effect_contribution + fin_ads_energy_term4
        

        G = mu = term1_backward + term2 + term3 + term4_backward
        formation_energy.append(term1_backward + term4_backward)
        energy_vector.append([term1_backward, term2, term3, term4_backward, mu, G])

        if species_name in json_species_list:
            frequencies.append(ts_vibration_data[json_species_list.index(
                species_name)]['frequencies'])
            references.append(ts_vibration_data[json_species_list.index(
                species_name)]['reference'])
        else:
            frequencies.append([])
            references.append('')

    df3 = pd.DataFrame(list(zip(surface, site, species, raw_energy,
                                charge_extrapolated_constant_potential_barriers,
                                dft_corr, zpe, enthalpy, entropy, rhe_corr,
                                formation_energy, energy_vector, frequencies,
                                references)),
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
