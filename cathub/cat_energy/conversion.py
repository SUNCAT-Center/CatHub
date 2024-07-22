'''
Module to conversion functions
'''
import re

import numpy as np
from ase.data import chemical_symbols
from ase.units import _e, _hplanck, _c


CM2M = 1E-02
CM2EV = _hplanck / _e * _c / CM2M
KB = 8.617333262145e-5
PHI_REF = 4.44  # experimental work function of SHE relative to vacuum at  0 V


def read_reaction_expression_data(rxn_expression):
    '''
    Read reaction expression data and return its components.

    This function parses a reaction expression string to extract its reactants, 
    products, transition state, and beta value. It handles both reversible and 
    irreversible reactions and processes the species names to conform to CatMAP 
    style.

    Parameters:
    -----------
    rxn_expression : str
        The reaction expression string to be parsed. It may include a beta value 
        separated by a semicolon.

    Returns:
    --------
    tuple
        A tuple containing:
        - reactant_dict (dict): Dictionary of reactant species with their quantities.
        - product_dict (dict): Dictionary of product species with their quantities.
        - ts_state (str): Transition state species, if present.
        - beta (float): Beta value for the reaction, default is 0.0 if not specified.
    '''
    discard_species_list = ['*_t', '_t', '_g']
    reactant_discard_species_list = ['_t', '_g']
    if ';' in rxn_expression:
        (rxn, beta) = rxn_expression.split(';')
    else:
        rxn = rxn_expression
        beta = 0.0  # chemical reaction if beta not specified

    if isinstance(beta, str):
        if '=' in beta:
            beta = float(beta.split('=')[1])

    rxn_no_spaces = re.sub(' ', '', rxn)
    split_rxn = re.split('->|<->', rxn_no_spaces)

    reactant_term = split_rxn[0]
    product_term = split_rxn[-1]

    if '+' in reactant_term:
        reactant_species = reactant_term.split('+')
        for discard_species in reactant_discard_species_list:
            if discard_species in reactant_species:
                reactant_species.remove(discard_species)
        reactant_list = reactant_species
    else:
        reactant_list = [reactant_term]

    reactant_dict = {}
    for reactant in reactant_list:
        if 'ele_g' not in reactant:
            (catmap_species, num_species) = get_catmap_style_species(reactant)
            if catmap_species in reactant_dict:
                reactant_dict[catmap_species] += num_species
            else:
                reactant_dict[catmap_species] = num_species

    if '+' in product_term:
        product_species = product_term.split('+')
        for discard_species in discard_species_list:
            if discard_species in product_species:
                product_species.remove(discard_species)
        product_list = product_species
    else:
        product_list = [product_term]

    product_dict = {}
    for product in product_list:
        if 'ele_g' not in product:
            (catmap_species, num_species) = get_catmap_style_species(product)
            if catmap_species in product_dict:
                product_dict[catmap_species] += num_species
            else:
                product_dict[catmap_species] = num_species

    if rxn_no_spaces.count('->') == 2:
        ts_term = split_rxn[1]
        if '+' in ts_term:
            ts_species = ts_term.split('+')
            for discard_species in discard_species_list:
                if discard_species in ts_species:
                    ts_species.remove(discard_species)
            ts_state = ts_species[0]
        else:
            ts_state = ts_term
    else:
        ts_state = ''

    for discard_species in discard_species_list:
        ts_state = ts_state.replace(discard_species, '')
    return (reactant_dict, product_dict, ts_state, beta)


def get_catmap_style_species(species):
    '''
    Return species name in CatMAP style.

    This function converts a species name to the CatMAP style format and determines 
    the number of species present. It handles different formats of species names, 
    including gas and adsorbed species.

    Parameters:
    -----------
    species : str
        The species name to be converted.

    Returns:
    --------
    tuple
        A tuple containing:
        - catmap_species (str): The species name in CatMAP style.
        - num_species (float): The number of species.
    '''
    if 'H_g' in species:
        split_species = re.split('(\\d+)', species)
        if '' in split_species:
            split_species.remove('')
        catmap_species = 'H2gas'
        if len(split_species) > 1:
            num_species = float(split_species[0]) * 0.5
        else:
            num_species = 0.5
    else:
        if species[0].isnumeric():
            split_species = re.split('(\\d+)', species)
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


def populate_chemical_symbols_dict(chemical_symbols_dict, last_chemical_symbol):
    '''
    Populate dictionary of chemical symbols.

    This function updates a dictionary of chemical symbols by incrementing the count 
    of a specified chemical symbol. If the chemical symbol is not already present 
    in the dictionary, it is added with a count of one.

    Parameters:
    -----------
    chemical_symbols_dict : dict
        Dictionary containing chemical symbols as keys and their counts as values.
    last_chemical_symbol : str
        The chemical symbol to be added or incremented in the dictionary.

    Returns:
    --------
    dict
        The updated dictionary with the specified chemical symbol's count incremented.
    '''
    if last_chemical_symbol in chemical_symbols_dict:
        chemical_symbols_dict[last_chemical_symbol] += 1
    else:
        chemical_symbols_dict[last_chemical_symbol] = 1
    return chemical_symbols_dict


def formula_to_chemical_symbols(formula):
    '''
    Return dictionary mapping chemical symbols to number of atoms.

    This function parses a chemical formula string and returns a dictionary 
    mapping each chemical symbol to the number of atoms present in the formula.

    Parameters:
    -----------
    formula : str
        The chemical formula as a string.

    Returns:
    --------
    dict
        Dictionary containing chemical symbols as keys and their counts as values.
    '''
    chemical_symbols_dict = {}

    # split chemical formula string into alpha and numeric characters
    regex = re.compile('(\\d+|\\s+)')
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
                chemical_symbols_dict = populate_chemical_symbols_dict(
                                    chemical_symbols_dict, last_chemical_symbol)

    # multiply number of atoms for each chemical symbol with
    # number of formula units
    for key in chemical_symbols_dict:
        chemical_symbols_dict[key] = (formula_unit_count
                                      * chemical_symbols_dict[key])
    return chemical_symbols_dict


def u_she_to_u_rhe(u_she, temp, pH):
    '''
    Convert SHE potential to RHE potential.

    This function converts the Standard Hydrogen Electrode (SHE) potential to 
    the Reversible Hydrogen Electrode (RHE) potential based on temperature and pH.

    Parameters:
    -----------
    u_she : float
        Standard Hydrogen Electrode potential.
    temp : float
        Temperature in Kelvin.
    pH : float
        pH value of the system.

    Returns:
    --------
    float
        The converted Reversible Hydrogen Electrode potential.
    '''
    u_rhe = u_she + KB * temp * np.log(10) * pH
    return u_rhe


def u_she_to_field(u_she, u_pzc, d):
    '''
    Convert SHE potential to electric field.

    This function converts the Standard Hydrogen Electrode (SHE) potential to 
    the electric field based on the potential of zero charge (PZC) and the 
    distance over which the potential difference is applied.

    Parameters:
    -----------
    u_she : float
        Standard Hydrogen Electrode potential.
    u_pzc : float
        Potential of zero charge.
    d : float
        Distance between the positively (center of the Helmholtz plane) and 
        negatively charged planes.

    Returns:
    --------
    float
        The computed electric field.
    '''
    field = (u_she - u_pzc) / d
    return field


def u_rhe_to_field(u_rhe, pH, u_pzc, d, temp):
    '''
    Convert RHE potential to electric field.

    This function converts the Reversible Hydrogen Electrode (RHE) potential 
    to the Standard Hydrogen Electrode (SHE) potential and then calculates 
    the electric field based on the potential of zero charge (PZC) and the 
    distance between the positively (center of the Helmholtz plane) and 
    negatively charged planes.

    Parameters:
    -----------
    u_rhe : float
        Reversible Hydrogen Electrode potential.
    pH : float
        pH value of the system.
    u_pzc : float
        Potential of zero charge.
    d : float
        Distance between the positively (center of the Helmholtz plane) and 
        negatively charged planes.
    temp : float
        Temperature in Kelvin.

    Returns:
    --------
    tuple
        A tuple containing:
        - u_she (float): The converted Standard Hydrogen Electrode potential.
        - field (float): The computed electric field.
    '''
    u_she = u_rhe - KB * temp * np.log(10) * pH
    field = (u_she - u_pzc) / d
    return (u_she, field)


def field_to_voltage(field, pH, u_pzc, d, temp):
    '''
    Convert electric field to SHE/RHE potential.

    This function converts an electric field to the Standard Hydrogen Electrode (SHE) 
    potential and the Reversible Hydrogen Electrode (RHE) potential based on the 
    potential of zero charge (PZC), the distance between the positively (center of the 
    Helmholtz plane) and negatively charged planes, the pH value, and the temperature.

    Parameters:
    -----------
    field : float
        The electric field.
    pH : float
        pH value of the system.
    u_pzc : float
        Potential of zero charge.
    d : float
        Distance between the positively (center of the Helmholtz plane) and 
        negatively charged planes.
    temp : float
        Temperature in Kelvin.

    Returns:
    --------
    tuple
        A tuple containing:
        - u_she (float): The Standard Hydrogen Electrode potential.
        - u_rhe (float): The Reversible Hydrogen Electrode potential.
    '''
    u_she = field * d + u_pzc
    u_rhe = u_she + KB * temp * np.log(10) * pH
    return (u_she, u_rhe)


def get_rhe_contribution(u_rhe, species_value, reference_gases,
                         reactants=None, beta=None):
    '''
    Compute the RHE-scale contribution towards free energy.

    This function calculates the RHE-scale contribution to the free energy change 
    based on the Reversible Hydrogen Electrode (RHE) potential, the species being 
    analyzed, and the reference gases. It considers whether the species is a 
    reactant or a product in the reaction mechanism.

    Parameters:
    -----------
    u_rhe : float
        Reversible Hydrogen Electrode potential.
    species_value : str
        The species for which the RHE contribution is being computed.
    reference_gases : list of str
        List of reference gas species.
    reactants : dict, optional
        Dictionary of reactant species with their quantities (default is None).
    beta : float, optional
        Beta value for the reaction (default is None).

    Returns:
    --------
    float
        The computed RHE-scale contribution to the free energy.
    '''
    # x CO + y (H++e-) = CxH(y-2x+2z)Oz + (x-z) H2O
    # x CO + y/2 H2 = CxH(y-2x+2z)Oz + (x-z) H2O
    # Based on computational hydrogen electrode, n should be twice the number
    # of H2 gas molecules that are required for the reduction reaction
    if reactants:
        if 'H2gas' in reactants:
            n = 2 * reactants['H2gas']
        else:
            n = 0
    else:
        if '_ref' in species_value or species_value == 'H2':
            chemical_symbols_dict = {}
        else:
            chemical_symbols_dict = formula_to_chemical_symbols(species_value)

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
            n_h2_gas = (2 * x - z + y / 2)
        # CO Reduction Reaction
        elif set(reference_gases) == set(['CO', 'H2_ref', 'H2O']):
            n_h2_gas = (x - z + y / 2)

        n = 2 * n_h2_gas

    # u_rhe derived from epsilon
    # zero-charge potential of the metal electrode in V vs. SHE for Cu(100)
    # UM_PZC = -0.54 V
    # d = 1.2  # thickness in A (angstrom)
    # u_she = UM_PZC + d * epsilon
    # u_rhe = u_she + 0.059 * pH

    if beta:
        rhe_energy_contribution = - (1 - beta) * u_rhe
    else:
        rhe_energy_contribution = n * u_rhe

    return rhe_energy_contribution


def get_electric_field_contribution(
        field_effects, species_value, reference_gases, reactants=None,
        beta=None):
    '''
    Compute the contribution of electric field towards the free energies of adsorption energies.

    This function calculates the contribution of the electric field to the free energies 
    of adsorption for a given species. It considers both RHE-scale and SHE-scale dependencies.

    Parameters:
    -----------
    field_effects : dict
        Dictionary containing field effect parameters such as:
            'she_voltage' (float): Standard Hydrogen Electrode voltage.
            'pzc_voltage' (float): Potential of zero charge.
            'd' (float): Distance between the positively (center of the Helmholtz plane) 
                         and negatively charged planes.
            'U_RHE' (float): Reversible Hydrogen Electrode potential.
            'mu' (dict): Dipole moment values for species.
            'alpha' (dict): Polarizability values for species.
    species_value : str
        The species for which the electric field contribution is being computed.
    reference_gases : list of str
        List of reference gas species.
    reactants : dict, optional
        Dictionary of reactant species with their quantities (default is None).
    beta : float, optional
        Beta value for the reaction (default is None).

    Returns:
    --------
    tuple
        A tuple containing:
        - rhe_energy_contribution (float): The RHE-scale contribution to the free energy.
        - she_energy_contribution (float): The SHE-scale contribution to the free energy.
    '''
    she_voltage = field_effects['she_voltage']
    pzc_voltage = field_effects['pzc_voltage']
    d = field_effects['d']
    u_rhe = field_effects['U_RHE']
    mu = field_effects['mu']
    alpha = field_effects['alpha']

    # u_rhe-scale dependency
    rhe_energy_contribution = get_rhe_contribution(
                        u_rhe, species_value, reference_gases, reactants, beta)

    # u_she-scale dependency
    if reactants is None:
        if species_value + '_g' in mu:
            species_value = species_value + '_g'
    if species_value in mu:
        epsilon = (she_voltage - pzc_voltage) / d
        she_energy_contribution = (mu[species_value] * epsilon
                                     - 0.5 * alpha[species_value] * epsilon**2)
    else:
        she_energy_contribution = 0.0

    return (rhe_energy_contribution, she_energy_contribution)
