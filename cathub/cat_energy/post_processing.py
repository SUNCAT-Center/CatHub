'''
Module with function definitions for post-processing
'''
import numpy as np
from matplotlib import pyplot as plt

from cathub.cat_energy.conversion import read_reaction_expression_data


def get_free_energy_change_species(df, species_name):
    '''
    Compute free energy change for a given species
    '''
    if 'gas' in species_name:
        noncatmap_style_species = species_name.replace('gas', '')
        idx1 = df.index[df['site_name'] == 'gas']
        idx2 = df.index[df['species_name'] == noncatmap_style_species]
        idx = idx1.intersection(idx2)
        if len(idx) == 1:
            free_energy_change = df['energy_vector'][idx[0]][-1]
    elif 'star' in species_name:
        noncatmap_style_species = species_name.replace('star', '')
        if noncatmap_style_species:
            idx1 = df.index[df['site_name'] != 'gas']
            idx2 = df.index[df['species_name'] == noncatmap_style_species]
            idx = idx1.intersection(idx2)
            if len(idx) == 1:
                free_energy_change = df['energy_vector'][idx[0]][-1]
        else:
            free_energy_change = 0
    return free_energy_change

def get_free_energy_change_reaction(df, rxn_expression):
    '''
    Compute free energy change for a given reaction
    '''
    (reactant_dict, product_dict, _, _) = read_reaction_expression_data(
                                                                rxn_expression)
    free_energy_change_reactants = 0
    for reactant, num_reactants in reactant_dict.items():
        free_energy_change_reactants += (
                num_reactants * get_free_energy_change_species(df, reactant))

    free_energy_change_products = 0
    for product, num_products in product_dict.items():
        free_energy_change_products += (num_products
                                * get_free_energy_change_species(df, product))

    free_energy_change = (free_energy_change_products
                          - free_energy_change_reactants)
    return free_energy_change

def get_free_energy_change_rxn_mechanism(df, rxn_mechanism):
    '''
    Compute free energy change for a given reaction mechanism
    '''
    free_energy_change = 0
    for rxn_expression in rxn_mechanism:
        free_energy_change += get_free_energy_change_reaction(df,
                                                              rxn_expression)
    return free_energy_change

def plot_free_energy_diagram(df, rxn_mechanisms, labels):
    '''
    Function to plot free energy diagram for a given reaction mechanism
    '''
    # title_size = 12
    font_size = 10
    # label_size = 8
    # tick_size = 6
    line_width = 2
    # marker_size = 2
    color_list = ['#000000', '#a6611a', '#018571', '#ca0020', '#2c7bb6']
    # fig, ax = plt.subplots()
    _, ax = plt.subplots()
    for rxn_mechanism_index, rxn_mechanism in enumerate(rxn_mechanisms):
        color_value = color_list[rxn_mechanism_index]
        num_steps = len(rxn_mechanism)
        energy_data = np.zeros((num_steps, 3))
        for step_index, rxn_expression in enumerate(rxn_mechanism):
            (reactant_dict, product_dict,
             ts_state, _) = read_reaction_expression_data(rxn_expression)
            for reactant, num_reactants in reactant_dict.items():
                energy_data[step_index, 0] += (
                                num_reactants
                                * get_free_energy_change_species(df, reactant))
            for product, num_products in product_dict.items():
                energy_data[step_index, 1] += (
                    num_products * get_free_energy_change_species(df, product))
            if ts_state:
                ts_state += 'star'
                energy_data[step_index, 2] += get_free_energy_change_species(
                                                                df, ts_state)
            if step_index == 0:
                ax.plot([step_index + 0.75, step_index + 1.25],
                        [energy_data[step_index, 0]] * 2, color=color_value,
                        linewidth=line_width, label=labels[rxn_mechanism_index])
            ax.plot([step_index + 1.75, step_index + 2.25],
                    [energy_data[step_index, 1]] * 2, linewidth=line_width,
                    color=color_value)
            if energy_data[step_index, 2]:
                p = np.polyfit(
                    [step_index + 1.25, step_index + 1.50, step_index + 1.75],
                    [energy_data[step_index, 0], energy_data[step_index, 2],
                     energy_data[step_index, 1]], 2)
                x_data = np.linspace(step_index + 1.25, step_index + 1.75, 100)
                y_data = (np.poly1d(p))(x_data)
                ax.plot(x_data, y_data, linewidth=line_width, color=color_value)
            else:
                ax.plot(
                    [step_index + 1.25, step_index + 1.75],
                    [energy_data[step_index, 0], energy_data[step_index, 1]],
                    linewidth=line_width, color=color_value)

    figure_name = 'Free Energy Diagram.png'
    plt.xlabel('Step', fontsize=font_size)
    plt.ylabel('Energy (eV)', fontsize=font_size)
    plt.legend()
    plt.tight_layout()
    plt.savefig(figure_name, dpi=600)
    return None
