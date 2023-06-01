#!/usr/bin/env python

from cathub.cathubsql import CathubSQL
import numpy as np
import matplotlib.pyplot as plt
import yaml


def plot_vshe_energy(surface, facet, product, site, vshe_reaction_dict):
    data_dict = vshe_reaction_dict[surface][facet][product][site]
    vshe_values = data_dict['vshe']
    reaction_energies = data_dict['energy']

    # Fit a polynomial curve
    degree = 2  # Set the degree of the polynomial
    coefficients = np.polyfit(vshe_values, reaction_energies, degree)
    poly = np.poly1d(coefficients)

    # Calculate R-squared value
    y_mean = np.mean(reaction_energies)
    ss_total = np.sum((reaction_energies - y_mean) ** 2)
    ss_residual = np.sum((reaction_energies - poly(vshe_values)) ** 2)
    r_squared = 1 - (ss_residual / ss_total)

    # Generate a range of values for the x-axis
    x_range = np.linspace(min(vshe_values), max(vshe_values), 100)

    # Plot the data points and the fitted curve
    color = 'blue'
    plt.plot(vshe_values, reaction_energies, 'o', color=color, label='Data')
    plt.plot(x_range, poly(x_range), color=color, label='Curve Fit')

    # Determine the coordinates for the annotation
    x_min, x_max = plt.xlim()
    y_min, y_max = plt.ylim()
    x_pos = x_max - 0.05 * (x_max - x_min)
    y_pos = y_min + 0.05 * (y_max - y_min)

    # Add the R-squared value as a text annotation on the plot
    plt.text(x_pos, y_pos, f'R-squared = {r_squared:.2f}', ha='right', va='bottom')

    plt.xlabel('$V_{SHE}$')
    plt.ylabel('Reaction Energy')
    plt.title(f'{product}@{site} on {surface}({facet})')
    plt.show()


def print_available_sites(surface, facet, product, vshe_reaction_dict):
    available_sites = list(vshe_reaction_dict[surface][facet][product].keys())
    print(f"Available Sites for Facet: {facet}, Product: {product}: {available_sites}")
    print()


# Convert tuples to lists in vshe_reaction_dict
def convert_tuples_to_lists(data):
    if isinstance(data, dict):
        for key, value in data.items():
            data[key] = convert_tuples_to_lists(value)
    elif isinstance(data, tuple):
        data = list(data)
    return data


def generate_yaml_from_pub_id(pub_id_name):
    pub_id_name = 'PasumarthiFacetDependence2023'

    # To get data on catalysis-hub.org
    db = CathubSQL()

    # Get reactions in pandas dataframe
    dataframe = db.get_dataframe(pub_id=pub_id_name)

    vshe_reaction_dict = {}
    surface_list = dataframe['surface_composition'].unique()
    for surface in surface_list:
        vshe_reaction_dict[surface] = {}  # Create surface-specific dictionary
        surface_rows = dataframe[dataframe['surface_composition'] == surface]
        facet_list = dataframe['facet'].unique()
        for facet in facet_list:
            vshe_reaction_dict[surface][facet] = {}  # Create facet-specific dictionary

            facet_rows = surface_rows[surface_rows['facet'] == facet]
            unique_products = list({list(row.keys())[0] for row in facet_rows['products']})

            # Convert 'products' column to strings and fill NaN values in-place
            facet_rows.loc[:, 'products'] = facet_rows['products'].astype(str).fillna({i: {} for i in facet_rows.index})
            for product in unique_products:
                vshe_reaction_dict[surface][facet][product] = {}  # Create product-specific dictionary

                product_rows = facet_rows[facet_rows['products'].str.contains(r'\b{}\b'.format(product), regex=True)]
                unique_sites = list({list(row.values())[0] for row in product_rows['sites']})

                for site in unique_sites:
                    vshe_reaction_dict[surface][facet][product][site] = {'vshe': [], 'energy': []}  # Create site-specific dictionary with empty lists

                    site_rows = product_rows[product_rows['sites'].apply(lambda x: x.get(list(x.keys())[0])) == site]

                    vshe_values = [float(value) for value in site_rows['dft_functional'].str.extract(r'(-?\d+\.\d+)VSHE', expand=False)]
                    reaction_energies = list(site_rows['reaction_energy'])

                    sorted_pairs = sorted(zip(vshe_values, reaction_energies))
                    vshe_reaction_dict[surface][facet][product][site]['vshe'], vshe_reaction_dict[surface][facet][product][site]['energy'] = zip(*sorted_pairs)

    # Convert tuples to lists in vshe_reaction_dict
    vshe_reaction_dict = convert_tuples_to_lists(vshe_reaction_dict)

    # Save the vshe_reaction_dict as a YAML file
    with open('vshe_reaction_dict.yaml', 'w') as file:
        yaml.dump(vshe_reaction_dict, file)


def dict_to_yaml_str(data, indent=0):
    lines = []
    for key, value in sorted(data.items()):
        if isinstance(value, dict):
            sub_lines = dict_to_yaml_str(value, indent + 2).split('\n')
            lines.append(f"{' ' * indent}{key}:")
            lines.extend([f"{' ' * (indent + 2)}{line}" for line in sub_lines])
        else:
            lines.append(f"{' ' * indent}{key}: {str(value)}")
    return '\n'.join(lines)


def analyze_vshe_reaction_dict(vshe_reaction_dict):
    analysis_dict = {}
    surface_list = list(vshe_reaction_dict.keys())
    for surface in surface_list:
        analysis_dict[surface] = {}
        facet_list = list(vshe_reaction_dict[surface].keys())
        for facet in facet_list:
            analysis_dict[surface][facet] = {}
            adsorbate_list = list(vshe_reaction_dict[surface][facet].keys())
            for adsorbate in adsorbate_list:
                site_list = list(vshe_reaction_dict[surface][facet][adsorbate].keys())
                analysis_dict[surface][facet][adsorbate] = site_list

    yaml_str = dict_to_yaml_str(analysis_dict)

    # Save the vshe_reaction_dict as a YAML file
    with open('available_site_info.yaml', 'w') as file:
        file.write(yaml_str)


def main():
    pub_id_name = 'PasumarthiFacetDependence2023'
    generate_yaml_from_pub_id(pub_id_name)

    # Load the YAML file
    with open('vshe_reaction_dict.yaml', 'r') as file:
        vshe_reaction_dict = yaml.safe_load(file)

    # analyzes to generate a yaml with list of available surfaces, facets, products, and sites
    analyze_vshe_reaction_dict(vshe_reaction_dict)

    # Print available sites for a given facet and product
    surface = 'Cu'
    facet = '100'
    product = 'COHstar'
    print_available_sites(surface, facet, product, vshe_reaction_dict)

    # Plot example
    surface = 'Cu'
    facet = '100'
    product = 'OCCOHstar'
    site = '4foldhollow'
    plot_vshe_energy(surface, facet, product, site, vshe_reaction_dict)

if __name__ == '__main__':
    main()
