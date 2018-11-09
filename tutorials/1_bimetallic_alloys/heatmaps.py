import sys
import json
import re
import numpy as np
import matplotlib
import pylab as p
from ase.db import connect
from ase.visualize import view
import matplotlib.transforms

from cathub.query import get_reactions

from tools import ordered_metals, get_AB_from_formula, references, site2int, site_labels

"""
Tutorial 1.1

Plotting adsorption energy - and site - heatmaps for bimetallic alloys.

Choose the adsorbate and site below.

Adsorbates includes: H, N, C, O, S, CH, NH, CH2, CH3, OH and H2O

sites can be:
'~', '~hollow', '~bridge, '~top', hollow|A_A_A|HCP, hollow|A_A_A|FCC, hollow|A_A_B|HCP, hollow|A_A_B|FCC, 
bridge|A_A|A, bridge|A_A|B, bridge|A_B|A, bridge|B_B|B, top|A, top|B

Also you will have to specift the Structure Bericht symbol (SB_symbol) of the alloy, which can be L12 or L10

"""

adsorbate = 'C'
adsorption_site = '~'
SB_symbol = 'L12' # 'L10' or 'L12'

data = get_reactions(n_results='all',
                     pubId='WintherOpen2018',
                     sites = adsorption_site,
                     reactants=references[adsorbate],
                     products=adsorbate,
                     columns=['surfaceComposition, reactionEnergy', 'sites', 'products'])

data = data['reactions']
totalCount = data['totalCount']
edges = data['edges']

site_points = (np.array(range(1, 14)) - 0.5) * 20 / 12



dim = len(ordered_metals)
EADS = np.zeros([dim, dim])
SITES = np.zeros([dim, dim])
EADS.fill(None)
SITES.fill(None)
for edge in edges:
    result = edge['node']

    adsorbates = list(json.loads(result['products']).keys())
    prefactor_adsorbate = list(json.loads(result['products']).values())[0]

    # Only include results with one adsorbate
    if len(adsorbates) > 1 or prefactor_adsorbate > 1:  
        continue
    
    formula = result['surfaceComposition']
    E = result['reactionEnergy']
    sites = json.loads(result['sites'])
    site = list(sites.values())[0]
    if 'tilt' in site:
        continue
    site = site2int(site)
    
    A, B, SB = get_AB_from_formula(formula)

    if SB != 'A1' and SB != SB_symbol:
        continue

    iA = ordered_metals.index(A)
    iB = ordered_metals.index(B)

    if np.isnan(EADS[iA, iB]) or EADS[iA, iB] > float(E):
        EADS[iA, iB] = E
        SITES[iA, iB] = site
        if SB_symbol == 'L10':
            EADS[iB, iA] = E
            SITES[iB, iA] = site



""" Plot adsorption energies"""
fig = p.figure(figsize=(6, 6))
cmap = matplotlib.cm.cool
cmap.set_bad('gray', 1.)
p.imshow(EADS, cmap=cmap)

metalstrx = ['\n{}'.format(o) if i in range(0,37,2) else o for i, o in enumerate(ordered_metals)]
metalstry = ['     {}'.format(o) if i in range(1,37,2) else o for i, o in enumerate(ordered_metals)]
p.xticks(range(dim), metalstrx)
p.yticks(range(dim), metalstry)

offset = matplotlib.transforms.ScaledTranslation(-0.45, 0, fig.dpi_scale_trans)

ax = p.gca()
for label in ax.yaxis.get_majorticklabels():
    label.set_horizontalalignment('left')
    label.set_transform(label.get_transform() + offset)

p.xlabel('Metal B')
p.ylabel('Metal A')
p.colorbar(fraction=0.046, pad=0.035).set_label(label='E$_\mathrm{ads}$ (eV)', size=12)
p.title('Adsorption energies: {} @ {}$_{}$'.format(adsorbate, SB_symbol[:-1], SB_symbol[-1]))

p.subplots_adjust(left=0.15, bottom=0.1, top=0.9)

p.savefig('{}_{}_{}_adsorption_energies.pdf'.format(adsorption_site, adsorbate, SB_symbol))    
p.show()

""" Plot adsorption sites"""
p.figure(figsize=(7,6))
cmap = matplotlib.cm.tab20c
import matplotlib.colors as colors
cmap = colors.ListedColormap(cmap(list(range(4,6)) + list(range(8, 18))))

cmap.set_bad('white',1.)
p.imshow(SITES, vmin= 0, vmax = 20, cmap=cmap)

p.xticks(range(dim), metalstrx)
p.yticks(range(dim), metalstry)

offset = matplotlib.transforms.ScaledTranslation(-0.45, 0, fig.dpi_scale_trans)

ax = p.gca()
for label in ax.yaxis.get_majorticklabels():
    label.set_horizontalalignment('left')
    label.set_transform(label.get_transform() + offset)

p.xlabel('Metal B')
p.ylabel('Metal A')
p.title('Adsorption sites: {} @ {}$_{}$'.format(adsorbate,  SB_symbol[:-1], SB_symbol[-1]))

cbar = p.colorbar(ticks=site_points, fraction=0.046, pad=0.035)
cbar.set_label(label='site', size=12)
cbar.ax.set_yticklabels(site_labels, size=12)
p.subplots_adjust(right=0.75, bottom=0.1, top=0.9)

p.savefig('{}_{}_{}_adsorption_sites.pdf'.format(adsorption_site, adsorbate, SB_symbol))
p.show()
