import json
import numpy as np
import pylab as p
from cathub.query import get_reactions
from tools import ordered_metals, get_AB_from_formula, references


""" 
Tutorial 1.2

Plotting adsorption energy relation (scaling) plot for bimetallic alloys.

Choose adsorbate1, adsorbate2 and site

Adsorbates includes: H, N, C, O, S, CH, NH, CH2, CH3, OH and H2O

sites can be:
'~', '~hollow', '~bridge, '~top', hollow|A_A_A|HCP, hollow|A_A_A|FCC, hollow|A_A_B|HCP, hollow|A_A_B|FCC, 
bridge|A_A|A, bridge|A_A|B, bridge|A_B|A, bridge|B_B|B, top|A, top|B
"""

adsorbate1 = 'C'
adsorbate2 = 'O'
adsorption_site = '~'  # get all sites

edges = []
for adsorbate in [adsorbate1, adsorbate2]:
    data = get_reactions(n_results='all',
                         pubId='WintherOpen2018',
                         sites=adsorption_site,
                         reactants=references[adsorbate],
                         products=adsorbate,
                         columns=['surfaceComposition',
                                  'products',
                                  'reactionEnergy',
                                  'sites'])

    data = data['reactions']
    totalCount = data['totalCount']

    edges += data['edges']
    
dim = len(ordered_metals)
E_ads =  {'L12': {},
          'L10': {},
          'A1': {}}

for ads in [adsorbate1, adsorbate2]:
    for SB_symbol in ['L12', 'L10', 'A1']:
        E_ads[SB_symbol].update({ads: np.zeros([dim, dim])})
        E_ads[SB_symbol][ads].fill(None)

    
for edge in edges:
    result = edge['node']

    adsorbates = list(json.loads(result['products']).keys())
    prefactor_adsorbate = list(json.loads(result['products']).values())[0]

    # Only include results with one adsorbate
    if len(adsorbates) > 1 or prefactor_adsorbate > 1:  
        continue
    formula = result['surfaceComposition']
    adsorbate = adsorbates[0].replace('star', '') 
    e_ads = result['reactionEnergy']
    sites = json.loads(result['sites'])
    
    # get metal A, B and SB_symbol from formula
    A, B, SB_symbol = get_AB_from_formula(formula)
    
    iA = ordered_metals.index(A)
    iB = ordered_metals.index(B)

    if np.isnan(E_ads[SB_symbol][adsorbate][iA, iB]) or \
    E_ads[SB_symbol][adsorbate][iA, iB] > float(e_ads):
        E_ads[SB_symbol][adsorbate][iA, iB] = e_ads


""" Plot results """

d_metal_index = np.array(range(3, 24))
non_d_metal_index = np.array([i for i in range(37) if i not in d_metal_index])

p.figure(figsize=(6, 6))

p.scatter(E_ads['L12'][adsorbate1],
          E_ads['L12'][adsorbate2], c='skyblue', label='non-d alloys')
p.scatter(E_ads['L10'][adsorbate1],
          E_ads['L10'][adsorbate2], c='skyblue')#, label='non-d alloys')


p.scatter(E_ads['L12'][adsorbate1]\
          [tuple([d_metal_index])].T[tuple([d_metal_index])].T,
          E_ads['L12'][adsorbate2]\
          [tuple([d_metal_index])].T[tuple([d_metal_index])].T,
          c='lightseagreen', label='d alloys')
p.scatter(E_ads['L10'][adsorbate1]\
          [tuple([d_metal_index])].T[tuple([d_metal_index])].T,
          E_ads['L10'][adsorbate2]\
          [tuple([d_metal_index])].T[tuple([d_metal_index])].T,
          c='lightseagreen')


p.scatter(E_ads['A1'][adsorbate1].diagonal()[non_d_metal_index],
          E_ads['A1'][adsorbate2].diagonal()[non_d_metal_index], c='magenta',
          label='non-d metals')


p.scatter(E_ads['A1'][adsorbate1].diagonal()[d_metal_index],
          E_ads['A1'][adsorbate2].diagonal()[d_metal_index],
          c='orange', label='d metals')


p.legend()

p.title('{} vs {} adsorption energy'.format(adsorbate1, adsorbate2, va='bottom'))

p.xlabel('E$_\mathrm{{ads}}$({}) (eV)'.format(adsorbate1))
p.ylabel('E$_\mathrm{{ads}}$({}) (eV)'.format(adsorbate2))

for i, txt in enumerate(np.array(ordered_metals)):
    color = 'magenta'
    if i in d_metal_index:
        color = 'orange'
    sgn1 = 0.4
    sgn2 = 0.5 
    if txt in ['Ir', 'Sc', 'Nb']:
        sgn1 = -2.5
        sgn2 = -1.9
    elif txt in ['Pb','Sn', 'Ru', 'La', 'Y']:
        sgn1 = 0.4
        sgn2 = -2 
    p.annotate(txt,
               xy=(E_ads['A1'][adsorbate1][i,i], E_ads['A1'][adsorbate2][i,i]),
               xytext=(sgn1 * 5, sgn2 * 5),
               textcoords='offset points',
               color=color)

p.savefig('scaling_{}_{}.pdf'.format(adsorbate1, adsorbate2))
p.show()
