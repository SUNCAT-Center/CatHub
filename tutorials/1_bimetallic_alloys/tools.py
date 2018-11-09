import re

""" Modified Pettifor scale """
ordered_metals = ['Y', 'La', 'Sc',  # Group 3 metals
                  'Zr', 'Hf', 'Ti', 'Ta', 'Nb',  # D band metals
                  'V',  'Cr', 'Mo', 'W', 'Re',
                  'Tc', 'Os', 'Ru', 'Ir','Rh',
                  'Ni', 'Co', 'Fe', 'Mn', 'Pt', 'Pd',
                  'Au', 'Ag', 'Cu',  # Filled d band 
                  'Zn', 'Cd', 'Hg', 'Al', 'Ga', 'In',
                  'Tl', 'Pb', 'Sn', 'Bi']   # Post transition metals


references = {}
references['H'] = 'H2'
references['N'] = 'N2'
references['C'] = 'CH4+H2'
references['O'] = 'H2O+H2'
references['S'] = 'H2S+H2'
references['CH'] = 'CH4+H2'
references['NH'] = 'N2+H2'
references['CH2'] = 'CH4+H2'
references['CH3'] = 'CH4+H2'
references['OH'] = 'CH4+H2'
references['H2O'] = 'CH4+H2'


def get_AB_from_formula(formula):
    if '3' in formula:
        A, B = formula.split('3')
        return A, B, 'L12'

    AB = sorted(re.findall('([A-Z][^A-Z]*)', formula))
    if len(AB) > 1:
        A, B = AB
        return A, B, 'L10'
    else:
        A = AB[0]
        B = A
        return A, B, 'A1'


def site2int(site):
    a = None
    if site == 'top' or site == 'top|A':
        a=1
    elif site == 'top|B':
        a=2    
    elif site == 'bridge|A_A|A':
        a=3
    elif site == 'bridge|A_A|B':
        a=4
    elif site == 'bridge|A_B|A':
        a=5
    elif site == 'bridge|B_B|B':
        a=6       
    elif site == 'hollow|A_A_A|HCP' or site == 'HCP':
        a=7
    elif site == 'hollow|A_A_A|FCC' or site == 'FCC':
        a=8
    elif site == 'hollow|A_A_B|HCP' or site == 'hollow|A_B_B|HCP':
        a=9
    elif site == 'hollow|A_A_B|FCC' or site == 'hollow|A_B_B|FCC':
        a=10
    elif '4fold' in site:
        a=11
    elif site == 'Subsurface|':
        a=12
    if not a:
        return a
    return (a-0.5) * 20 / 12 #+ 6 #  6

site_labels = ['top|A',
               'top|B',
               'bridge AA|A',
               'bridge AA|B',
               'bridge AB|A',
               'bridge BB|B',
               'hollow AAA|HCP',
               'hollow AAA|FCC',
               'hollow AAB|HCP',
               'hollow AAB|FCC',
               '4fold',
               'subsurface']
