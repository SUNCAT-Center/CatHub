#!/usr/bin/env python

# builtin imports
from .ase_tools import gas_phase_references, get_chemical_formula, \
    get_reduced_chemical_formula, symbols, collect_structures, \
    compare_parameters
import cathub.ase_tools
import ase.atoms
import ase.utils
import ase.io
import numpy as np
import pickle
import json
import yaml
from yaml import Dumper
import re
import pprint
import collections
from random import randint
from pathlib import Path
Path().expanduser()

# local imports


def fuzzy_match(structures, options):
    # filter out cell with ill-defined unit cells
    structures = [structure[-1] for structure in structures
                  if structure[-1].number_of_lattice_vectors == 3
                  ]
    # sort by density
    structures = sorted(structures,
                        key=lambda x: len(x) / x.get_volume()
                        )
    adsorbate_numbers = [sorted(ase.atoms.Atoms(ads).numbers)
                         for ads in options.adsorbates]  # .split(',')]
    # group in to bulk, surface, or bulk
    molecules, surfaces, bulks = [], [], []
    gas_phase_candidates = []
    reference_energy = {}
    collected_energies = {}
    key_count = {}
    collected_structures = {}
    check_parameters = not options.skip_parameters

    if options.verbose:
        print("\nGroup By Densities")
        print("===================")
    for structure in structures:
        # add more info from filename
        facet_match = re.search(
            '(?<=[^0-9])?[0-9]{3,3}(?=[^0-9])?', structure.info['filename'])
        site_match = [site_name for site_name in
                      ['top', 'bridge', 'hollow'] if site_name in structure.info['filename']]

        if facet_match and options.facet_name == 'facet':
            structure.info['facet'] = facet_match.group()
        else:
            structure.info['facet'] = options.facet_name or ''
        if site_match:
            structure.info['site'] = site_match[0]

        density = len(structure) / structure.get_volume()
        if options.verbose:
            print("  {density:7.3f} {filename}".format(
                density=density,
                filename=structure.info['filename'],
            ))
        if density < options.max_density_gas:
            structure.info['state'] = 'molecule'
            molecules.append(structure)
            collected_structures \
                .setdefault(options.dft_code
                            or structure.info['filetype'], {}) \
                .setdefault(options.xc_functional or 'XC_FUNCTIONAL', {}) \
                .setdefault('gas', {}) \
                .setdefault(get_chemical_formula(structure), structure)

            formula = get_chemical_formula(structure)
            if formula not in options.exclude_reference.split(','):
                gas_phase_candidates.append(get_chemical_formula(structure))
                reference_energy[formula] = structure.get_potential_energy()
                if formula in options.energy_corrections.keys():
                    reference_energy[formula] += \
                        options.energy_corrections[formula]
                if options.verbose:
                    print("        GAS", formula,
                          structure.info['filename'])

        elif density < options.max_density_slab:
            structure.info['state'] = 'surface'
            formula = get_chemical_formula(structure)
            surfaces.append(structure)
            if options.verbose:
                print("        SURFACE", formula,
                      structure.info['filename'])
        else:
            structure.info['state'] = 'bulk'
            bulks.append(structure)

            if options.verbose:
                print("        BULK", formula, structure.info['filename'])

    # sort surfaces by volume to get difference facets
    surfaces = sorted(surfaces,
                      key=lambda x: x.get_volume()
                      )

    # Get minimal set of gas phase candidates
    gas_phase_candidates = list(
        sorted(
            set(gas_phase_candidates)
        )
    )
    if options.verbose:
        print("\n\nGas phase candidates: {gas_phase_candidates}".format(
            gas_phase_candidates=gas_phase_candidates,
        ))
    if len(gas_phase_candidates) == 0:
        raise UserWarning(
            "No gas phase references found."
            "\nInclude additional folders with 'cathub organize -d foldername/'\n")

    volume_groups = {}
    tolerance = 1e-5
    if options.verbose:
        print("\nGroup by volume")
        print("==================")
    for surface in sorted(surfaces,
                          key=lambda x: x.get_volume(),):
        formula = symbols(surface)

        for volume in volume_groups:
            if abs(volume - surface.get_volume()) < tolerance:
                volume_groups[volume].append(surface)
                break
        else:
            volume_groups[surface.get_volume()] = [surface]

    for volume in volume_groups:
        if options.verbose:
            print("\nInspecting volume={volume}".format(
                volume=volume.round(2),
            ))
            print('====================')
        surfaces = volume_groups[volume]
        N = len(surfaces)
        # Order surfaces by number of atoms and energy
        surface_size = [len(s) for s in surfaces]
        sorted_surfaces = []
        for ss in sorted(set(surface_size)):
            idx = [i for i, s in enumerate(surface_size) if s == ss]
            subset = [surfaces[i] for i in idx]
            energies = [s.get_potential_energy() for s in subset]
            formulas = [s.get_chemical_formula() for s in subset]
            idx = np.argsort(energies)
            energies = np.sort(energies)
            subset = [subset[i] for i in idx]
            formulas = [formulas[i] for i in idx]
            if options.keep_all_energies or options.keep_all_slabs:
                subset = [s for i, s in enumerate(subset) if not
                          energies[i] in energies[:i]]
            else:
                subset = [s for i, s in enumerate(subset) if not
                          formulas[i] in formulas[:i]]
            sorted_surfaces += subset

        surfaces = sorted_surfaces
        if options.keep_all_slabs:
            n_empty = len(surface)
        else:
            n_empty = 1  # Only consider lowest energy empty slab for now
        for i, surf_empty in enumerate(surfaces[:n_empty]):
            for j, surf_ads in enumerate(surfaces[i+1:]):
                if surf_empty.get_chemical_formula() == surf_ads.get_chemical_formula():
                    continue
                if options.verbose:
                    print('\n    {} vs {}'.format(get_chemical_formula(surf_empty),
                                                  get_chemical_formula(surf_ads)))
                    print('    -------------------')

                # Check for calculator parameter consistency
                if check_parameters:
                    param_check = compare_parameters(surf_empty,
                                                     surf_ads)
                    if param_check == 2 and options.verbose:
                        print("        -Warning: No calculator information detected for"
                              " {} vs {}".format(surf_empty.info['filename'],
                                                 surf_ads.info['filename']))

                    elif not param_check:
                        if options.verbose:
                            print("\        n-Warning: Not included."
                                  " different calculator parameters detected for"
                                  " {} vs {}".format(surf_empty.info['filename'],
                                                     surf_ads.info['filename']))
                        continue
                if not options.skip_constraints:
                    constraints_empty = [c.todict()['kwargs']['indices']
                                         for c in surf_empty.constraints
                                         if c.todict()['name'] == 'FixAtoms']
                    constraints_ads = [c.todict()['kwargs']['indices']
                                       for c in surf_ads.constraints
                                       if c.todict()['name'] == 'FixAtoms']

                    c_flag = 0
                    if len(constraints_empty) > 0:
                        c_indices = constraints_empty[0]
                        constrained_positions = surf_empty[c_indices].positions
                        ads_positions = surf_ads.positions
                        for cp in constrained_positions:
                            if not np.any(np.all(np.isclose(cp, ads_positions,
                                                            atol=1e-4), axis=1)):
                                c_flag = 1
                                break
                        c_indices = constraints_ads[0]
                        constrained_positions = surf_ads[c_indices].positions
                        empty_positions = surf_empty.positions
                        for cp in constrained_positions:
                            if c_flag:
                                break
                            if not np.any(np.all(np.isclose(cp, empty_positions,
                                                            atol=1e-4), axis=1)):
                                c_flag = 1
                                break

                    if c_flag:
                        if options.verbose:
                            print("\        n-Warning: Not included."
                                  " different constraint settings detected for"
                                  " {} vs {}".format(surf_empty.info['filename'],
                                                     surf_ads.info['filename']))
                        continue

                atomic_num_ads = sorted(surf_ads.numbers)
                atomic_num_surf = sorted(surf_empty.numbers)

                diff_numbers = atomic_num_ads.copy()
                equal_numbers = []
                for n in atomic_num_surf:
                    if n in diff_numbers:
                        diff_numbers.remove(n)
                        equal_numbers += [n]
                if not sorted(equal_numbers) ==\
                   sorted(atomic_num_surf):
                    continue
                equal_formula = get_reduced_chemical_formula(
                    ase.atoms.Atoms(equal_numbers))

                adsorbate = get_reduced_chemical_formula(
                    ase.atoms.Atoms(diff_numbers))

                red_diff_numbers, rep = \
                    cathub.ase_tools.get_reduced_numbers(diff_numbers)

                if not red_diff_numbers in adsorbate_numbers:
                    #index = adsorbate_numbers.index(red_diff_numbers)
                    if options.verbose:
                        print("        -Adsorbate {} detected.".format(adsorbate),
                              "Include with 'cathub organize -a {}'".format(adsorbate))
                    continue

                dE = surf_ads.get_potential_energy() \
                    - surf_empty.get_potential_energy()

                dE /= rep

                references, prefactors = \
                    gas_phase_references \
                    .construct_reference_system(adsorbate,
                                                gas_phase_candidates)

                if not references:
                    print(
                        "        -Warning: Gas phase references could not be constructed for adsorbate {}.".format(adsorbate))
                    continue

                equation = ''
                for i, ref in enumerate(references):
                    dE -= prefactors[i] * reference_energy[ref]
                    pf = prefactors[i]
                    if pf == 1.0:
                        pf = ''
                    elif str(pf).isdigit():
                        pf = int(pf)
                    equation += '{}{}gas_'.format(pf, ref)

                equation += 'star__'

                facet = options.facet_name if options.facet_name != 'facet' \
                    else surf_empty.info['facet']

                site = surf_ads.info.get('site', None)

                if not abs(dE) < options.max_energy:
                    if options.verbose:
                        print("\n        -Adsorbate {} detected with adsorption energy: {} eV.".format(adsorbate, dE),
                              "This above current threshold of {} eV:".format(
                                  options.max_energy),
                              "Increase --max-energy to include")
                    continue
                surface_ads = ("{equal_formula}"
                               "({surface_facet})"
                               "+{adsorbate}").format(
                                   equal_formula=equal_formula,
                                   surface_facet=surf_empty.info['facet'],
                                   adsorbate=adsorbate)

                energy = dE
                key = equal_formula

                if not options.keep_all_energies:
                    if energy > np.min(list(collected_energies.get(
                            key, {}).get(facet, {}).get(adsorbate, {}).values()) +[np.inf]):
                        print('FOUND:', surface_ads)
                        continue
                else:
                    n_energies = len(collected_energies.get(
                        key, {}).get(facet, {}).get(adsorbate, {}).values())
                    if not site:
                        site = 'site{}'.format(n_energies + 1)

                if options.keep_all_slabs:
                    key += '_Epot={}'.format(surf_empty.get_potential_energy())

                if site:
                    equation += '{}@{}'.format(adsorbate, site)
                else:
                    equation += '{}star'.format(adsorbate)

                if options.verbose:
                    print("        -Adsorption energy found for {surface_ads}"
                          .format(**locals()))

                if options.interactive:
                    print(' ')
                    include = input('Include reaction: {}({}) | {} | dE={} ?\n  return(yes) / n(no) / u(update) '.
                                    format(key.split('_')[0], facet, equation.replace('__', '->'), round(dE, 3)))
                    if include == 'n':
                        continue
                    if include == ('u' or 'update'):
                        print('\nUpdating info for structure: "{}"'.format(
                            surf_ads.info['filename']))
                        print('  Provide updated info or press "return" to accept.')
                        for update_key in ['site', 'facet', 'adsorbate']:
                            update_value = input("    {} ({}): ".format(
                                update_key, locals().get(update_key)))
                            if not update_value:
                                continue
                            if update_key == 'site':
                                site = update_value
                                if '@' in equation:
                                    equation = equation.split('@')[0]
                                else:
                                    equation = 'star'.join(
                                        equation.split('star')[:-1])
                                equation += '@' + site

                            elif update_key == 'facet':
                                facet = update_value

                            elif update_key == 'adsorbate':
                                old_symbols = sorted(
                                    ase.atoms.Atoms(adsorbate).symbols)
                                new_symbols = sorted(
                                    ase.atoms.Atoms(update_value).symbols)
                                if not old_symbols == new_symbols:
                                    print(
                                        '    Error: new adsorbate must contain same atoms as {}'.format(adsorbate))
                                else:
                                    equation = equation.replace('{}star'.format(
                                        adsorbate), '{}star'.format(update_value))

                dft_code = options.dft_code or structure.info['filetype']
                dft_functional = options.xc_functional

                collected_structures \
                    .setdefault(dft_code, {}) \
                    .setdefault(dft_functional, {}) \
                    .setdefault(key, {}) \
                    .setdefault(facet, {}) \
                    .setdefault('empty_slab', surf_empty)

                key_count[key] = key_count.get(key, 0) + 1
                if rep > 1:
                    adsorbate = str(rep) + adsorbate
                if not key in collected_energies:
                    collected_energies[key] = {}
                if not facet in collected_energies[key]:
                    collected_energies[key][facet] = {}
                if not adsorbate in collected_energies[key][facet]:
                    collected_energies[key][facet][adsorbate] = {}
                if not energy in list(collected_energies[key][facet][adsorbate].values()):
                    collected_energies[key][facet][adsorbate][site] = energy
                    collected_structures[dft_code][dft_functional][key][facet]\
                        .setdefault(equation, {})[adsorbate] = surf_ads

    # if options.verbose:
    #    print("\n\nCollected Adsorption Energies Data")
    #    print("====================================")
    #    pprint.pprint(collected_structures, compact=True, depth=4)
    print("\n\nCollected Adsorption Energies")
    print("=============================")
    if len(collected_energies) == 0:
        print("Warning: no energies collected. Some ways to fix this:")
        print("  * raise the allowed maximum reaction energy (default: 10 eV)")
        print("     --max-energy 50")
        print("  * make sure you have gas phase molecules in the directory")
        print("  * raise the maximum density for gas-phase molecule")
        print("     --max-density-gas 0.004")
        print("  * raise the maximum density for slab structures")
        print("    --max-density-slab 0.03 ")

    for key, facets in collected_energies.items():
        for facet, adsorbates in facets.items():
            for ads, sites in adsorbates.items():
                for site, e in sites.items():            
                    print("{key:15s}: {energy:.3f} eV".format(
                        key='{}({}) + {}@{}'.format(key, facet, ads, site),
                        energy=e,
                    ))

    return collected_structures


def dict_representer(dumper, data):
    return dumper.represent_dict(data.items())


def create_folders(options, structures, publication_template, root=''):
    out_format = 'json'
    Dumper.add_representer(collections.OrderedDict, dict_representer)

    for key in structures:
        if isinstance(structures[key], dict):
            d = Path(root).joinpath(key)
            Path(d).mkdir(parents=True, exist_ok=True)
            if Path(root).parent.as_posix() == '.':
                # Have to explicitly convert Path to str
                # to work under python 3.4
                with open(str(Path(root).joinpath('publication.txt')),
                          'w') as outfile:
                    yaml.dump(publication_template,
                              outfile,
                              indent=4,
                              Dumper=Dumper)

                if options.energy_corrections:
                    with open(str(Path(root).joinpath('energy_corrections.txt')),
                              'w') as outfile:
                        yaml.dump(options.energy_corrections,
                                  outfile
                                  )

            create_folders(options, structures[key], publication_template={},
                           root=d)
        else:
            ase.io.write(
                str(Path(root).joinpath(key + '.' + out_format)),
                structures[key],
                format=out_format,
            )


def main(options):
    print("Running Organize script")
    print("=======================")
    pickle_file = options.foldername.strip().rstrip(
        '/').strip('.').rstrip('/') + '.cache.pckl'

    if Path(pickle_file).exists() \
            and Path(pickle_file).stat().st_size \
            and options.use_cache:
        with open(pickle_file, 'rb') as infile:
            structures = pickle.load(infile)
    else:
        structures = collect_structures(options.foldername,
                                        options.verbose,
                                        options.include_pattern,
                                        options.exclude_pattern,
                                        level='**/*')

        if options.gas_dir:
            structures.extend(
                collect_structures(
                    options.gas_dir,
                    options.verbose,
                    level='**/*')
            )
        if options.use_cache:
            with open(pickle_file, 'wb') as outfile:
                pickle.dump(structures, outfile)

    publication_template = cathub.ase_tools.PUBLICATION_TEMPLATE
    structures = fuzzy_match(structures, options)
    if not structures:
        return
    root = options.foldername.strip('/') + '.organized/'
    create_folders(options, structures,
                   root=root,
                   publication_template=publication_template,
                   )

    print("\nFile organization complete!"
          " Files placed in '{}'".format(root))

    print('\nInstructions:')
    print('=============')
    print("    1) Update the file '{root}publication.txt' with your publication info and email."
          "\n\n    2) Make sure DFT-CODE and XC-FUNCTIONAL folder names + FACET and @site extensions are changed"
          " with the right information."
          "\n\n    3) Run 'cathub folder2db {root}'".format(root=root),
          "to create a local database of reaction energies.\n")
