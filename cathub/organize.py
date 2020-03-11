#!/usr/bin/env python

# builtin imports
import cathub.ase_tools
from .ase_tools import gas_phase_references, get_chemical_formula, \
    get_reduced_chemical_formula, symbols, collect_structures
import numpy as np
import ase.io
import ase.utils
import ase.atoms
import pickle
import json
import yaml
from yaml import Dumper
import re
import pprint
import collections

from pathlib import Path
Path().expanduser()

# other library imports


# local imports

np.set_printoptions(threshold=500, linewidth=1800, edgeitems=80)


def fuzzy_match(structures, options):
    # filter out cell with ill-defined unit cells
    structures = [structure[-1] for structure in structures
                  if structure[-1].number_of_lattice_vectors == 3
                  ]
    # sort by density
    structures = sorted(structures,
                        key=lambda x: len(x) / x.get_volume()
                        )
    # print(options.adsorbates)
    #print([a for a in options.adsorbates])
    adsorbate_numbers = [sorted(ase.atoms.Atoms(ads).numbers)
                         for ads in options.adsorbates]  # .split(',')]
    # group in to bulk, surface, or bulk
    molecules, surfaces, bulks = [], [], []
    gas_phase_candidates = []
    reference_energy = {}
    collected_energies = {}
    key_count = {}
    collected_structures = {}
    if options.verbose:
        print("Group By Densities")
        print("===================")
    for structure in structures:
        if options.include_pattern:
            if not re.search(
                    options.include_pattern, structure.info['filename']):
                continue
            if options.exclude_pattern \
                    and re.search(
                        options.exclude_pattern,
                        structure.info['filename']):
                continue

        elif options.exclude_pattern:
            if re.search(
                    options.exclude_pattern,
                    structure.info['filename']):
                continue

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
            print("  {density:10.3f} {filename}".format(
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
                    print("           GAS", formula,
                          structure.info['filename'])

        elif density < options.max_density_slab:
            structure.info['state'] = 'surface'
            formula = get_chemical_formula(structure)
            surfaces.append(structure)
            if options.verbose:
                print("           SURFACE", formula,
                      structure.info['filename'])
        else:
            structure.info['state'] = 'bulk'
            bulks.append(structure)

            if options.verbose:
                print("           BULK", formula, structure.info['filename'])

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
        print("\n\n  Gas phase candidates: {gas_phase_candidates}".format(
            gas_phase_candidates=gas_phase_candidates,
        ))

    volume_groups = {}
    tolerance = 1e-5
    if options.verbose:
        print("\nGroup by volume")
        print("===============")
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
        surfaces = volume_groups[volume]
        N = len(surfaces)
        # Order surfaces by number of atoms
        surface_size = [len(s) for s in surfaces]
        idx = np.argsort(surface_size)
        surfaces = [surfaces[i] for i in idx]

        for i, surf_empty in enumerate(surfaces):
            for j, surf_ads in enumerate(surfaces[i+1:]):
                if options.verbose:
                    print('\n    {} vs {}'.format(get_chemical_formula(surf_empty),
                                                  get_chemical_formula(surf_ads)))

                # Check for calculator parameter consistency
                if surf_empty.calc and surf_ads.calc:
                    if not surf_empty.calc.parameters == surf_ads.calc.parameters:
                        if options.verbose:
                            print("\nWarning: Not included."
                                  " different calculator parameters detected for"
                                  " {} vs {}".format(surf_empty.info['filename'],
                                                     surf_ads.info['filename']))
                        continue
                elif options.verbose:
                    print("\nWarning: No calculator information detected for"
                          " {} vs {}".format(surf_empty.info['filename'],
                                             surf_ads.info['filename']))
                    # " \nPlease include calculator information when possible")

                if not surf_ads.constraints == surf_empty.constraints:
                    # will lead to errors if adsorbate is constrained
                    if options.verbose:
                        print("\nWarning: Not included."
                              " different constraint settings detected for"
                              " {} vs {}".format(surf_empty.info['filename'],
                                                 surf_ads.info['filename']))
                    continue

                atomic_num_ads = sorted(surf_ads.numbers)
                atomic_num_surf = sorted(surf_empty.numbers)

                equal_numbers = [
                    n for n in atomic_num_ads if n in atomic_num_surf]
                diff_numbers = [n for n in atomic_num_ads
                                if not n in atomic_num_surf]

                # if not equal_numbers + diff_numbers == atomic_num_ads \
                #   or len(diff_numbers) == 0:
                #    continue

                equal_formula = get_reduced_chemical_formula(
                    ase.atoms.Atoms(equal_numbers))

                adsorbate = get_reduced_chemical_formula(
                    ase.atoms.Atoms(diff_numbers))

                red_diff_numbers, rep = \
                    cathub.ase_tools.get_reduced_numbers(diff_numbers)

                if not red_diff_numbers in adsorbate_numbers:
                    index = adsorbate_numbers.index(red_diff_numbers)
                    print("\nAdsorbate {} detected".format(options.adsorbates[index]),
                          " Include by setting the --adsorbates options.")

                dE = surf_ads.get_potential_energy() \
                    - surf_empty.get_potential_energy()

                dE /= rep

                references, prefactors = \
                    gas_phase_references \
                    .construct_reference_system(adsorbate,
                                                gas_phase_candidates)

                #stoich_factors = stoichiometry_factors[adsorbate]
                formula = ''

                for i, ref in enumerate(references):
                    dE -= prefactors[i] * reference_energy[ref]
                    pf = prefactors[i]
                    if pf == 1.0:
                        pf = ''
                    elif pf.isdigit():
                        pf = int(pf)
                    formula += '{}{}gas_'.format(pf, ref)

                formula += 'star__'

                site = surf_ads.info.get('site', None)
                if site:
                    formula += '{}@{}'.format(adsorbate, site)
                else:
                    formula += '{}star'.format(adsorbate)

                key = ("{equal_formula}"
                       "({surface_facet})"
                       "+{adsorbate}").format(
                           equal_formula=equal_formula,
                           surface_facet=surf_empty.info['facet'],
                           adsorbate=adsorbate)

                if abs(dE) < options.max_energy:
                    energy = dE
                    if options.verbose:
                        print("      Adsorption energy found for {key}"
                              .format(**locals()))

                    equation = (" {formula:30s}").format(
                        formula=formula,
                        equal_formula=equal_formula,
                        surface_facet=surf_empty.info['facet']) \
                        .replace(' ', '') \
                        .replace('+', '_') \
                        .replace('->', '__') \
                        .replace('*', 'star') \
                        .replace('(g)', 'gas')

                    # We keep the empty structure whether or not
                    # we keep all structures
                    collected_structures \
                        .setdefault(
                            options.dft_code
                            or structure.info['filetype'],
                            {}) \
                        .setdefault(options.xc_functional, {}) \
                        .setdefault(equal_formula + ('_' +
                                                     options.structure
                                                     or ''
                                                     ), {}) \
                        .setdefault(
                            options.facet_name
                            if options.facet_name != 'facet'
                            else surf_empty.info['facet'], {}) \
                        .setdefault('empty_slab', surf_empty)

                    collected_energies[key] = energy
                    key_count[key] = key_count.get(key, 0) + 1
                    # if options.verbose:
                    # print(key)
                    # print(collected_energies)
                    #print(key in collected_energies)
                    if not options.keep_all_energies:
                        if energy > collected_energies.get(
                                key, float("inf")):
                            continue

                    # persist adsorbate slab structures
                    ####################################
                    collected_energies[key] = energy
                    collected_structures.setdefault(
                        options.dft_code
                        or structure.info['filetype'],
                        {}).setdefault(
                        options.xc_functional,
                        {}).setdefault(
                        equal_formula + '_' + (
                            options.structure
                            or 'structure'), {}
                    ).setdefault(
                        options.facet_name
                        if options.facet_name != 'facet'
                        else surf_empty.info['facet'],
                        {}).setdefault(
                        equation,
                        {})[adsorbate] = surf_ads

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

    for key, energy in collected_energies.items():
        print("{key:40s}: {energy:.3f} eV".format(
            key=key,
            energy=energy,
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
    pickle_file = options.foldername.strip().rstrip(
        '/').strip('.').rstrip('/') + '.cache.pckl'

    if Path(pickle_file).exists() \
            and Path(pickle_file).stat().st_size \
            and options.use_cache:
        with open(pickle_file, 'rb') as infile:
            structures = pickle.load(infile)
    else:
        structures = collect_structures(options.foldername, options.verbose,
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
    create_folders(options, structures,
                   root=options.foldername.strip('/') + '.organized',
                   publication_template=publication_template,
                   )
