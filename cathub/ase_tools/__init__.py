import sys
import collections
from functools import reduce
from fractions import gcd
from ase import Atoms
from ase.io import read
import numpy as np
import ase
from ase.utils import formula_metal
import copy
from cathub.tools import clear_state, get_state, clear_prefactor, get_prefactor

from pathlib import Path
Path().expanduser()


PUBLICATION_TEMPLATE = collections.OrderedDict({
    'title': 'Fancy title',
    'authors': ['Doe, John', 'Einstein, Albert'],
    'journal': 'JACS',
    'volume': '1',
    'number': '1',
    'pages': '23-42',
    'year': '2017',
    'email': 'winther@stanford.edu',
    'publisher': 'ACS',
    'doi': '10.NNNN/....',
})

REACTION_TEMPLATE = collections.OrderedDict({
    'title': 'Fancy title',
    'authors': ['Doe, John', 'Einstein, Albert'],
    'journal': 'JACS',
    'volume': '1',
    'number': '1',
    'pages': '23-42',
    'year': '2017',
    'email': 'winther@stanford.edu',
    'publisher': 'ACS',
    'doi': '10.NNNN/....',
    'DFT_code': 'Quantum Espresso',
    'DFT_functionals': ['BEEF-vdW', 'HSE06'],
    'reactions': [
        collections.OrderedDict({'reactants':
                                 ['2.0H2Ogas', '-1.5H2gas', 'star'],
                                 'products': ['OOHstar@top']}),
        collections.OrderedDict({'reactants': ['CCH3star@bridge'],
                                 'products':
                                 ['Cstar@hollow', 'CH3star@ontop']}),
        collections.OrderedDict({'reactants':
                                 ['CH4gas', '-0.5H2gas', 'star'],
                                 'products': ['CH3star@ontop']})
    ],
    'bulk_compositions': ['Pt'],
    'crystal_structures': ['fcc', 'hcp'],
    'facets': ['111'],
    'energy_corrections': {},
})


def get_chemical_formula(atoms, mode='metal'):
    """
    Compatibility function, return mode=metal, when
    available, mode=hill, when not (ASE <= 3.13)
    """
    try:
        return atoms.get_chemical_formula(mode=mode)
    except ValueError:
        return atoms.get_chemical_formula(mode='hill')


def get_reduced_chemical_formula(atoms):
    numbers = atoms.numbers
    unique_numbers, counts = np.unique(numbers, return_counts=True)
    denominator = reduce(gcd, counts)
    reduced_numbers = []
    for i, atomic_number in enumerate(unique_numbers):
        reduced_count = int(counts[i] / denominator)
        reduced_numbers += [atomic_number] * reduced_count
    return formula_metal(reduced_numbers)


def symbols(atoms):
    formula = get_chemical_formula(atoms)
    symbols = ase.symbols.string2symbols(formula)
    return ''.join(symbols)


def collect_structures(foldername, verbose=False, level='*'):
    structures = []
    if verbose:
        print(foldername)
    for i, filename in enumerate(Path(foldername).glob(level)):
        posix_filename = str(filename.as_posix())
        if verbose:
            print(i, posix_filename)
        if posix_filename.endswith('publication.txt'):
            with open(posix_filename) as infile:
                global PUBLICATION_TEMPLATE
                PUBLICATION_TEMPLATE = infile.read()
            continue
        if posix_filename.endswith('traj.old'):
            continue
        elif Path(posix_filename).is_file():
            try:
                filetype = ase.io.formats.filetype(posix_filename)
            except Exception as e:
                continue
            if filetype:
                if filetype == 'db':
                    with ase.db.connect(posix_filename) as db:
                        count = db.count()
                        print('Processing ASE db with {} structures'.format(count))
                        for row in db.select('energy'):
                            structure = [row.toatoms()]
                            structure[-1].info['filename'] = posix_filename
                            structure[-1].info['filetype'] = ase.io.formats.filetype(
                                posix_filename)
                            structures += [structure]
                else:
                    try:
                        structure = ase.io.read(posix_filename, ':')
                        structure[-1].info['filename'] = posix_filename
                        structure[-1].info['filetype'] = ase.io.formats.filetype(
                            posix_filename)
                        try:
                            structure[-1].get_potential_energy()
                            # ensure that the structure has an energy
                            structures.append(structure)
                        except RuntimeError:
                            print("Did not add {posix_filename} since it has no energy"
                                  .format(
                                      posix_filename=posix_filename,
                                  ))
                    except TypeError:
                        print("Warning: Could not read {posix_filename}"
                              .format(
                                  posix_filename=posix_filename,
                              ))

                    except StopIteration:
                        print("Warning: StopIteration {posix_filename} hit."
                              .format(
                                  posix_filename=posix_filename,
                              ))
                    except IndexError:
                        print("Warning: File {posix_filename} looks incomplete"
                              .format(
                                  posix_filename=posix_filename,
                              ))
                    except OSError as e:
                        print("Error with {posix_filename}: {e}".format(
                            posix_filename=posix_filename,
                            e=e,
                        ))
                    except AssertionError as e:
                        print("Hit an assertion error with {posix_filename}: {e}".format(
                            posix_filename=posix_filename,
                            e=e,
                        ))
                    except ValueError as e:
                        print("Trouble reading {posix_filename}: {e}".format(
                            posix_filename=posix_filename,
                            e=e,
                        ))
                    except DeprecationWarning as e:
                        print("Trouble reading {posix_filename}: {e}".format(
                            posix_filename=posix_filename,
                            e=e,
                        ))
                    except ImportError as e:
                        print("Trouble reading {posix_filename}: {e}".format(
                            posix_filename=posix_filename,
                            e=e,
                        ))
                    except ase.io.formats.UnknownFileTypeError as e:
                        print("Trouble reading {posix_filename}: {e}".format(
                            posix_filename=posix_filename,
                            e=e,
                        ))
    return structures


def get_energies(atoms_list):
    """ Potential energy for a list of atoms objects"""
    if len(atoms_list) == 1:
        return atoms_list[0].get_potential_energy()
    elif len(atoms_list) > 1:
        energies = []
        for atoms in atoms_list:
            energies.append(atoms.get_potential_energy())
        return energies


def get_atomic_numbers(atoms):
    return list(atoms.get_atomic_numbers())


def get_formula_from_numbers(numbers, mode='all'):
    formula = Atoms(numbers).get_chemical_formula(mode=mode)
    return formula


def get_numbers_from_formula(formula):
    atoms = Atoms(formula)
    return get_atomic_numbers(atoms)


def check_in_ase(atoms, ase_db, energy=None):
    """Check if entry is allready in ASE db"""

    db_ase = ase.db.connect(ase_db)
    if energy is None:
        energy = atoms.get_potential_energy()
    formula = get_chemical_formula(atoms)
    rows = db_ase.select(energy=energy)
    n = 0
    ids = []
    for row in rows:
        if formula == row.formula:
            return row.id, row.unique_id
    return None, None


def _normalize_key_value_pairs_inplace(data):
    for key in data:
        if isinstance(data[key], np.int64):
            data[key] = int(data[key])


def write_ase(atoms, db_file, stdout=sys.stdout, user=None, data=None,
              **key_value_pairs):
    """Connect to ASE db"""
    db_ase = ase.db.connect(db_file)
    _normalize_key_value_pairs_inplace(key_value_pairs)
    id = db_ase.write(atoms, data=data, **key_value_pairs)
    stdout.write('  writing atoms to ASE db row id = {}\n'.format(id))
    unique_id = db_ase.get(id)['unique_id']
    return unique_id


def update_ase(db_file, identity, stdout, **key_value_pairs):
    """Connect to ASE db"""
    db_ase = ase.db.connect(db_file)

    _normalize_key_value_pairs_inplace(key_value_pairs)
    count = db_ase.update(identity, **key_value_pairs)
    stdout.write('  Updating {0} key value pairs in ASE db row id = {1}\n'
                 .format(count, identity))
    return


def get_reaction_from_folder(folder_name):
    reaction = {}
    assert '__' in folder_name, "Please use __ as reaction arrow"

    reaction.update({'reactants': folder_name.split('__')[0].split('_'),
                     'products': folder_name.split('__')[1].split('_')})

    sites = {}
    for key, mollist in reaction.items():
        for n, mol in enumerate(mollist):
            if '@' in mol:
                mol, site = mol.split('@')
                clean_mol = clear_prefactor(mol)
                if not clean_mol in sites:
                    sites.update({clean_mol: site})
                else:
                    site0 = sites.get(clean_mol)
                    if not site0 == site:  # sites is a list
                        diff_sites = [site0, site]
                        sites.update({clean_mol: diff_sites})
                reaction[key][n] = mol
            if not np.any([state in mol for state in ['gas', 'aq', 'star']]):
                reaction[key][n] = mol + 'star'

    for key, mollist in reaction.items():
        n_star = mollist.count('star')
        if n_star > 1:
            for n in range(n_star):
                mollist.remove('star')
            mollist.append(str(n_star) + 'star')

    return reaction, sites


def get_all_atoms(molecule):
    molecule = clear_state(molecule)
    molecule, prefactor = get_prefactor(molecule)

    atoms = Atoms(molecule)

    molecule = ''.join(sorted(atoms.get_chemical_formula(mode='all')))

    return molecule, prefactor


def debug_assert(expression, message, debug=False):
    if debug:
        try:
            assert expression, message
        except AssertionError as e:
            print(e)
            return False
    else:
        assert expression, message

    return True
