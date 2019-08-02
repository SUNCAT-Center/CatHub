from .cathubsqlite import CathubSQLite
from .tools import get_bases, clear_prefactor, clear_state, get_pub_id, extract_atoms
from .ase_tools import collect_structures
from . import ase_tools

import sys
from datetime import date
import numpy as np
import os
import copy
import json
import yaml


class FolderReader:
    """
    Class for reading data from organized folders and writing to local
    CathubSQLite database. Folders should be arranged with
    make_folders_template and are read in the order:

    level:

    0    folder_name
    1    |-- publication
    2        |-- dft_code
    3            |-- dft_functional
    4                |-- gas
    4                |-- metal1
    5                    |-- facet
    6                        |-- reaction

    Parameters
    ----------
    foldername: str
    debug: bool
        default is False. Choose True if the folderreader should  continue
        in spite of errors.
    update: bool
        Update data if allready present in database file. defalt is True
    energy_limit: float
        Limit for acceptable absolute reaction energies
    """

    def __init__(self, folder_name, debug=False, strict=True, verbose=False,
                 update=True, energy_limit=5, stdin=sys.stdin,
                 stdout=sys.stdout):
        self.debug = debug
        self.strict = strict
        self.verbose = verbose
        self.update = update
        self.energy_limit = energy_limit

        self.data_base, self.user, self.user_base \
            = get_bases(folder_name=folder_name)
        self.user_base_level = len(self.user_base.split("/"))

        self.pub_level = 1
        self.DFT_level = 2
        self.XC_level = 3
        self.reference_level = 4
        self.slab_level = 5
        self.reaction_level = 6
        self.final_level = 6

        self.stdin = stdin
        self.stdout = stdout

        self.cathub_db = None
        self.coverages = None
        self.omit_folders = []
        self.doi = None
        self.title = None
        self.authors = None
        self.year = None
        self.tags = None
        self.pub_id = None
        self.warnings = []

    def read(self, skip=[], goto_metal=None, goto_reaction=None):
        """
        Get reactions from folders.

        Parameters
        ----------
        skip: list of str
            list of folders not to read
        goto_reaction: str
            Skip ahead to this metal
        goto_reaction:
            Skip ahead to this reacion
        """
        if len(skip) > 0:
            for skip_f in skip:
                self.omit_folders.append(skip_f)

        """ If publication level is input"""
        if os.path.isfile(self.data_base + '/publication.txt'):
            self.user_base_level -= 1

        self.stdout.write('---------------------- \n')
        self.stdout.write('Starting folderreader! \n')
        self.stdout.write('---------------------- \n')
        found_reaction = False
        for root, dirs, files in os.walk(self.user_base):
            for omit_folder in self.omit_folders:  # user specified omit_folder
                if omit_folder in dirs:
                    dirs.remove(omit_folder)
            level = len(root.split("/")) - self.user_base_level

            if level == self.pub_level:
                self.read_pub(root)

            if level == self.DFT_level:
                self.DFT_code = os.path.basename(root)

            if level == self.XC_level:
                self.DFT_functional = os.path.basename(root)
                self.gas_folder = root + '/gas/'
                self.read_gas()

            if level == self.reference_level:
                if 'gas' in os.path.basename(root):
                    continue

                if goto_metal is not None:
                    if os.path.basename(root) == goto_metal:
                        goto_metal = None
                    else:
                        dirs[:] = []  # don't read any sub_dirs
                        continue
                self.read_bulk(root)

            if level == self.slab_level:
                self.read_slab(root)

            if level == self.reaction_level:
                if goto_reaction is not None:
                    if os.path.basename(root) == goto_reaction:
                        goto_reaction = None
                    else:
                        dirs[:] = []  # don't read any sub_dirs
                        continue

                self.read_reaction(root)

            if level == self.final_level:
                self.root = root
                self.read_energies(root)
                if self.key_value_pairs_reaction is not None:
                    yield self.key_value_pairs_reaction

    def write(self, skip=[], goto_reaction=None):
        for key_values in self.read(skip=skip, goto_reaction=goto_reaction):
            with CathubSQLite(self.cathub_db) as db:
                id = db.check(
                    key_values['chemical_composition'],
                    key_values['reaction_energy'])
                if id is None:
                    try:
                        id = db.write(key_values)
                        self.stdout.write(
                            '  Written to reaction db row id = {}\n'.format(id))
                    except BaseException as e:
                        self.raise_error(
                            'Writing to db: {}. {}'.format(e, self.root))

                elif self.update:
                    db.update(id, key_values)
                    self.stdout.write(
                        '  Updated reaction db row id = {}\n'.format(id))
                else:
                    self.stdout.write(
                        '  Already in reaction db with row id = {}\n'.format(id))
        assert self.cathub_db is not None, \
            'Wrong folder! No reactions found in {base}'\
            .format(base=self.user_base)
        self.print_warnings()
        self.get_summary()

    def get_summary(self):
        with CathubSQLite(self.cathub_db) as db:
            db.print_summary()

    def write_publication(self, pub_data):
        with CathubSQLite(self.cathub_db) as db:
            pid = db.check_publication(self.pub_id)
            if pid is None:
                pid = db.write_publication(pub_data)
                self.stdout.write(
                    'Written to publications db row id = {}\n'.format(pid))
        return pid

    def read_pub(self, root):
        pub_folder = os.path.basename(root)
        publication_keys = {}
        try:
            with open(root + '/publication.txt', 'r') as f:
                pub_data = yaml.load(f)
            if 'url' in pub_data.keys():
                del pub_data['url']
            self.title = pub_data['title']
            self.authors = pub_data['authors']
            self.year = pub_data['year']
            if 'doi' not in pub_data:
                pub_data.update({'doi': None})
                self.stdout.write('ERROR: No doi\n')
            else:
                self.doi = pub_data['doi']
            if 'tags' not in pub_data:
                pub_data.update({'tags': None})
                self.stdout.write('ERROR: No tags\n')

            for key, value in pub_data.items():
                if isinstance(value, list):
                    value = json.dumps(value)
                else:
                    try:
                        value = int(value)
                    except BaseException:
                        pass

        except Exception as e:
            self.stdout.write(
                'ERROR: insufficient publication info {e}\n'.format(
                    **locals()))
            pub_data = {'title': None,
                        'authors': None,
                        'journal': None,
                        'volume': None,
                        'number': None,
                        'pages': None,
                        'year': None,
                        'publisher': None,
                        'doi': None,
                        'tags': None
                        }

        try:
            with open(root + '/energy_corrections.txt', 'r') as f:
                self.energy_corrections = yaml.load(f)
                if self.energy_corrections:
                    self.stdout.write('----------------------------\n')
                    self.stdout.write('Applying energy corrections:\n')
                    for key, value in self.energy_corrections.items():
                        if value > 0:
                            sgn = '+'
                        else:
                            sgn = '-'
                            sys.stdout.write('  {key}: {sgn}{value}\n'
                                             .format(key=key, sgn=sgn,
                                                     value=value))
                    self.stdout.write('----------------------------\n')
                else:
                    self.energy_corrections = {}

        except BaseException:
            self.energy_corrections = {}

        if pub_data['title'] is None:
            self.title = os.path.basename(root)
            pub_data.update({'title': self.title})
        if pub_data['authors'] is None:
            self.authors = [self.user]
            pub_data.update({'authors': self.authors})
        if pub_data['year'] is None:
            self.year = date.today().year
            pub_data.update({'year': self.year})
        if pub_data['email']:
            self.user = pub_data['email']

        self.pub_id = get_pub_id(self.title, self.authors, self.year)
        self.cathub_db = '{}{}.db'.format(self.data_base, self.pub_id)
        self.stdout.write(
            'Writing to .db file {}:\n \n'.format(self.cathub_db))
        pub_data.update({'pub_id': self.pub_id})
        pid = self.write_publication(pub_data)

    def read_gas(self):
        gas_structures = collect_structures(self.gas_folder)
        self.ase_ids_gas = {}
        self.gas = {}

        for gas in gas_structures:
            gas = gas[-1]
            ase_id = None
            found = False

            chemical_composition = \
                ''.join(sorted(ase_tools.get_chemical_formula(
                    gas, mode='all')))
            chemical_composition_hill = ase_tools.get_chemical_formula(
                gas, mode='hill')
            energy = gas.get_potential_energy()
            key_value_pairs = {"name": chemical_composition_hill,
                               'state': 'gas',
                               'epot': energy}

            id, ase_id = ase_tools.check_in_ase(
                gas, self.cathub_db)

            if ase_id is None:
                ase_id = ase_tools.write_ase(gas, self.cathub_db,
                                             self.stdout,
                                             self.user,
                                             **key_value_pairs)
            elif self.update:
                ase_tools.update_ase(self.cathub_db, id,
                                     self.stdout, **key_value_pairs)

            self.ase_ids_gas.update({chemical_composition: ase_id})
            self.gas.update({chemical_composition: gas})

    def read_bulk(self, root):
        basename = os.path.basename(root)
        assert '_' in basename, \
            """Wrong folderstructure! Folder should be of format
            <metal>_<crystalstructure> but found {basename}""".format(
                basename=basename
            )
        self.metal, self.crystal = basename.split('_', 1)

        self.stdout.write(
            '------------------------------------------------------\n')
        self.stdout.write(
            '                    Surface:  {}\n'.format(self.metal))
        self.stdout.write(
            '------------------------------------------------------\n')

        self.ase_ids = {}

        bulk_structures = collect_structures(root)
        n_bulk = len(bulk_structures)
        if n_bulk == 0:
            return
        elif n_bulk > 1:
            self.raise_warning('More than one bulk structure submitted at {root}'
                               .format(root=root))
            return

        bulk = bulk_structures[0][-1]
        ase_id = None
        energy = ase_tools.get_energies([bulk])

        key_value_pairs = {"name": self.metal,
                           'state': 'bulk',
                           'epot': energy}

        id, ase_id = ase_tools.check_in_ase(
            bulk, self.cathub_db)
        if ase_id is None:
            ase_id = ase_tools.write_ase(bulk, self.cathub_db, self.stdout,
                                         self.user, **key_value_pairs)
        elif self.update:
            ase_tools.update_ase(self.cathub_db, id,
                                 self.stdout, **key_value_pairs)

        self.ase_ids.update({'bulk' + self.crystal: ase_id})

    def read_slab(self, root):
        self.facet = root.split('/')[-1].split('_')[0]
        self.ase_facet = 'x'.join(list(self.facet))

        empty_structures = collect_structures(root)
        n_empty = len(empty_structures)

        if n_empty == 0:
            self.raise_warning('No empty slab submitted at {root}'
                               .format(root=root))
            self.empty = None
            return
        elif n_empty > 1:
            self.raise_warning('More than one empty slab submitted at {root}'
                               .format(root=root))
            filename_collapse = ''.join([empty[-1].info['filename']
                                         for empty in empty_structures])
            if 'TS' not in filename_collapse:
                return

        self.empty = empty_structures[0][-1]

        ase_id = None
        energy = ase_tools.get_energies([self.empty])
        key_value_pairs = {"name": self.metal,
                           'state': 'star',
                           'epot': energy}

        key_value_pairs.update({'species': ''})

        id, ase_id = ase_tools.check_in_ase(
            self.empty, self.cathub_db)

        if ase_id is None:
            ase_id = ase_tools.write_ase(self.empty, self.cathub_db, self.stdout,
                                         self.user, **key_value_pairs)
        elif self.update:
            ase_tools.update_ase(self.cathub_db, id,
                                 self.stdout, **key_value_pairs)
        self.ase_ids.update({'star': ase_id})

    def read_reaction(self, root):
        folder_name = os.path.basename(root)

        self.reaction, self.sites = ase_tools.get_reaction_from_folder(
            folder_name)  # reaction dict

        self.stdout.write(
            '----------- REACTION:  {} --> {} --------------\n'
            .format('+'.join(self.reaction['reactants']),
                    '+'.join(self.reaction['products'])))

        self.get_reaction_atoms()
        """Create empty dictionaries"""
        r_empty = ['' for n in range(len(self.reaction_atoms['reactants']))]
        p_empty = ['' for n in range(len(self.reaction_atoms['products']))]
        self.structures = {'reactants': r_empty[:],
                           'products': p_empty[:]}

        key_value_pairs = {}

        """ Match reaction gas species with their atomic structure """
        for key, mollist in self.reaction_atoms.items():
            for i, molecule in enumerate(mollist):
                if self.states[key][i] == 'gas':
                    assert molecule in self.ase_ids_gas.keys(), \
                        """Molecule {molecule} is missing in folder {gas_folder}"""\
                        .format(molecule=clear_prefactor(self.reaction[key][i]),
                                gas_folder=self.gas_folder)
                    self.structures[key][i] = self.gas[molecule]
                    species = clear_prefactor(
                        self.reaction[key][i])
                    key_value_pairs.update(
                        {'species': clear_state(species)})
                    self.ase_ids.update({species: self.ase_ids_gas[molecule]})

        """ Add empty slab to structure dict"""
        for key, mollist in self.reaction_atoms.items():
            if '' in mollist:
                n = mollist.index('')
                self.structures[key][n] = self.empty

    def read_energies(self, root):
        self.key_value_pairs_reaction = None

        slab_structures = collect_structures(root)

        if len(slab_structures) == 0:
            self.raise_warning('No structure files in {root}: Skipping this folder'
                               .format(root=root))
            return

        # Remove old species from ase_ids
        all_reaction_species = list(self.reaction.values())[0] \
            + list(self.reaction.values())[1]
        all_reaction_species = [clear_prefactor(rs) for rs in
                                all_reaction_species]

        for ase_id in list(self.ase_ids.keys()):
            if ase_id == 'star' or 'bulk' in ase_id:
                continue
            if not ase_id in all_reaction_species:
                del self.ase_ids[ase_id]

        if 'TS' in self.structures:  # Delete old TS
            del self.structures['TS']
            del self.structures['TSempty']
            del self.prefactors['TS']

        for k in list(self.structures.keys()):
            if 'neb' in k:
                del self.structures[k]

        neb_indices = [i for i, slab in enumerate(slab_structures) if 'neb' in
                       slab[-1].info['filename']]
        neb_names = {}

        if len(neb_indices) == 1:
            index = neb_indices[0]
            slab = slab_structures[index]
            f = slab[-1].info['filename']
            del slab_structures[index]
            neb_indices = []
            for i, s in enumerate(slab):
                s.info['filename'] = f
                slab_structures.append(s)
                index = len(slab_structures) - 1
                neb_indices += [index]
                neb_names.update({str(i): 'neb' + str(i)})

        elif len(neb_indices) > 1:
            for i in neb_indices:
                f = slab_structures[i][-1].info['filename']
                neb_names.update({str(i): os.path.basename(f).split('.')[0]})

        for i, slab in enumerate(slab_structures):
            if isinstance(slab, list):
                slab_structures[i] = slab[-1]

        empty = self.empty

        reactant_entries = self.reaction['reactants'] + \
            self.reaction['products']

        if not empty:
            if 'star' in reactant_entries and len(neb_indices) == 0:
                message = 'Empty slab needed for reaction!'
                self.raise_error(message)
                return
            else:
                empty0 = slab_structures[0]
                chemical_composition = ase_tools.get_chemical_formula(empty0)
                self.raise_warning("Using '{}' as a reference instead of empty slab"
                                   .format(empty0.info['filename']))
                empty_atn = list(empty0.get_atomic_numbers())
        else:
            chemical_composition = ase_tools.get_chemical_formula(empty)
            empty_atn = list(empty.get_atomic_numbers())

        self.prefactor_scale = copy.deepcopy(self.prefactors)
        for key1, values in self.prefactor_scale.items():
            self.prefactor_scale[key1] = [1 for v in values]

        key_value_pairs = {}

        key_value_pairs.update({'name':
                                chemical_composition,
                                'facet': self.ase_facet,
                                'state': 'star'})

        """ Match adsorbate structures with reaction entries"""
        for i, slab in enumerate(slab_structures):
            f = slab.info['filename']
            atns = list(slab.get_atomic_numbers())
            if not (np.array(atns) > 8).any() and \
               (np.array(empty_atn) > 8).any():
                self.raise_warning("Only molecular species for structure: {}"
                                   .format(f))
                continue

            """Get supercell size relative to empty slab"""
            supercell_factor = 1
            if len(atns) > len(empty_atn) * 2:  # different supercells
                supercell_factor = len(atns) // len(empty_atn)

            """Atomic numbers of adsorbate"""
            ads_atn = []
            if 'star' in reactant_entries and len(neb_indices) == 0:
                ads_atn = copy.copy(atns)
                for atn in empty_atn * supercell_factor:
                    try:
                        ads_atn.remove(atn)
                    except ValueError as e:
                        self.raise_error(
                            'Empty slab: {} contains species not in: {}'
                            .format(empty.info['filename'], f))
                ads_atn = sorted(ads_atn)
                if ads_atn == []:
                    self.raise_warning("No adsorbates for structure: {}"
                                       .format(f))
                    continue

            ase_id = None
            id, ase_id = ase_tools.check_in_ase(slab, self.cathub_db)
            key_value_pairs.update({'epot': ase_tools.get_energies([slab])})

            if 'empty' in f and 'TS' in f:  # empty slab for transition state
                self.structures.update({'TSempty': [slab]})
                self.prefactors.update({'TSempty': [1]})
                self.prefactor_scale.update({'TSempty': [1]})
                key_value_pairs.update({'species': ''})
                if ase_id is None:
                    ase_id = ase_tools.write_ase(slab, self.cathub_db,
                                                 self.stdout,
                                                 self.user,
                                                 **key_value_pairs)
                elif self.update:
                    ase_tools.update_ase(self.cathub_db, id, self.stdout,
                                         **key_value_pairs)
                self.ase_ids.update({'TSemptystar': ase_id})
                continue

            elif 'TS' in f:  # transition state
                self.structures.update({'TS': [slab]})
                self.prefactors.update({'TS': [1]})
                self.prefactor_scale.update({'TS': [1]})
                key_value_pairs.update({'species': 'TS'})
                if ase_id is None:
                    ase_id = ase_tools.write_ase(slab, self.cathub_db,
                                                 self.stdout,
                                                 self.user,
                                                 **key_value_pairs)
                elif self.update:
                    ase_tools.update_ase(self.cathub_db, id, self.stdout,
                                         **key_value_pairs)
                self.ase_ids.update({'TSstar': ase_id})
                continue

            if i in neb_indices:
                self.structures.update({neb_names[str(i)]: [slab]})
                key_value_pairs.update({'species': 'neb'})
                if ase_id is None:
                    ase_id = ase_tools.write_ase(slab,
                                                 self.cathub_db,
                                                 self.stdout,
                                                 self.user,
                                                 **key_value_pairs)
                elif self.update:
                    ase_tools.update_ase(self.cathub_db, id, self.stdout,
                                         **key_value_pairs)
                self.ase_ids.update({neb_names[str(i)]: ase_id})
                continue

            found = False
            for key, mollist in self.reaction.items():
                if found:
                    break
                for n, molecule in enumerate(mollist):
                    if found:
                        break
                    if not self.states[key][n] == 'star':  # only slab stuff
                        continue
                    if not self.structures[key][n] == '':  # allready found
                        continue
                    molecule = clear_state(clear_prefactor(molecule))
                    if molecule == '':
                        continue
                    molecule_atn = sorted(
                        ase_tools.get_numbers_from_formula(molecule))
                    for n_ads in range(1, 5):
                        mol_atn = sorted(molecule_atn * n_ads)
                        if (ads_atn == mol_atn or len(ads_atn) == 0):
                            found = True
                            reaction_side = key
                            mol_index = n
                            break

            if found:
                self.structures[reaction_side][mol_index] = slab
                species = clear_prefactor(
                    self.reaction[reaction_side][mol_index])
                id, ase_id = ase_tools.check_in_ase(slab, self.cathub_db)
                key_value_pairs.update(
                    {'species': clear_state(species),
                     'n': n_ads,
                     'site': str(self.sites.get(species, ''))})
                if ase_id is None:
                    ase_id = ase_tools.write_ase(slab,
                                                 self.cathub_db,
                                                 self.stdout,
                                                 self.user,
                                                 **key_value_pairs)
                elif self.update:
                    ase_tools.update_ase(self.cathub_db,
                                         id,
                                         self.stdout,
                                         **key_value_pairs)
                self.ase_ids.update({species: ase_id})

            # For high coverage, re-balance chemical equation
            if n_ads > 1 and not \
               self.prefactors[reaction_side][mol_index] == n_ads:
                for key1, states in self.states.items():
                    indices_gas = [i for i, s in enumerate(states) if
                                   s == 'gas']
                    for i in indices_gas:
                        self.prefactor_scale[key1][i] = n_ads

                    indices_ads = [i for i, s in enumerate(states) if
                                   s == 'star' and not
                                   self.reaction_atoms[key1][i] == '']
                    for i in indices_ads:
                        if key1 == reaction_side and i == mol_index:
                            continue
                        self.prefactor_scale[key1][i] = n_ads

                    if len(indices_ads) > 0:
                        n = len(indices_ads)
                        self.add_empty_slabs(key1, - (n_ads - 1) * n)
            # Assume adsorbate is balanced in equation, and balance empty slabs
            elif n_ads > 1 and self.prefactors[reaction_side][mol_index] == n_ads:
                self.prefactor_scale[reaction_side][mol_index] = 1/n_ads
                self.add_empty_slabs(reaction_side, n_ads - 1)

            if supercell_factor > 1:
                for key1, states in self.states.items():
                    indices = [i for i, r in enumerate(self.reaction) if
                               r == 'star']
                    for i in indices:
                        self.prefactor_scale[key1][mol_i] *= supercell_factor

        # Check that all structures have been found
        self.clear_extra_empty_slabs()
        structurenames = [s for s in list(self.structures.keys())
                          if s not in ['reactants', 'products']]
        for k in ['reactants', 'products']:
            structurenames += [s for s in self.structures[k] if s != ''
                               and s is not None]
        only_neb = np.any(['neb' in s for s in structurenames])
        surface_composition = self.metal

        if only_neb:
            if not self.empty:
                for ads in self.reaction_atoms['reactants']:
                    ads_atn = ase_tools.get_numbers_from_formula(ads)
                    for atn in ads_atn:
                        empty_atn.remove(atn)
                chemical_composition = \
                    ase_tools.get_formula_from_numbers(empty_atn, mode='metal')
            neb_numbers = []
            neb_energies = []
            for key, structure in self.structures.items():
                if key in ['reactants', 'products']:
                    continue

                neb_no = int(key.split('.')[0].replace('neb', ''))
                neb_numbers += [neb_no]
                neb_energies += [structure[0].get_potential_energy()]

            initial = neb_energies[np.argmin(neb_numbers)]
            final = neb_energies[np.argmax(neb_numbers)]
            TS = np.max(neb_energies)

            reaction_energy = final - initial
            activation_energy = TS - initial

            if activation_energy == 0 or activation_energy == reaction_energy:
                activation_energy = None
        else:
            for key, structurelist in self.structures.items():
                if '' in structurelist:
                    index = structurelist.index('')
                    molecule = clear_state(
                        clear_prefactor(self.reaction[key][index]))
                    if self.states[key][index] == 'star':
                        message = "Adsorbate '{}' not found for any structure files in '{}'."\
                            .format(molecule, root) + \
                            " Please check your adsorbate structures and the empty slab."
                    if self.states[key][index] == 'gas':
                        message = "Gas phase molecule '{}' not found for any structure files in '{}'."\
                            .format(molecule, self.gas_folder) + \
                            " Please check your gas phase references."
                    self.raise_error(message)
                    return

            original_prefactors = copy.deepcopy(self.prefactors)
            for key in self.prefactors:
                for i, v in enumerate(self.prefactors[key]):
                    self.prefactors[key][i] = original_prefactors[key][i] * \
                        self.prefactor_scale[key][i]

            reaction_energy = None
            activation_energy = None
            try:
                reaction_energy, activation_energy = \
                    self.get_reaction_energy()
            except BaseException as e:
                message = "reaction energy failed for files in '{}'"\
                    .format(root)
                self.raise_error(message + '\n' + str(e))

        if not -self.energy_limit < reaction_energy < self.energy_limit:
            self.raise_error('reaction energy is very large ({} eV)'
                             .format(reaction_energy) +
                             'for folder: {}. \n  '.format(root) +
                             'If the value is correct, you can reset the limit with cathub folder2db --energy-limit <value>. Default is --energy-limit=5 (eV)'
                             )
        if activation_energy is not None:
            if activation_energy < reaction_energy:
                self.raise_warning('activation energy is smaller than reaction energy: {} vs {} eV \n  Folder: {}'.format(
                    activation_energy, reaction_energy, root))
            if not activation_energy < self.energy_limit:
                self.raise_error(' Very large activation energy: {} eV \n  Folder: {}'
                                 .format(activation_energy, root))

        reaction_info = {'reactants': {},
                         'products': {}}

        for key in ['reactants', 'products']:
            for i, r in enumerate(self.reaction[key]):
                r = clear_prefactor(r)
                reaction_info[key].update({r: original_prefactors[key][i]})
        self.key_value_pairs_reaction = {
            'chemical_composition': chemical_composition,
            'surface_composition': surface_composition,
            'facet': self.facet,
            'sites': self.sites,
            'coverages': self.coverages,
            'reactants': reaction_info['reactants'],
            'products': reaction_info['products'],
            'reaction_energy': float(reaction_energy),
            'activation_energy': activation_energy,
            'dft_code': self.DFT_code,
            'dft_functional': self.DFT_functional,
            'pub_id': self.pub_id,
            'doi': self.doi,
            'year': int(self.year),
            'ase_ids': self.ase_ids,
            'energy_corrections': self.energy_corrections,
            'username': self.user}

    def raise_error(self, message):
        if self.debug:
            self.stdout.write('--------------------------------------\n')
            self.stdout.write('Error: ' + message + '\n')
            self.stdout.write('--------------------------------------\n')
            self.warnings.append('Error: ' + message)
        else:
            self.print_warnings()
            raise RuntimeError(message)

    def raise_warning(self, message):
        self.stdout.write('Warning: ' + message + '\n')
        self.warnings.append('Warning: ' + message)

    def print_warnings(self):
        self.stdout.write('-------------------------------------------\n')
        self.stdout.write('All errors and warnings: ' + '\n')
        for warning in self.warnings:
            self.stdout.write('    ' + warning + '\n')
        self.stdout.write('-------------------------------------------\n')

    def get_reaction_atoms(self):
        self.reaction_atoms = {'reactants': [],
                               'products': []}

        self.prefactors = {'reactants': [],
                           'products': []}

        self.states = {'reactants': [],
                       'products': []}

        for key, mollist in self.reaction.items():
            for molecule in mollist:
                atoms, prefactor = ase_tools.get_all_atoms(molecule)
                self.reaction_atoms[key].append(atoms)
                self.prefactors[key].append(prefactor)
                state = ase_tools.get_state(molecule)
                self.states[key].append(state)

        self.prefactors_TS = copy.deepcopy(self.prefactors)

        self.balance_slabs()

    def balance_slabs(self):
        """Balance the number of slabs on each side of reaction"""
        n_r, n_p = self.get_n_slabs()

        diff = n_p - n_r

        if abs(diff) > 0:
            if diff > 0:  # add empty slabs to left-hand side
                n_r += diff
                side = 'reactants'
            elif diff < 0:  # add to right-hand side
                diff *= -1  # diff should be positive
                n_p += diff
                side = 'products'

            if '' not in self.reaction_atoms[side]:
                self.append_reaction_entry(side, prefactor=diff,
                                           prefactor_TS=1)
            else:
                index = self.reaction_atoms[side].index('')
                self.prefactors[side][index] += diff
                if side == 'reactants':
                    self.prefactors_TS[side][index] += diff

        if n_r > 1:  # Balance slabs for transition state
            count_empty = 0
            if '' in self.reaction_atoms['reactants']:
                index = self.reaction_atoms['reactants'].index('')
                count_empty = self.prefactors_TS['reactants'][index]
                self.prefactors_TS['reactants'][index] = - \
                    (n_r - count_empty - 1)
            else:
                self.append_reaction_entry('reactants', prefactor=0,
                                           prefactor_TS=-n_r + 1)
        else:
            if '' in self.reaction_atoms['reactants']:
                index = self.reaction_atoms['reactants'].index('')
                self.prefactors_TS['reactants'][index] = 1

    def get_n_slabs(self):
        n_star = {'reactants': 0,
                  'products': 0}

        for key, statelist in self.states.items():
            indices = [i for i, s in enumerate(statelist) if s == 'star']
            for i in indices:
                n_star[key] += self.prefactors[key][i]

        return n_star['reactants'], n_star['products']

    def get_n_empty_slabs(self):
        n_star = {'reactants': 0,
                  'products': 0}

        for key, species in self.reaction_atoms.items():
            indices = [i for i, s in enumerate(species) if s == '']
            for i in indices:
                n_star[key] += self.prefactors[key][i]

        return n_star['reactants'], n_star['products']

    def append_reaction_entry(self, reaction_side, prefactor,
                              prefactor_TS=1, entry='star',
                              atoms='', state='star'):

        self.reaction[reaction_side].append(entry)
        self.reaction_atoms[reaction_side].append(atoms)
        self.prefactors[reaction_side].append(prefactor)
        if reaction_side == 'reactants':
            self.prefactors_TS[reaction_side].append(prefactor_TS)
        self.states[reaction_side].append(state)

    def delete_reaction_entry(self, reaction_side, index):

        assert len(self.reaction[reaction_side]) == \
            len(self.reaction_atoms[reaction_side])
        del self.reaction[reaction_side][index]
        del self.reaction_atoms[reaction_side][index]
        del self.states[reaction_side][index]
        del self.prefactors[reaction_side][index]
        del self.prefactors_TS[reaction_side][index]
        del self.structures[reaction_side][index]

    def add_empty_slabs(self, reaction_side, n_slabs):
        found_slab = False
        for key, atoms in self.reaction_atoms.items():
            if '' in atoms:
                found_slab = True
                i = atoms.index('')
                if key == reaction_side:
                    self.prefactors[key][i] += n_slabs
                else:  # substract from other side
                    self.prefactors[key][i] -= n_slabs
        if not found_slab:
            self.append_reaction_entry(reaction_side,
                                       n_ads - 1)

    def clear_extra_empty_slabs(self):
        n_r, n_p = self.get_n_empty_slabs()
        if n_r > 0 and n_p > 0:
            diff = min([n_r, n_p])
            for side, species in self.reaction_atoms.items():
                i = species.index('')
                self.prefactors[side][i] -= diff
                if self.prefactors[side][i] == 0:
                    self.delete_reaction_entry(side, i)

    def get_reaction_energy(self):
        energies = {}
        for key in self.structures.keys():
            energies.update(
                {key: ['' for n in range(len(self.structures[key]))]})
        for reaction_side, atoms_list in self.structures.items():
            for i, atoms in enumerate(atoms_list):
                Ecor = self.get_energy_correction(reaction_side, i)
                energies[reaction_side][i] = self.prefactors[reaction_side][i] * \
                    (atoms.get_potential_energy() + Ecor)

        # Reaction energy:
        energy_reactants = np.sum(energies['reactants'])
        energy_products = np.sum(energies['products'])

        reaction_energy = energy_products - energy_reactants

        # Activation energy
        if 'TS' in self.structures.keys():
            # Is a different empty surface used for the TS?
            if 'TSempty' in self.structures.keys():
                for key in reaction_atoms.keys():
                    if '' in reaction_atoms[key]:
                        index = self.reaction_atoms[key].index('')
                        empty = self.structures[key][index]
                tsempty = self.structures['TSempty'][0]
                tsemptydiff = tsempty.get_potential_energy - \
                    empty.get_potential_energy()

            for i, structure in enumerate(self.structures['reactants']):
                Ecor = self.get_energy_correction('reactants', i)
                energies['reactants'][i] = self.prefactors_TS['reactants'][i]\
                    * structure.get_potential_energy() + Ecor
                if 'TSempty' in self.structures.keys() and \
                   self.states['reactants'][i] == 'star':
                    energies['reactants'][i] += self.prefactors_TS['reactants'][i]\
                        * tsemptydiff
            energy_reactants = np.sum(energies['reactants'])
            energy_TS = energies['TS'][0]
            activation_energy = energy_TS - energy_reactants
        else:
            activation_energy = None

        return reaction_energy, activation_energy

    def get_energy_correction(self, reaction_side, index):
        try:
            name = clear_prefactor(self.reaction[reaction_side][index])
        except BaseException:
            name = None
        Ecor = self.energy_corrections.get(name, 0)

        return Ecor
