import subprocess
import os
import yaml
from yaml import Dumper
import click
import six
import collections
from tabulate import tabulate
from ase.symbols import string2symbols
from ase.cli import main
from . import query
from . import make_folders_template
from . import psql_server_connect
from . import folder2db as _folder2db
from . import db2server as _db2server
from . import organize as _organize
from . import folderreader
from . import ase_tools
from . import tools
from .cathubsqlite import CathubSQLite
from .postgresql import CathubPostgreSQL


@click.group()
@click.version_option()
def cli():
    pass


@cli.command()
@click.argument('dbfile')
def show_reactions(dbfile):
    """Extract and print reactions from sqlite3 (.db) file"""
    db = CathubSQLite(dbfile)
    db.print_summary()


@cli.command()
@click.argument('args',  default='', type=str)
@click.option('--dbuser', default='apiuser', type=str)
@click.option('--dbpassword', default='ubDwfqPw', type=str)
@click.option('--gui', default=False, show_default=True, is_flag=True,
              help='show structures in ase gui')
def ase(dbuser, dbpassword, args, gui):
    """Connection to atomic structures on the Catalysis-Hub
       server with ase db cli.
       Arguments to the the ase db cli client must be enclosed in one string.
       For example: <cathub ase 'formula=Ag6In6H -s energy -L 200'>.
       To see possible ase db arguments run <ase db --help>"""
    if dbuser == 'upload':
        dbpassword = 'cHyuuQH0'
    db = CathubPostgreSQL(user=dbuser, password=dbpassword)
    db._connect()
    server_name = db.server_name
    subprocess.call(
        ("ase db {} {}".format(server_name, args)).split())
    if gui:
        args = args.split('-')[0]
        subprocess.call(
            ('ase gui {}@{}'.format(server_name, args)).split())


@cli.command()
@click.argument('folder_name')
@click.option('--debug',
              is_flag=True,
              show_default=True,
              default=False)
@click.option(
    '--skip-folders',
    default='',
    show_default=True,
    help="""subfolders not to read, given as the name of a single folder,
    or a string with names of more folders seperated by ', '""")
@click.option(
    '--energy-limit',
    default=10.0,
    show_default=True,
    help="""Bounds for accepted absolute reaction energies in eV""")
@click.option('--goto-reaction',
              help="""name of reaction folder to skip ahead to""")
def folder2db(folder_name, debug, energy_limit, skip_folders,
              goto_reaction):
    """Read folder and collect data in local sqlite3 database"""

    folder_name = folder_name.rstrip('/')
    skip = []
    for s in skip_folders.split(', '):
        for sk in s.split(','):
            skip.append(sk)
    pub_id = _folder2db.main(folder_name, debug, energy_limit,
                             skip, goto_reaction)
    if pub_id:
        print('')
        print('')
        print('Ready to release the data?')
        print(
            "  Send it to the Catalysis-Hub server with 'cathub db2server {folder_name}/{pub_id}.db'.".format(**locals()))
        print("  Then log in at www.catalysis-hub.org/upload/ to verify and release. ")


@cli.command()
@click.argument('dbfile')
@click.option('--block-size', default=1000, type=int,
              help="Number of atomic structure to transfer per transaction",
              show_default=True)
@click.option('--dbuser', default='upload', type=str)
@click.option('--dbpassword', default='cHyuuQH0', type=str)
def db2server(dbfile, block_size, dbuser, dbpassword):
    """Transfer data from local database to Catalysis Hub server"""

    _db2server.main(dbfile,
                    write_reaction=True,
                    write_ase=True,
                    write_publication=True,
                    write_reaction_system=True,
                    block_size=block_size,
                    start_block=0,
                    user=dbuser,
                    password=dbpassword)


reaction_columns = [
    'chemicalComposition',
    'surfaceComposition',
    'facet',
    'sites',
    'coverages',
    'reactants',
    'products',
    'Equation',
    'reactionEnergy',
    'activationEnergy',
    'dftCode',
    'dftFunctional',
    'username',
    'pubId',
    'reactionSystems',
    'systems',
    'publication']
publication_columns = [
    'pubId',
    'title',
    'authors',
    'journal',
    'year',
    'doi',
    'tags']


@cli.command()
@click.option('--columns', '-c',
              default=('chemicalComposition', 'Equation', 'reactionEnergy'),
              type=click.Choice(reaction_columns),
              show_default=True,
              multiple=True)
@click.option('--n-results', '-n', default=10, show_default=True)
@click.option('--write-db', '-w', is_flag=True, default=False,
              show_default=True)
@click.option(
    '--queries',
    '-q',
    default={},
    multiple='True',
    show_default=True,
    help="""Make a selection on one of the columns:
    {0}\n Examples: \n -q chemicalComposition=~Pt for surfaces containing Pt
    \n -q reactants=CO for reactions with CO as a reactants"""
    .format(reaction_columns))
# Keep {0} in string.format for python2.6 compatibility
def reactions(columns, n_results, write_db, queries):
    """Search for reactions"""
    if not isinstance(queries, dict):
        query_dict = {}
        for q in queries:
            key, value = q.split('=')
            if key == 'distinct':
                if value in ['True', 'true']:
                    query_dict.update({key: True})
                    continue
            try:
                value = int(value)
                query_dict.update({key: value})
            except BaseException:
                query_dict.update({key: '{0}'.format(value)})
                # Keep {0} in string.format for python2.6 compatibility
    if write_db and n_results > 1000:
        print("""Warning: You're attempting to write more than a 1000 rows
        with geometries. This could take some time""")
    data = query.get_reactions(columns=columns,
                               n_results=n_results,
                               write_db=write_db,
                               **query_dict)

    if write_db:
        return
    table = []
    headers = []
    for row in data['reactions']['edges']:
        table += [list(row['node'].values())]

    headers = list(row['node'].keys())

    print(tabulate(table, headers) + '\n')


@cli.command()
@click.option('--columns', '-c',
              default=('pubId', 'title', 'authors', 'journal', 'year'),
              type=click.Choice(publication_columns),
              show_default=True,
              multiple=True)
@click.option('--n-results', '-n', default=10)
@click.option(
    '--queries',
    '-q',
    default={},
    multiple=True,
    show_default=True,
    help="""Make a selection on one of the columns:
    {0}\n Examples: \n -q: \n title=~Evolution \n authors=~bajdich
    \n year=2017""".format(publication_columns))
def publications(columns, n_results, queries):
    """Search for publications"""
    if not isinstance(queries, dict):
        query_dict = {}
        for q in queries:
            key, value = q.split('=')
            if key == 'distinct':
                if value in ['True', 'true']:
                    query_dict.update({key: True})
                    continue
            try:
                value = int(value)
                query_dict.update({key: value})
            except BaseException:
                query_dict.update({key: '{0}'.format(value)})
    if 'sort' not in query_dict:
        query_dict.update({'order': '-year'})
    data = query.query(table='publications',
                       columns=columns,
                       n_results=n_results,
                       queries=query_dict)
    table = []
    headers = []
    for row in data['publications']['edges']:
        value = list(row['node'].values())
        for n, v in enumerate(value):
            if isinstance(v, str) and len(v) > 20:
                splited = v.split(' ')
                size = 0
                sentence = ''
                for word in splited:
                    if size < 20:
                        size += len(word)
                        sentence += ' ' + word
                    else:
                        sentence += '\n' + word
                        size = 0
                sentence += '\n'
                value[n] = sentence

        table += [value]

    headers = list(row['node'].keys())
    print(tabulate(table, headers, tablefmt="grid") + '\n')


@cli.command()
@click.argument('template',
                default='template')
@click.option('--custom-base',
              help='Generate folders in a custom directory ')
def make_folders(template, custom_base):
    """Create a basic folder tree for dumping DFT calculcations for reaction energies.

    Dear all

    Use this command make the right structure for your folders
    for submitting data for Catalysis Hub's Surface Reactions.

    Start by creating a template file by calling:

    $ cathub make_folders  <template_name>

    Then open the template and modify it to so that it contains information
    about your data. You will need to enter publication/dataset information,
    and specify the types of surfaces, facets and reactions.

    The 'reactions' entry should include two lists for each reaction;
    'reactants' and 'products', corresponding to left- and right hand side of
    each chemical equation respectively.
    Remember to balance the equation by including a prefactor or minus sign
    in the name when relevant. For example:

    reactions:

    -    reactants: ['CCH3star@ontop']

         products: ['Cstar@hollow', 'CH3star@ontop']

    -    reactants: ['CH4gas', '-0.5H2gas', 'star']

         products:  ['CH3star']


    Please include the phase of the species as an extension:

        'gas' for gas phase (i.e. CH4 -> CH4gas)

        'star' for empty slab or adsorbed phase. (i.e. OH -> OHstar)

    The site of adsorbed species is also included as an extension:

        '@site' (i.e. OHstar in bridge-> OHstar@bridge)

    Energy corrections to gas phase molecules can be included as:

        energy_corrections: {H2: 0.1, CH4: -0.15}

    Then, save the template and call:

    $ cathub make_folders <template_name>

    And folders will be created automatically.

    You can create several templates and call make_folders again
    if you, for example, are using different functionals or are
    doing different reactions on different surfaces.

    After creating your folders, add your output files from the
    electronic structure calculations at the positions.
    Accepted file formats include everything that can be read by ASE
    and contains the total potential energy of the calculation, such
    as .traj or .OUTCAR files.

    After dumping your files, run `cathub folder2db <your folder>`
    to collect the data.
    """

    def dict_representer(dumper, data):
        return dumper.represent_dict(data.items())

    Dumper.add_representer(collections.OrderedDict, dict_representer)

    if custom_base is None:
        custom_base = os.path.abspath(os.path.curdir)
    template = custom_base + '/' + template

    template_data = ase_tools.REACTION_TEMPLATE
    if not os.path.exists(template):
        with open(template, 'w') as outfile:
            outfile.write(
                yaml.dump(
                    template_data,
                    indent=4,
                    Dumper=Dumper) +
                '\n')
            print("Created template file: {template}\n".format(**locals()) +
                  '  Please edit it and run the script again to create your folderstructure.\n' +
                  '  Run cathub make_folders --help for instructions')
            return

    with open(template) as infile:
        template_data = yaml.load(infile)
        title = template_data['title']
        authors = template_data['authors']
        journal = template_data['journal']
        volume = template_data['volume']
        number = template_data['number']
        pages = template_data['pages']
        year = template_data['year']
        email = template_data['email']
        publisher = template_data['publisher']
        doi = template_data['doi']
        dft_code = template_data['DFT_code']
        dft_functionals = template_data['DFT_functionals']
        reactions = template_data['reactions']
        crystal_structures = template_data['crystal_structures']
        bulk_compositions = template_data['bulk_compositions']
        facets = template_data['facets']
        energy_corrections = template_data['energy_corrections']

    make_folders_template.main(
        title=title,
        authors=eval(authors) if isinstance(
            authors, six.string_types) else authors,
        journal=journal,
        volume=volume,
        number=number,
        pages=pages,
        year=year,
        email=email,
        publisher=publisher,
        doi=doi,
        DFT_code=dft_code,
        DFT_functionals=dft_functionals,
        reactions=eval(reactions) if isinstance(
            reactions, six.string_types) else reactions,
        custom_base=custom_base,
        bulk_compositions=bulk_compositions,
        crystal_structures=crystal_structures,
        facets=facets,
        energy_corrections=energy_corrections
    )
    pub_id = tools.get_pub_id(title, authors, year)
    print(
        "Now dump your DFT output files into the folder, and run 'cathub folder2db {pub_id}'".format(**locals()))


@cli.command()
@click.argument('user', default='apiuser')
def connect(user):
    """Direct connection to PostreSQL server."""
    psql_server_connect.main(user)


@cli.command()
@click.argument(
    'foldername',
)
@click.option(
    '-a', '--adsorbates',
    type=str,
    default='C,O,N,H,S,OH,OOH,CH,CH2,CH3,CO,COH,NH,NH2,NH3,SH,SH2',
    show_default=True,
    help="Specify adsorbates that are to be included. (E.g. -a CO,O,H )")
@click.option(
    '-c', '--dft-code',
    default='DFT-CODE',
    type=str,
    show_default=True,
    help="Specify DFT Code used to calculate"
    " If not specified it will be generated from"
    " filetype the processed files.")
@click.option(
    '-e', '--exclude-pattern',
    type=str,
    default='',
    show_default=True,
    help="Regular expression that matches"
    " file (paths) are should be ignored.")
@click.option(
    '-f', '--facet-name',
    type=str,
    default='facet',
    show_default=True,
    help="Manually specify a facet names.")
@click.option(
    '-d', '--gas-dir',
    type=str,
    default='',
    show_default=True,
    help="Specify a folder where gas-phase molecules"
    " for calculating adsorption energies are located."
)
@click.option(
    '-g', '--max-density-gas',
    type=float,
    default=0.002,
    show_default=True,
    help="Specify the maximum density (#atoms/A^3)"
    " below which the structures are"
    " considered gas-phase molecules.")
@click.option(
    '-i', '--include-pattern',
    type=str,
    default='',
    show_default=True,
    help="Expressions that match"
         " only those files that are included.",)
@click.option(
    '-k', '--keep-all-energies',
    type=bool,
    is_flag=True,
    help="When multiple energies for the same facet and adsorbate"
    "are found keep all energies"
    "not only the most stable."
)
@click.option(
    '-m', '--max-energy',
    type=float,
    default=100.,
    show_default=True,
    help="Maximum absolute energy (in eV) that is considered.",)
@click.option(
    '-n', '--no-hydrogen',
    type=bool,
    is_flag=True,
    help="By default hydrogen is included as a gas-phase species"
         "to avoid using typically less accurate gas-phase references."
         "Use this flag to avoid using hydrogen."
)
@click.option(
    '-r', '--exclude-reference',
    type=str,
    default='',
    show_default=True,
    help="Gas phase reference molecules"
    " that should not be considered.")
@click.option(
    '-S', '--structure',
    default='',
    type=str,
    show_default=True,
    help='Bulk structure from which slabs where generated.'
    'E.g. fcc or A_a_225 for the general case.'
)
@click.option(
    '-s', '--max-density-slab',
    type=float,
    default=0.08,
    show_default=True,
    help="Specify the maximum density (#atoms/A^3) "
    " below which the structure are considered slabs and not bulk")
@click.option(
    '-t', '--traj-format',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help="Store intermediate filetype as traj"
    "instead of json files")
@click.option(
    '-u', '--use-cache',
    type=bool,
    is_flag=True,
    default=False,
    show_default=True,
    help="When set the script will cache"
    " structures between runs in a file named"
    " <FOLDER_NAME>.cache.pckl")
@click.option(
    '-v', '--verbose',
    is_flag=True,
    default=False,
    show_default=True,
    help="Show more debugging messages.")
@click.option(
    '-x', '--xc-functional',
    type=str,
    default='XC-FUNCTIONAL',
    show_default=True,
    help="Set the DFT exchange-correlation functional"
    " used to calculate total energies.")
@click.option(
    '--energy-corrections',
    default={},
    type=str,
    help="Energy correction to gas phase molecules.")
@click.option(
    '-p', '--skip-parameters',
    is_flag=True,
    default=False,
    show_default=True,
    help="Skip calculator parameter check.")
@click.option(
    '-sc', '--skip-constraints',
    is_flag=True,
    default=False,
    show_default=True,
    help="Skip constraint check.")


def organize(**kwargs):
    """Read reactions from non-organized folder"""

    # do argument wrangling  before turning it into an obect
    # since namedtuples are immutable
    if len(kwargs['adsorbates']) == 0:
        print("""Warning: no adsorbates specified,
        can't pick up reaction reaction energies.""")
        print("         Enter adsorbates like --adsorbates CO,O,CO2")
        print("         [Comma-separated list without spaces.]")
    kwargs['adsorbates'] = list(map(
        lambda x: (''.join(sorted(string2symbols(x)))),
        kwargs['adsorbates'].split(','),
    ))
    if kwargs['energy_corrections']:
        e_c_dict = {}
        for e_c in kwargs['energy_corrections'].split(','):
            key, value = e_c.split('=')
            e_c_dict.update({key: float(value)})
        kwargs['energy_corrections'] = e_c_dict

    options = collections.namedtuple(
        'options',
        kwargs.keys()
    )(**kwargs)
    _organize.main(options=options)
