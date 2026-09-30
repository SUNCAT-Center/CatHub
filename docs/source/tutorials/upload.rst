Submitting data
---------------

Please submit your electronic structures calculations for surface reactions! Data submissions will be part of the Surface Reactions app at http://www.catalysis-hub.org/energies and is open to all institutions. Furthermore, the atomic structures that are part of your reaction can be utilized for other apps.

Your publication/dataset will be listed on the http://www.catalysis-hub.org/publications page, with a link to the publishers homepage (if already published). Your data will be easily accessible to other researchers, who will be able to browse through reaction energies and atomic structures from your publication.


Installing CatHub
...........................
CatHub is a python module that is used to interface the Surface Reactions database of Catalysis-Hub, directly from a python script of the command line. CatHub will be used to arrange your data into folders and submit your data to the server. To install CatHub, use pip::

  pip3 install git+https://github.com/SUNCAT-Center/CatHub.git --upgrade --user

which will install CatHub and all their dependencies.

To test that the cathub cli is working, start by typing::

  cathub

And you should see a list of subcommands. If it's not working you probably have to add the installation path to :code:`PATH` in your ~/.bashrc. This would typically be :code:`export PATH=~/.local/bin:${PATH}` for Linux, and :code:`export PATH~/Library/PythonX.Y/bin:${PATH}` for Mac.

Organizing data
....................

You have two options for organizing your data:

* cathub organize: For larger systematic datasets without reaction barriers, this approach will create folders and and arrange your data-files in the right location for you.

* cathub make-folders: For smaller or more complicated datasets with reaction barriers, this method will only create your folders, and you will have to drop the files in the right location yourself.

In either case no data will be uploaded to `catalysis-hub.org/publications <https://www.catalysis-hub.org/publications>`_ before the data is submitted to the server. **Note:** the ``cathub db2server`` command is only available to SLAC affiliates. If you are not affiliated with SLAC, please send your ``.db`` file by email to winther@slac.stanford.edu instead of running ``cathub db2server``.
Once you upload data it will be held in a moderation stage which you can inspect yourself at `catalysis-hub.org/upload <https://www.catalysis-hub.org/upload>`_
and delete yourself to iterate until all structures and energies look as expected. Once you are satisfied with your uploaded dataset there will be a "Release"
button that will notify the platform administrator that the dataset is ready for release.

cathub organize
................
This tool will take all your structure files from a general folder and organize them in the right folder-structure that can be used for data submission. Note: this approach does not work for transition states / barrier calculations. And it will still need a lot of manual file organization for co-adsorbate configurations.  While we are working on this cathub organize might still give you a nice head start with file organization.

To learn about the organize command, type::

  cathub organize --help

To read the data from a general folder, type::

  cathub organize <FOLDER> -a ADS1,ADS2 -c <dft-code> -x <xc-functional> -f <facet> -S <crystal structure>

Use the ``-a`` option to specify which adsorbates to look for. Also, please use the ``-c`` and ``-x`` options to specify the DFT code and xc-functional respectively. Furthermore, you are highly encurraged to use the ``-f`` and ``-S`` options to specify the surface facet and crystal structure when applicable.

This will generate an organized folder named ``<FOLDER>.organized``. Please open the .txt file ``<FOLDER>.organized/publication.txt``, and update it with info about your publication. It should look something like this::

    volume: 8
    publisher: Wiley
    doi: 10.1002/cssc.201500322
    title: "The Challenge of Electrochemical Ammonia Synthesis: A New Perspective on the Role of Nitrogen Scaling Relations"
    journal: ChemSusChem
    authors: [Montoya, Joseph H., Tsai, Charlie, Vojvodic, Aleksandra, Norskov, Jens K.]
    year: 2015
    email: winther@stanford.edu
    number: 13
    tags: []
    pages: 2140-2267


Remember your ``email`` since it will be used to log in at http://www.catalysis-hub.org/upload. Note that authors should be a list, with names in the form "lastname, firstname M.".

Please go through the created folder and rename folders to make your data easier to localize later. For example, a structure folder like Pt16_structure, could be changed to Pt16_fcc or Pt16_bcc respectively. Please do not use spaces in folders or file-names.

If you, for example, have calculations with different facets, you can also split them into separate folders, run ``cathub organize -f <facet>``, and them merge the organized folders together afterwards with :code:`cp -R organized1 organized2`.


cathub make-folders
...................
An alternative to ``cathub organize``. This tool will create the right folder structure for you, but you must dump your files yourself.

To learn about the make-folders command type::

  cathub make-folders --help

Then create a folder in your user-name, 'cd' into it and type::

  cathub make-folders --create-template TEMPLATE_NAME

This will create a template (txt/yaml) file, that you should update with your publication and reaction info. The template should look similar to this::

    reactions:
    -   reactants: [2.0H2Ogas, -1.5H2gas, star]
        products: [OOHstar@top]
    -   reactants: [CCH3star@bridge]
        products: [Cstar@hollow, CH3star@ontop]
    -   reactants: [CH4gas, -0.5H2gas, star]
        products: [CH3star@ontop]
    journal: JACS
    year: 2017
    email: winther@stanford.edu
    number: 1
    crystal_structures: [fcc, hcp]
    volume: 1
    DFT_functionals: [BEEF-vdW, HSE06]
    authors: ['Doe, John', 'Einstein, Albert']
    pages: 23-42
    publisher: ACS
    doi: 10.NNNN/....
    title: "Fancy title"
    bulk_compositions: [Pt]
    DFT_code: Quantum Espresso
    facets: ['111']


Consult :code:`cathub make-folders --help` again for detailed instructions on how to specify the types of reactions and surfaces.

Then type::

   cathub make-folders <TEMPLATE>

And your folders will be created. You can check that they look right with :code:`tree -F <FOLDER>`. The template above will produce the following folder tree::

  $tree -F MontoyaChallenge2015/

  MontoyaChallenge2015
  ├── Quantum\ Espresso/
  │   └── BEEF-vdW/
  │       ├── Co_fcc/
  │       │   ├── 111/
  │       │   │   ├── 0.5H2gas_star__Hstar@bridge/
  │       │   │   │   ├── MISSING:H_slab
  │       │   │   │   └── MISSING:TS?
  │       │   │   ├── 0.5H2gas_star__Hstar@fcc/
  │       │   │   │   ├── MISSING:H_slab
  │       │   │   │   └── MISSING:TS?
  │       │   │   ├── 0.5H2gas_star__Hstar@hollow/
  │       │   │   │   ├── MISSING:H_slab
  │       │   │   │   └── MISSING:TS?
  │       │   │   ├── 0.5H2gas_star__Hstar@ontop/
  │       │   │   │   ├── MISSING:H_slab
  │       │   │   │   └── MISSING:TS?
  │       │   │   ├── 0.5N2gas_0.5H2gas_star__NHstar@bridge/
  │       │   │   │   ├── MISSING:NH_slab
  │       │   │   │   └── MISSING:TS?
  │       │   │   ├── 0.5N2gas_0.5H2gas_star__NHstar@hollow/
  │       │   │   │   ├── MISSING:NH_slab
  │       │   │   │   └── MISSING:TS?
  │       │   │   ├── 0.5N2gas_star__Nstar@hollow/
  │       │   │   │   ├── MISSING:N_slab
  │       │   │   │   └── MISSING:TS?
  │       │   │   └── MISSING:empty_slab
  │       │   └── MISSING:Co_fcc_bulk
  │       └── gas/
  │           ├── MISSING:H2_gas
  │           └── MISSING:N2_gas
  └── publication.txt


Then add your atomic structure output files to the right folders. The files can be in any format that ASE can read, and must contain the total potential energy from the calculation. ASE trajectory (.traj) files are generally preferred. If you're using Vasp, please add your OUTCAR files as ``<name>.OUTCAR``. Your structures will include the adsorbed atoms/molecules, empty slabs, and gas phase species for your reactions. Also, if you have done calculations for the bulk geometries, please include them as well. All gas phase species involved must be added to the ``<publication>/<dft code>/<dft functional>/gas/`` folder. Also, notice that dummy files named ``MISSING:..`` have been placed in the folders, to help you determine the right location for your files. It's recommended to write a script to transfer files from one folder-structure to another in a systematic way, for example using :code:`shutils.copyfile('/path/to/initial/file', '/path/to/final/file')`.


Reading into database
......................
After adding all your structure files (or after running cathub organize), read your structures into a local database file with the command::

  cathub folder2db <FOLDER>

If anything is wrong with your files, or anything is missing, you should receive appropriate error messages. When reading of the folder is complete, a table with a summary with reaction energies will be printed in you terminal. Please verify that everything looks right. Also, a database file has been written at ``<FOLDER>/<DBNAME>.db``.

**SLAC affiliates** can upload directly to the server by typing::

  cathub db2server <DBNAME>.db

and follow the feedback in the terminal.

**Contributors outside SLAC** should instead send the ``<DBNAME>.db`` file by email to winther@slac.stanford.edu and the Catalysis-Hub team will handle the upload on your behalf.

In both cases, your data will not be made accessible from catalysis-hub.org until you have approved it. Please include the email address you defined above in your message.
