## Introduction

CatHub provides an interface to the Surface Reactions database on [Catalysis-Hub.org](http://www.catalysis-hub.org).

The module includes a command line interface that can be used to access and upload data. A short guide is given below. We refer to the [catalysis-hub documentation](http://docs.catalysis-hub.org/en/latest/tutorials/upload.html) for details on how to submit data.

## Using the cathub cli

Run `cathub` from the command line:

    cathub --help

or with any of its sub-commands:

    cathub reactions --help

## Examples

Querying the Surface Reactions database in Python:

    from cathub.cathubsql import CathubSQL

    # To get data on catalysis-hub.org
    db = CathubSQL()

    # Data from local cathub .db file
    db = CathubSQL('filename.db')

Get reactions in pandas dataframe:

    dataframe = db.get_dataframe(pub_id='PengRole2020',
                                 include_atoms=False,
                                 include_atoms=True,  # include atoms in dataframe
                                 #include_atoms='PengRole2020.db',  # save atoms to local db
                                 reactants=['COgas'],
                                 products=['COstar'],
                                 elements=['Cu', 'Al'],
                                 #surface_composition='Cu', # match specific composition
                                 facet = '100'
                                 )

Get atomic structure separately:

    # Get atoms for one reaction_id taken from dataframe
    atoms_list = db.get_atoms_for_reaction(reaction_id)

    # Get atoms for entire dataset
    atoms_list = db.get_atoms_for_publication(pub_id='PengRole2020')


Quick view of atomic structures on Catalysis Hub with ase db CLI:

    cathub ase 'CuAg pub_id=PengRole2020'

## Uploading data

Organizing a general folder into a structured folder:

    cathub organize <foldername> -a <ads1,ads2> -c <dft-code> -x <xc-functional> -f <facet> -S <crystal structure>

New: organize in interactive manner to update adsorbate name, site and facet on the run:

    cathub organize < foldername > -I ...

As an alternative to cathub organize, create an empty organized folderstructure for dropping files yourself. First create a template and edit it, then create the folders.

    cathub make_folders --create-template <template>
    cathub make_folders <template>

Reading folders into a local .db file:

    cathub folder2db <foldername>

Sending the data to the Catalysis Hub server:

    cathub db2server <dbfile>
