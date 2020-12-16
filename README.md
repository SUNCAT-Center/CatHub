## Introduction

CatHub provides an interface to the Surface Reactions database on [Catalysis-Hub.org](http://www.catalysis-hub.org).

The module includes a command line interface that can be used to access and upload data. A short guide is given below. We refer to the [catalysis-hub documentation](http://docs.catalysis-hub.org/en/latest/tutorials/upload.html) for details on how to submit data.

## Using the cathub cli

Run `cathub`, like so

    cathub --help

or with any of its sub-commands, like so

    cathub reactions --help

## Examples

Querying the Surface Reactions database in Python:

    from cathub.cathubsql import CathubSQL
    db = CathubSQL() # All data on catalysis-hub.org
    db = CathubSQL('filename.db') # Data from local cathub .db file
    # Get reactions and structures
    dataframe = db.get_dataframe(pub_id='PengRole2020',
                                 include_atoms=False,
                                 #include_atoms='PengRole2020.db', # save to local db
                                 reactants=['COgas'],
                                 products=['COstar'],
                                 elements=['Cu', 'Al'],
                                 #surface_composition='Cu', # match specific composition
                                 facet = '100'
                                 )


    # Get atoms for one reaction_id taken from dataframe
    atoms_list = db.get_atoms_for_reaction(reaction_id)
    # Get atoms for entire dataset
    atoms_list = db.get_atoms_for_publication(pub_id='PengRole2020')

Querying the Surface Reactions database with the CLI:

    cathub reactions -q reactants=CO -q chemicalComposition=~Pt

    cathub publications -q title=~Evolution -q year=2017

Querying atomic structures on Catalysis Hub with ase db:

    cathub ase 'AgSr' --gui

## Uploading data

Organizing a general folder into a structured folder:

    cathub organize <folderame> -a <ads1,ads2> -c <dft-code> -x <xc-functional> -f <facet> -S <crystal structure>

As an alternative to cathub organize - create an empty organized folderstructure for dropping files yourself. First create a template and edit it, then create the folders.
    cathub make_folders --create-template <template>
    cathub make_folders <template>

Reading folders into a local .db file:

    cathub folder2db <foldername>

Sending the data to the Catalysis Hub server:

    cathub db2server <dbfile>
