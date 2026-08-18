## Introduction

CatHub provides an interface to the Surface Reactions database on [Catalysis-Hub.org](http://www.catalysis-hub.org).

The module includes a command line interface that can be used to access and upload data. A short guide is given below. We refer to the [catalysis-hub documentation](http://docs.catalysis-hub.org/en/latest/tutorials/upload.html) for details on how to submit data.

## Using the cathub cli

Run `cathub` from the command line:

    cathub --help

or with any of its sub-commands:

    cathub reactions --help

## Authentication

Queries to `api.catalysis-hub.org/graphql` require an API key. Pass it directly or set the `CATHUB_API_KEY` environment variable:

    export CATHUB_API_KEY=your_key_here

## Examples

Querying the Surface Reactions database in Python:

    from cathub.query import CathubQuery

    client = CathubQuery(api_key='your_key_here')
    dataframe = client.get_reactions(pub_id='PengRole2020')

If `CATHUB_API_KEY` is set in the environment, the `api_key` argument can be omitted:

    client = CathubQuery()
    dataframe = client.get_reactions(pub_id='PengRole2020')

Filtering results:

    dataframe = client.get_dataframe(
        pub_id='PengRole2020',
        surface_composition='Cu',
        facet='100',
    )

The returned DataFrame has snake_case columns (`chemical_composition`, `reaction_energy`, `pub_id`, etc.), an `equation` column, and `atoms_name`/`atoms_id` list columns from the associated structures.


## Uploading data

Organizing a general folder into a structured folder:

    cathub organize <foldername> -a <ads1,ads2> -c <dft-code> -x <xc-functional> -f <facet> -S <crystal structure>

New: organize in interactive manner to update adsorbate name, site and facet on the run:

    cathub organize <foldername> -I ...

As an alternative to cathub organize, create an empty organized folderstructure for dropping files yourself. First create a template and edit it, then create the folders.

    cathub make-folders --create-template <template>
    cathub make-folders <template>

Reading folders into a local .db file:

    cathub folder2db <foldername>

Inspecting local local .db file:

    cathub show-reactions <dbfile>

Sending the data to the Catalysis Hub server:

    cathub db2server <dbfile>
