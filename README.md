## Introduction

CatHub provides an interface to the Surface Reactions database on [Catalysis-Hub.org](http://www.catalysis-hub.org).

The module includes a command line interface that can be used to access and upload data. A short guide is given below. We refer to the [catalysis-hub documentation](http://docs.catalysis-hub.org/en/latest/tutorials/upload.html) for details on how to submit data.

## Using the cathub cli

Run `cathub`, like so

    cathub --help

or with any of its sub-commands, like so

    cathub reactions --help

## Examples

Querying the Surface Reactions database:

    cathub reactions -q reactants=CO -q chemicalComposition=~Pt

    cathub publications -q title=~Evolution -q year=2017

Querying atomic structures on Catalysis Hub with ase db:

    cathub ase --args AgSr --gui

## Uploading data

Organizing a general folder into a structured folder:

    cathub organize <folderame> -a <ads1,ads2> -c <dft-code> -x <xc-functional> -f <facet> -S <crystal structure>

As an alternative to cathub organize - create an empty organized folderstructure for dropping files yourself. First create a template and edit it, then create the folders.
    cathub make_folders --create-template <template>
    cathub make_folders <template>

Reading folders into a local .db file:

    cathub folder2db <foldername> --userhandle <slack-username or gmail-address>

Sending the data to the Catalysis Hub server:

    cathub db2server <dbfile>
