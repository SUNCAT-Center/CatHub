## Introduction

CatHub provides an interface to the Surface Reactions database on [Catalysis-Hub.org](http://www.catalysis-hub.org][Catalysis-Hub.org).

The module includes a command line interface that can be used to access and upload data. A short guide is given below. We refer to the [catalysis-hub documentation](http://docs.catalysis-hub.org/en/latest/tutorials/upload.html) for details on how to submit data.

## Using the cathub cli

Run `cathub`, like so

    cathub --help

or with any of its sub-commands, like so

    cathub reactions --help

## Examples

Querying the Catalysis Hub database:

    cathub reactions -q reactants=CO -q chemicalComposition=~Pt

    cathub publications -q title=~Evolution -q year=2017

Reading folders into sqlite3 db file:

    cathub folder2db <foldername>

Sending the data to the Catalysis Hub server:

    cathub db2server <dbfile>
