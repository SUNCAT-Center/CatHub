## CatHub Tutorial

This tutorial is developed for the ACS Fall 2022 Workshop "Open-Source Software for Kinetics, Chemical Networks, & Reactor Modeling". Here, you will learn how to access catalysis-hub.org data via the Python interface. CatHub provides an interface to the Surface Reactions database on [Catalysis-Hub.org](http://www.catalysis-hub.org). The module includes a command line interface (in your terminal) as well as a Python interface to access and upload data.

# Installing cathub
To install CatHub use pip:

    pip3 install git+https://github.com/SUNCAT-Center/CatHub.git --upgrade --user

which will install CatHub and all their dependencies.

To test that the cathub cli is working, start by typing in your terminal:

    $ cathub --help

and you should see a list of subcommands. If it’s not working you probably have to add the installation path to PATH in your `~/.bashrc`. This would typically be `export PATH=~/.local/bin:${PATH}` for Linux, and `export PATH~/Library/PythonX.Y/bin:${PATH}` for Mac.

