#CatHub Tutorial

This tutorial is developed in connection to the ACS Fall 2022 Symposium "Open-Source Software for Kinetics, Chemical Networks, & Reactor Modeling".  CatHub provides an interface to the Surface Reactions database on [Catalysis-Hub.org](http://www.catalysis-hub.org). The module includes a command line interface (in your terminal) as well as a Python interface to access and upload data. In this tutorial you will learn how to access catalysis-hub.org data via the Python interface.


## Installing cathub
To install CatHub use pip:

    pip3 install git+https://github.com/SUNCAT-Center/CatHub.git --upgrade --user

which will install CatHub and all their dependencies.

To test that the cathub cli is working, start by typing in your terminal:

    $ cathub --help

and you should see a list of subcommands. If it’s not working you probably have to add the installation path to PATH in your `~/.bashrc`. This would typically be `export PATH=~/.local/bin:${PATH}` for Linux, and `export PATH~/Library/PythonX.Y/bin:${PATH}` for Mac.

## Catalysis-hub background info
Familiarize yourself with the datasets on the Catalysis-hub webpage http://www.catalysis-hub.org/publications

Data is divided into distinct datasets belonging to a particular publication, that can be queries based on the "pub_id"  a unique dataset id constructed from title, first author name and publication year.

# Exercise 1: Scaling relations
Here you will learn how to fetch adsorption energies and plot scaling relations in python. We will use the example of O vs OH scaling from our 2022 transition metal oxides dataset https://www.catalysis-hub.org/publications/ComerUnraveling2022, publication link: https://pubs.acs.org/doi/10.1021/acs.jpcc.2c02381


To connect to the catalysis-hub.org server in your script, start by importing the cathub SQL interface, creating a database connection to the catalysis-hub server.

    from cathub.cathubsql import CathubSQL
    db = CathubSQL()

Then use the `get_dataframe()` method to query adsorption energy data into a pandas dataframe. Basic Syntax for Pandas is found [here](https://pandas.pydata.org/pandas-docs/stable/getting_started/intro_tutorials/03_subset_data.html#min-tut-03-subset):

    dataframe = db.get_dataframe(pub_id='ComerUnraveling2022')
    print(dataframe)

Inspect the dataframe by printing it to your terminal. Main columns consists of the chemicalComposition (chemical formula of the total slab), surface_composition (reduced chemical composition with surface specific tags), equation (equation for the reaction), reaction_energy (which can also be an adsorption energy)

To continue the analysis, please save the dataframe into a pickle file on your local workspace:

    dataframe.to_pickle('ComerUnraveling2022.pickle')

Now you can examine the your local file without pulling from the server:

    # db = CathubSQL()
    # dataframe = db.get_dataframe(pub_id='ComerUnraveling2022')
    # dataframe.to_pickle('ComerUnraveling2022.pickle')

    dataframe = pandas.read_pickle('ComerUnraveling2022.pickle')

Use your favorite python plotting module to plot the OH vs. O scaling (i.e. plotting OH vs O adsorption energies) Start by examining the unique chemical reactions and facets for the dataset, for example:

    print(dataframe["equation"].unique().tolist()) # Unique reactions
    print(dataframe["facet"].unique().tolist()) # unique facets


Using matplotlib/pylab you can plot the scaling relation like this:

    import pylab as p

    O_110 = dataframe[(dataframe["equation"] =='H2O(g) - H2(g) + * -> O*') & (dataframe['facet']=='110')]
    OH_110 = dataframe[(dataframe["equation"] =='H2O(g) - 0.5H2(g) + * -> HO*' ) &(dataframe['facet']=='110')]
    dataframe_together = O_110.merge(OH_110, on='surface_composition',
                                      suffixes=('_O', '_OH'))

    x_data = dataframe_together['reaction_energy_OH']
    y_data = dataframe_together['reaction_energy_OH']
    p.scatter(x_data, y_data)
    for i, txt in enumerate(dataframe_together['surface_composition']):
      p.gca().annotate(txt,
        (x_data[i],
         y_data[i]))
    p.title('O-OH Scaling relation')
    p.show()


Now try to repeat the plot choosing another facet from the dataset.

### Exercise 1 challenge: Plot the a scaling relationship for another dataset of your choice.
tip: remember to update chemical equation in the script if you choose other adsorbates.

# Exercise 2: Atomic structures

Use Python to query atomic structure for a dataset of choice (Please choose a smaller dataset with Nreactions < 500) to save time). In the example below I have chosen (https://www.catalysis-hub.org/publications/AraComputational2022):

    from cathub.cathubsql import CathubSQL
    from ase.visualize import view

    pub_id = 'AraComputational2022'

    db = CathubSQL()
    atoms_list = db.get_atoms_for_publication(pub_id=pub_id)
    dblocal = connect(pub_id + '.db')
    for atoms in atoms_list:
        dblocal.write(atoms)
    view(atoms_list)

You should now see ase gui open with several atomic structures. Also, inspect your db with the ASE from the command line:

    $ ase db <dbfile>
    $ ase gui <dbfile>

Next, try to query structures into the pandas dataframe where atomic structures and reaction energies are connected.

    db = CathubSQL()
    dataframe = db.get_dataframe(pub_id=pub_id,
                                 include_atoms=True)
    dataframe.to_pickle(pub_id + '_with_atoms.pickle')


You can view atoms for a specific reaction row by choosing a reaction row in the script below. Atoms objects include empty surface, surface with adsorbate, gas phase molecules and (for some datasets) the bulk geometry.

    dataframe = pandas.read_pickle(pub_id + '_with_atoms.pickle')
    print(dataframe[['chemical_composition', 'equation', 'atoms_name']])
    view(dataframe['atoms']["row_id"])

Notice that the "atoms_name" column contain names of geometries, so that specific types of structures can be fetch systematically. For example, to query only OH adsorption geometries, try this:

    atoms_list_OH = []
    for id, row in dataframe.iterrows():
        if not 'HOstar' in row['atoms_name']:
            continue
        index = row['atoms_name'].index('HOstar')
        atoms_list_OH += [row['atoms'][index]]

    view(atoms_list_OH)


### Exercise 2 challenge:
Query structures for a specific reaction directly from the server based
on the "reaction_id":

        atoms_list = db.get_atoms_for_reaction(reaction_id=???)

tip: The reaction_id's can be found in the pandas dataframe.


# Exercise 3: Experimental data

In this exercise we will access experimental datasets on the catalysis-hub server
(only available through Python API so far)

For this exercise you will only need the cathub cli. Go to your terminal and type:

    $ cathub exp

And you should get a pop-up opening in your browser. The table lists the
experimental datasets, where you can select the pub_id see more details:

    $ cathub exp <pub_id>
