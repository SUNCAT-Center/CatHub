# CatHub Tutorial

This tutorial is developed in connection to the SUNCAT Summer Institute 2023.  CatHub provides an interface to the Surface Reactions database on [Catalysis-Hub.org](http://www.catalysis-hub.org). The module includes a command line interface (in your terminal) as well as a Python interface to access and upload data. In this tutorial you will learn how to access catalysis-hub.org data via the Python interface.

## Installing cathub
To install CatHub use pip:

    $ pip3 install git+https://github.com/SUNCAT-Center/CatHub.git --upgrade --user

which will install CatHub and all their dependencies.

To test that the cathub cli is working, start by typing in your terminal:

    $ cathub --help

and you should see a list of subcommands. If it’s not working you probably have to add the installation path to PATH in your `~/.bashrc`. This would typically be `export PATH=~/.local/bin:${PATH}` for Linux, and `export PATH~/Library/PythonX.Y/bin:${PATH}` for Mac.


## Demo 1: Querying adsorption energetics from Catalysis-hub
In this demonstration you will learn how to fetch adsorption energies from Catalysis-hub.org

Start by familiarizing yourself with the datasets on the main webpage http://www.catalysis-hub.org/publications

Data is divided into distinct datasets, usually belonging to a particular publication. A dataset is queried based on the "pub_id"  which is a a unique dataset id constructed from title, first author name and publication year. In this example we will start by using our recent dataset https://www.catalysis-hub.org/publications/SainiElectronic2022


### Fetching data in Python
To connect to the catalysis-hub.org server in your Python script, start by importing the cathub SQL interface, creating a database connection to the catalysis-hub server:


```python
from cathub.cathubsql import CathubSQL
db = CathubSQL()
```

Then use the `get_dataframe()` method to query adsorption energy data into a pandas dataframe ( Basic syntax for Pandas is found [here](https://pandas.pydata.org/pandas-docs/stable/getting_started/intro_tutorials/03_subset_data.html#min-tut-03-subset) )

Inspect the dataframe by printing it to your terminal/notebook. Main columns consists of the <b> chemicalComposition </b> (chemical formula of the total slab), <b>surface_composition</b> (reduced chemical composition with surface specific tags), <b> equation </b> (equation for the reaction), and <b>reaction_energy </b> (which can also be an adsorption energy).


```python
pub_id = "SainiElectronic2022"
dataframe = db.get_dataframe(pub_id=pub_id)
print(dataframe)
```

### Inspecting the data
To continue the analysis, please save the dataframe into a pickle file on your local workspace. Now you can examine the your local file without pulling from the server.

Select specific columns for a better visualization of the data


```python
dataframe.to_pickle(pub_id + '.pickle')

import pandas
dataframe = pandas.read_pickle(pub_id + '.pickle')

print(dataframe.columns)
print(dataframe[['chemical_composition', 'surface_composition','facet', 'equation', 'reaction_energy']].to_markdown())
```

### Filtering the data

Data can be filtered using pandas syntax, for example selecting a specific chemical reaction like this:



```python
equation = '0.5O2(g) + * -> O*'  

loc1 = (dataframe["equation"] == equation) & (dataframe['surface_composition'].str.contains('Au5'))
loc2 = (dataframe["equation"] == equation) & (dataframe['surface_composition'].str.contains('Cu5'))
dataframe_1 = dataframe[loc1]
dataframe_2 = dataframe[loc2]

import pylab as p
p.figure(figsize=(10,8))
p.hist(dataframe_1['reaction_energy'].values, bins=20, alpha=0.5, label='Au5M')
p.hist(dataframe_2['reaction_energy'].values, bins=20, alpha=0.5, label='Cu5M')
p.title('0.5O2(g) + * -> O*')
p.ylabel('Occurence')
p.xlabel('$\Delta$E (eV)')
p.legend()
p.show()
```

### Challenge: Refine query based on chemical composition, adsorbates and facet
It is also posible to filter the data already on the level of the cathub query, and to search across publications


```python
dataframe2 = db.get_dataframe(#reactants={'COgas': 1},
                              #products={'COstar': 1},
                              elements=['Cu', 'Al'], #contains Cu and Al
                              #surface_composition='Cu', # match specific composition
                              facet = '100'
                                )
print(dataframe2[['pub_id', 'surface_composition', 'reaction_energy']])
```
