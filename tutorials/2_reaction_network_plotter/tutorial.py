import matplotlib
import os
import sys
from cathub.reaction_networks import *

# %matplotlib inline

db_file = 'TangRevised2018.db'
df = db_to_df(db_file)

# CREATE REACTION NETWORK AND PLOT
op = ReactionNetwork(db_file)
ints = ['COgas','COstar','COHstar','Cstar','CHstar', 'CH2star']
plot = op.plot_network(intermediate_list=ints)
plot.show()
