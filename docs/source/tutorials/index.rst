Tutorials
=========

Searching for reaction energies
--------------------------------

Here you will learn how to use the webpage to search for chemical reactions. Go to http://www.catalysis-hub.org/energies to use the build in search function. This feature corresponds to `CatApp v2.0` where you can access reaction energies and barriers for different surfaces.

1) Type something in the reactants field and press the search button (for example CH3CH2OH*). How much is returned? Try to decrease the number of reactions by filling out the :code:`Products`, :code:`Surface` and :code:`Facet` fields.

2) Look at the atomic structures associated with a reaction of choice. Hint: press the cube icon to the left of the reaction. (Note: a few of the reactions doesn't have structures associated with them).

3) Interested in knowing how the data is pulled from the database backend? Click the :code:`GRAPHQL QUERY` button below the list of reactions, and see how you search looks from the backend interface.


GraphQL
-------

These tutorials will focus on how to use the GraphQL interface to the database. You might also want to read the documentation at http://catalysis-hub.readthedocs.io/en/latest/topics .


Go to http://api-catalysis-hub.slac.stanford.edu/graphql to get started.

A simple query
...............

Type your query in the left panel. In order to perform queries on the `reactions` table try this::

   {reactions(first:2)}

And type `command + return` to see the result in the right panel. This should return the `id` of the first two reactions in the database. Notice that the left hand side is updated as well.

A general note: Always include one of the 'first' or 'last' field in your query to limit the number of results. (Otherwise things could get slow!).


Now try to add more columns after `id` and see what happens. For example::

  {reactions(first:2) {
    edges {
      node {
        id
	Equation
        chemicalComposition
     	reactionEnergy
      }
    }
  }}


For a complete list of all the tables in the database, and the associated columns, see the `Docs` tab on the top right of the GraphiQL page. There is also a schema overview posted at  http://docs.catalysis-hub.org/reference/schema.html .

In order to make selections on the result, add more fields after :code:`(first:2)`. For example::

   {reactions(first:2, reactants: "CO", chemicalComposition: "Pt")}

Notice that it's possible to construct queries for all the existing columns in a table.


Searching for publications
..........................

1) Find all titles and DOI's from publications with year=2017.

2) How many publications are there with year>2015?

3) How many publications are authored by Thomas Bligaard? Hint: use the pubtexsearch field.
   You can list the total number of results using the `totalCount` field::

     {publications {
       totalCount
       edges {
         node {
           id
           authors
         }
       }
     }}


   Verify that you get the same result by using :code:`(authors: "~bligaard")`

4) Find Michal Bajdich's paper with "Oxygen Evolution" in the title.



Using the reactions table
.........................
1) Order all the reactions with respect to increasing reaction energy and print out the first 100 results.

2) Find the reactions with the lowest activation energy. Hint: take care of 'null' results by requesting that the activation energy should be > 0.


3) Find the number of reactions with H2O on the left hand side and OH on the right hand side.
   How many distinct reactions does that give rise to?
   What happens when you add the state (star or gas) to your query?
   What happens when you add :code:`chemicalComposition: "~"`


4) Chose a few of the reactions from the query before, and get all the chemical formula of the atomic structures beloning to them
   Hint: you can call the 'systems' table inside the 'reactions' table.


5) Find the publication year for the first 10 reactions in exercise 1)


Combining tables
....................
1) Find Julia Schuman's recent paper, and list all the reactions belonging to the paper. Hint: you can either go through the publication table::

     {publications(first:10) {
       edges {
         node {
           id
	   title
	   pubId
	   authors
	   reactions
	   }
	 }
     }}

   or use the `pubId` field to query directly on the reactions table to make additional queries::

     reactions(pubId: "")}

2) Using the :code:`(pubID:)` solution suggested above, list all the distinct
   a) reactions
   b) surfaces
   from Julia's publication.

3) Chose one of Julia's reactions and find the `aseId` of the empty slab. Hint: It has :code:`"name"="star"` in the `reactionSystems` table.
   Copy the aseId and use it to find all the reactions that are linked to that particular empty slab.

Custom GraphQL clients
----------------------

Calling the Backend from the Command Line
.........................................


This is easy using `curl` or `wget`::

    curl -XPOST https://api-catalysis-hub.slac.stanford.edu/graphql --data 'query={systems(last: 10 ) {
      edges {
        node {
          energy Cifdata
              }
      }
    }}'

Check out `jq <https://stedolan.github.io/jq/manual/>`_, `yaml2json <https://www.npmjs.com/package/yaml2json>`_, `json2yaml <https://www.npmjs.com/package/json2yaml>`_, or `glom <http://glom.readthedocs.io/en/latest/tutorial.html>`_ for terse json processing on the command-line. Try this::

    curl -XPOST https://api-catalysis-hub.slac.stanford.edu/graphql --data 'query={systems(last: 10 ) {
      edges {
        node {
          energy Cifdata
              }
      }
    }}' | jq '.data.systems.edges[].node.Cifdata' | sed -e 's/"//g' | split - -l 1 structure_ -d 

    sed -i 's/\\n/\n/g' structure_*

to write structures into many files. Or try this::

    curl -XPOST https://api-catalysis-hub.slac.stanford.edu/graphql --data 'query={
      reactions( reactants:"CO") {
        totalCount
        edges {
          node {
            chemicalComposition
            reactionEnergy
            sites
          }
        }
      }
    }
    ' | jq -r '.data.reactions.edges[].node | [.reactionEnergy,.chemicalComposition, .sites ] | @csv' 


for creating a CSV output.


Calling the Backend from a Python Script
........................................

Write a short python script with a GraphQL query of your choice. The script should look something like this::

     import requests
     import pprint

     root = 'http://api-catalysis-hub.slac.stanford.edu/graphql'

     query = \
     """
     {}
     """

     data = requests.post(root, {'query': query}).json()
     pprint.pprint(data)


And see the result printed in the terminal. How would you like to save the data?

Calling the Backend from a Perl Script
......................................

Write a short Perl script with a GraphQL query of your choice. This result could look something like this

.. code-block:: perl

   #!/usr/bin/env perl

   require LWP::UserAgent;

   my $uri = 'http://api-catalysis-hub.slac.stanford.edu/graphql';
   my $json = <<'EOF';
   {"query": "{
     reactions(first: 10) {
       edges {
         node {
           Equation
           chemicalComposition
           reactionEnergy
         }
       }
     }
   } "}
   EOF

   # remove newlines
   $json =~ s/(\n|\r)//g;

   my $req = HTTP::Request->new( 'POST', $uri );
   $req->header( 'Content-Type' => 'application/json' );
   $req->content( $json );

   my $lwp = LWP::UserAgent->new;
   my $res = $lwp->request( $req );

   print $res->content . "\n";

Calling the Backend from JavaScript
...................................

.. code-block:: js

   #!/usr/bin/env node

   var axios = require('axios');

   axios.post('http://api-catalysis-hub.slac.stanford.edu/graphql',
     {query:`
     {
     reactions(first: 10) {
       edges {
         node {
           Equation
           chemicalComposition
           reactionEnergy
         }
       }
     }
   }
       `}
   ).then(function(response){
     console.log(JSON.stringify(response.data))
   })


Calling the Backend from Coffee Script
......................................

.. code-block:: coffee

    #!/usr/bin/env coffee

    axios = require 'axios'

    axios.post 'http://api-catalysis-hub.slac.stanford.edu/graphql', {query:"{
      reactions(first: 10) {
        edges {
          node {
            Equation
            chemicalComposition
            reactionEnergy
          }
        }
      }
    }
        "}
    .then (response) ->
      console.log JSON.stringify response.data


Connecting to the database server with psql
.............................................

This exercise requires that you have postgreSQL installed, so you can use the `psql` terminal client.
Also you need the password for the `catvisitor` user, or optionally your own user account. Contact Kirsten Winther at winther@stanford.edu for question.

Type into the terminal::

  psql --host=catalysishub.c8gwuc8jwb7l.us-west-2.rds.amazonaws.com
  --port=5432 --username=catvisitor --dbname=catalysishub

And write the password when prompted.

Now you can start writing SQL statements directly against the database server. Try for example::

  SELECT title, year from publication LIMIT 10;

and see the output. Please use the :code:`LIMIT` clause to limit the number of results, or specify :code:`id=int`. See https://www.postgresql.org/docs/9.6/static/index.html for documentation on the SQL language and postgres.


Connecting to the database with ASE db
--------------------------------------
For this exercise you need to have a recent version of ASE installed. See https://wiki.fysik.dtu.dk/ase/install.html .

1) Now use the ASE cli to connect. Type this in the terminal (with an updated DB_PASSWORD)::

     ase db postgresql://catvisitor:$DB_PASSWORD@catalysishub.
     c8gwuc8jwb7l.us-west-2.rds.amazonaws.com:5432/catalysishub Pt3Co


(Note: this query is probably going to take some time. We're still working on optimizing the ASE database part.)


2) Write a python script to connect via ase.db.connect. Hint: the connect() function will take the same server URL as used in the previous exercise.

   You can now use the select() function to make queries against the database. See https://wiki.fysik.dtu.dk/ase/ase/db/db.html for documentation.



Using the CatHub cli
--------------------

This CLI will be used for collecting data from remote sources.

.. toctree::

   :maxdepth: 1
   cathubcli.md


IPython Notebook tutorials for API usage
----------------------------------------

The following noteboks demonstrate how make interactive use of the `Catalysis-Hub.Org API <http://api-catalysis-hub.slac.stanford.edu/>`_. The recommend of using these Notebook is `Jupyter Lab <http://jupyterlab.readthedocs.io/>`_.

`GraphQL Querying <http://nbviewer.jupyter.org/gist/mhoffman/556332aaac0e7e11769ce28848b6b721>`_
................................................................................................


`Retrieve ASE Atoms object through GraphQL <https://gist.github.com/mhoffman/c418acb6b3f928eb4ef71f4001d120d9>`_
................................................................................................................


`Prototype Search <http://nbviewer.jupyter.org/gist/mhoffman/2f680eb90a5531fe6e3e588f2c835775>`_
.................................................................................................


`Search Structures and Cut Slabs <http://nbviewer.jupyter.org/gist/mhoffman/e0e9edf6771c3c0c5838043a022e0857>`_
................................................................................................................


`Get Spacegroup and Wyckoff Sites from Structure <http://nbviewer.jupyter.org/gist/mhoffman/c483dece941930fbb785f019893a1c15>`_
................................................................................................................................


Partner Projects
-------------------------

The catalysis-hub.org database is large and its possibilities are manifold.
Fortunately there are two partner projects that are best digested together:
`CatKit <http://catkit.readthedocs.io/en/latest/>`_ and
`CatLearn <http://catlearn.readthedocs.io/en/latest/>`_.
