Apps
=====

HTTP Request to the Apps backend can be made in Python using e.g.
the `requests` library for GET requests, like so::

    #!/usr/bin/env python

    import json
    import pprint

    import requests

    url = 'http://api-catalysis-hub.slac.stanford.edu/apps/prototypeSearch/facet_search/'
    r = requests.get(url,
                      params={
                          'search_terms': [
                              'hollandite',
                          ],
                          'facet_filters': '["spacegroup:87"]'
                      }
                      )
    pprint.pprint(json.loads(r.content))

or for POST requests, like so::

    #!/usr/bin/env python

    import json
    import pprint

    import requests

    url = 'http://api-catalysis-hub.slac.stanford.edu/apps/prototypeSearch/get_structure/'

    r = requests.post(url,
        json={
            'parameters': '[3.1]',
            'species': '["S"]',
            'spacegroup': 221,
            })
    pprint.pprint(json.loads(r.content))

That is, every function below that is declared either as a GET or a POST request can
be translated into a corresponding HTTP request by replacing
every dot ('.') with a slash ('/') and passing in arguments either as params=... (GET)
or json=... (POST) .

Please refer to the following API documentation for details.


ActivityMaps
------------
.. automodule:: apps.activityMaps
    :members:
    :undoc-members:
    :show-inheritance:

AtoMl
-----
.. automodule:: apps.AtoML
    :members:
    :undoc-members:
    :show-inheritance:

.. automodule:: apps.AtoML.run_atoml
    :members:
    :undoc-members:
    :show-inheritance:

Bulk Enumerator
---------------
.. automodule:: apps.bulkEnumerator
    :members:
    :undoc-members:
    :show-inheritance:

CatKitDemo
----------
.. automodule:: apps.catKitDemo
    :members:
    :undoc-members:
    :show-inheritance:

Pourbaix diagrams
-----------------

.. automodule:: apps.pourbaix.run_pourbaix
    :members:
    :undoc-members:
    :show-inheritance:

Prototype Search
-----------------

.. automodule:: apps.prototypeSearch
    :members:
    :undoc-members:
    :show-inheritance:

Utilities
---------
.. automodule:: apps.utils
    :members:
    :undoc-members:
    :show-inheritance:
