Installation
============

Requirements
------------

* Python 3.10 or newer
* `NumPy <https://numpy.org>`_ and `matplotlib <https://matplotlib.org>`_
* `ASE <https://ase-lib.org/>`_ 3.23 or newer – reads and writes POSCAR files.
  Older ASE versions silently ignore negative and per-axis POSCAR scale factors, so
  elastic_vasp refuses to run with them.
* `spglib <https://spglib.readthedocs.io/>`_ 2.4 or newer – finds the crystal symmetry for
  the orientation and symmetry checks and for :doc:`guide/auto`.

All of these are installed automatically. VASP itself is not needed by the package: it
only prepares the input structures and analyses the results.

.. important::

   This documentation describes the current code on GitHub (``main``). The release on PyPI
   is older and lacks most of it (command line, automatic mode, symmetry checks, …), so
   until the next release is published, install from GitHub:

   .. code-block:: bash

      pip install "git+https://github.com/prnvrvs/elastic_vasp.git"

From PyPI
---------

.. code-block:: bash

   pip install elastic-vasp

From source
-----------

.. code-block:: bash

   git clone https://github.com/prnvrvs/elastic_vasp.git
   cd elastic_vasp
   pip install -e .

To also run the tests or build this documentation:

.. code-block:: bash

   pip install -e ".[test,doc]"

Check the installation
----------------------

.. code-block:: bash

   python -c "import elastic_vasp; print(elastic_vasp.__version__)"

Branches
--------

The ``main`` branch holds the current code. The original code (version 0.0.4) is kept on
the ``legacy`` branch and on PyPI (``pip install "elastic-vasp==0.0.4"``) for reproducing
earlier results; see :doc:`legacy` for how to install it and its known issues.
