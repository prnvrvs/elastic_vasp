elastic_vasp
============

.. raw:: html

   <section class="ev-hero">
     <div class="ev-hero__eyebrow">elastic_vasp documentation</div>
     <h1>Elastic constants from VASP in four steps.</h1>
     <p>
       Give it a relaxed POSCAR. elastic_vasp tells you which calculations to run,
       writes the strained structures, and turns the finished VASP runs into the
       6&times;6 stiffness matrix &ndash; with checks at every step.
     </p>
     <div class="ev-hero__actions">
       <a class="ev-btn ev-btn--primary" href="quickstart.html">Start here</a>
       <a class="ev-btn" href="installation.html">Install</a>
       <a class="ev-btn" href="tutorial.html">Watch the video</a>
     </div>
   </section>

The whole method
----------------

.. list-table::
   :widths: 8 42 50
   :header-rows: 0

   * - **1**
     - **Check** your structure
     - ``elastic-vasp check`` – space group, and which workflow to use
   * - **2**
     - **Generate** the strained structures
     - ``elastic-vasp pre cubic`` (or ``hexagonal``, …, or ``auto``)
   * - **3**
     - **Run VASP** in every ``strain_*`` folder
     - your usual INCAR (``ISIF = 2``), KPOINTS and POTCAR
   * - **4**
     - **Analyse**
     - ``elastic-vasp post cubic`` – stiffness matrix, stability, averages

Every step is also a Python function (``check_symmetry()``, ``cubic_pre()``,
``cubic_post()``, …).

Where to go
-----------

* **New here?** :doc:`installation`, then the :doc:`quickstart` – one complete run in a few
  minutes of reading. Prefer watching? See the :doc:`tutorial`.
* **Want to see real numbers?** The :doc:`worked_example` runs fcc Al through VASP.
* **Your own material?** Follow the :doc:`User guide <guide/check>` page by page – it is
  written in the order you work.
* **Something specific?** The :doc:`command line <guide/cli>`, the
  :doc:`Python API <api/index>`, or the :doc:`theory`.

.. toctree::
   :maxdepth: 1
   :caption: Getting started
   :hidden:

   installation
   quickstart
   tutorial
   worked_example
   examples

.. toctree::
   :maxdepth: 1
   :caption: User guide – in the order you work
   :hidden:

   guide/check
   guide/systems
   guide/inputs
   guide/workflow
   guide/postprocessing
   guide/visualise

.. toctree::
   :maxdepth: 1
   :caption: Special cases
   :hidden:

   guide/orientation
   guide/auto

.. toctree::
   :maxdepth: 1
   :caption: Reference
   :hidden:

   guide/cli
   api/index
   theory
   release_notes
   legacy
   credit
   maintainer
