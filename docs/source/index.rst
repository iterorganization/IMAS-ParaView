.. _`index`:

.. 
   Main "index". This will be converted to a landing index.html by sphinx. We
   define TOC here, but it'll be put in the sidebar by the theme

====================
IMAS-ParaView Manual
====================

IMAS-ParaView is a tool to convert GGD (Generalized Grid Description) structures to VTK formats, and back.
This is complemented by a number of Paraview plugins that can visualise both GGD and non-GGD IDS data in Paraview.


README
------

The README is best read on the `git page <https://github.com/IMAS-ParaView/>`_.

Manual
------

.. toctree::
   :caption: Getting Started
   :maxdepth: 2

   self
   installing
   
.. toctree::
   :caption: How To Use
   :maxdepth: 2

   usage
   cli

.. toctree::
   :caption: API docs
   :maxdepth: 1

   api

.. toctree::
   :caption: Development
   :maxdepth: 1

   code_style
   ci_config


LICENSE
-------

.. literalinclude:: ../../LICENSE.txt
   :language: text


Sitemap
-------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
