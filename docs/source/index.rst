.. _`index`:

====================
IMAS-ParaView Manual
====================

**IMAS-ParaView** is a comprehensive tool for visualizing
`IMAS <https://imas-data-dictionary.readthedocs.io/en/latest/>`_ data within `ParaView <https://www.paraview.org/>`_. 
It provides a collection of ParaView plugins that can visualize different kinds of
IMAS data, ranging from multi-dimensional 
`GGD (Generalized Grid Description) <https://imas-data-dictionary.readthedocs.io/en/latest/ggd_guide/doc.html>`_ 
meshes to non-GGD structures like machine descriptions and diagnostic geometries.

Get Started
===========

.. grid:: 1 2 3 3
   :gutter: 3
   :margin: 2 2 2 2

   .. grid-item-card::
      :link: installing
      :link-type: ref

      :si-icon:`material/download` **Installation**

      Get started by installing **IMAS-ParaView**.

   .. grid-item-card::
      :link: usage
      :link-type: ref

      :si-icon:`material/puzzle` **Plugins**

      Learn which **ParaView plugins** are available and how to use them.

   .. grid-item-card::
      :link: training
      :link-type: ref

      :si-icon:`material/school` **Training**

      Step-by-step tutorials on how to use **IMAS-ParaView**.

Gallery
=======

Explore what others have created with IMAS-ParaView.

.. jinja:: gallery_ctx

   .. div::
      :name: all-gallery-source
      :style: display:none

      .. grid:: 1 2 3 3
         :gutter: 2

      {% for entry in entries %}
         .. grid-item-card::
            :img-top: /gallery/examples/{{ entry.dir_name }}/{{ entry.image_name }}
            :link: {{ entry.name }}
            :link-type: ref

            {{ entry.title }}
      {% endfor %}

.. button-ref:: gallery
   :color: primary
   :outline:

   See all gallery examples →

.. toctree::
   :caption: Getting Started
   :maxdepth: 2
   :hidden:

   self
   installing

.. toctree::
   :caption: Gallery
   :maxdepth: 2
   :hidden:

   gallery/gallery
   
.. toctree::
   :caption: How To Use
   :maxdepth: 2
   :hidden:

   usage
   cli
   training/training

.. toctree::
   :caption: API docs
   :maxdepth: 1
   :hidden:

   api

.. toctree::
   :caption: Development
   :maxdepth: 1
   :hidden:

   code_style
   ci_config
   dev_guide
   license
