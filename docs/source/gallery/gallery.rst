.. _`gallery`:

=======
Gallery
=======

The gallery showcases the capabilities of what you can do with the IMAS-ParaView plugins.
You can click on the images to enlarge them and to display their descriptions.

.. jinja:: gallery_ctx

   .. grid:: 1 2 3 3
      :gutter: 2

   {% for entry in entries %}
      .. grid-item-card::
         :img-top: /gallery/examples/{{ entry.dir_name }}/{{ entry.image_name }}
         :class-img-top: fixed-150height
         :link: {{ entry.name }}
         :link-type: ref

         {{ entry.title }}
   {% endfor %}

.. toctree::
   :maxdepth: 1
   :glob:
   :hidden:

   examples/*/index

Contributing to the Gallery
---------------------------

Have a cool visualisation made with IMAS-ParaView? We'd love to feature it here!
Contributions are welcome via a pull request on the
`IMAS-ParaView GitHub repository <https://github.com/iterorganization/IMAS-ParaView>`_.

Please follow the steps below to contribute:

#. Create a `fork <https://github.com/iterorganization/IMAS-ParaView/fork>`_ of the IMAS-ParaView repository.
#. Create a new directory in ``docs/source/gallery/examples``.
#. Create a ``description.yaml`` file in that directory (see template below).
#. Add an image showcasing the visualisation to the same directory (only ``.png``, ``.jpeg``, ``.jpg``, ``.gif``, or ``.webp`` file formats are supported)
#. Optionally add a ParaView state file (``.pvsm``) to the same directory, so other users can easily load your example.
#. Open a pull request from your fork back to the main repository on the ``develop`` branch.
   A maintainer will review and merge it.

Example ``description.yaml``:

.. code-block:: yaml

   title: My Visualization
   author: Your Name
   description: |
     A description of what the visualization shows.
     Can be multiple lines.
   uri: imas:hdf5?path=/path/on/sdcc # May also be a list of URIs
   link: www.example.com # Optional link (e.g. to dataset on Zenodo)
   imas_paraview_version: 2.3.0

