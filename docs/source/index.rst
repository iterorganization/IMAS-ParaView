.. _`index`:

====================
IMAS-ParaView Manual
====================

IMAS-ParaView is a comprehensive tool for visualizing and analyzing plasma physics
simulation data within ParaView. It provides a collection of ParaView plugins that
can load and visualize both GGD (Generalized Grid Description) and non-GGD IDS data.


.. grid:: 1 2 3 3
   :gutter: 3
   :margin: 2 2 2 2

   .. grid-item-card::
      :link: installing
      :link-type: ref

      .. raw:: html


         <div style="display: flex; align-items: center; gap: 0.75rem;">
            <svg xmlns="http://www.w3.org/2000/svg" width="32" height="32" viewBox="0 0 24 24"><!-- Icon from Google Material Icons by Material Design Authors - https://github.com/material-icons/material-icons/blob/master/LICENSE --><path fill="currentColor" d="M20 17H4V5h8V3H4c-1.11 0-2 .89-2 2v12a2 2 0 0 0 2 2h4v2h8v-2h4c1.1 0 2-.9 2-2v-3h-2z"/><path fill="currentColor" d="m17 14l5-5l-1.41-1.41L18 10.17V3h-2v7.17l-2.59-2.58L12 9z"/></svg>
           <strong>Installation</strong>
         </div>

      Get started by installing IMAS-ParaView on your system.

   .. grid-item-card::
      :link: usage
      :link-type: ref

      .. raw:: html

         <div style="display: flex; align-items: center; gap: 0.75rem;">
            <svg xmlns="http://www.w3.org/2000/svg" width="32" height="32" viewBox="0 0 24 24"><!-- Icon from Google Material Icons by Material Design Authors - https://github.com/material-icons/material-icons/blob/master/LICENSE --><path fill="currentColor" d="M21 14c0-.55-.45-1-1-1h-2v2h2c.55 0 1-.45 1-1m-1 3h-2v2h2c.55 0 1-.45 1-1s-.45-1-1-1m-4-5h-2c-1.1 0-2 .9-2 2h-1c-.55 0-1 .45-1 1v2c0 .55.45 1 1 1h1c0 1.1.9 2 2 2h2c.55 0 1-.45 1-1v-6c0-.55-.45-1-1-1M5 13c0-1.1.9-2 2-2h1.5c1.93 0 3.5-1.57 3.5-3.5S10.43 4 8.5 4H5c-.55 0-1 .45-1 1s.45 1 1 1h3.5c.83 0 1.5.67 1.5 1.5S9.33 9 8.5 9H7c-2.21 0-4 1.79-4 4s1.79 4 4 4h2v-2H7c-1.1 0-2-.9-2-2"/></svg>
           <strong>Plugins</strong>
         </div>

      Learn how to use the ParaView plugins

   .. grid-item-card::
      :link: training
      :link-type: ref

      .. raw:: html

         <div style="display: flex; align-items: center; gap: 0.75rem;">
         <svg xmlns="http://www.w3.org/2000/svg" width="32" height="32" viewBox="0 0 24 24"><!-- Icon from Google Material Icons by Material Design Authors - https://github.com/material-icons/material-icons/blob/master/LICENSE --><path fill="currentColor" d="M21 5c-1.11-.35-2.33-.5-3.5-.5c-1.95 0-4.05.4-5.5 1.5c-1.45-1.1-3.55-1.5-5.5-1.5S2.45 4.9 1 6v14.65c0 .25.25.5.5.5c.1 0 .15-.05.25-.05C3.1 20.45 5.05 20 6.5 20c1.95 0 4.05.4 5.5 1.5c1.35-.85 3.8-1.5 5.5-1.5c1.65 0 3.35.3 4.75 1.05c.1.05.15.05.25.05c.25 0 .5-.25.5-.5V6c-.6-.45-1.25-.75-2-1m0 13.5c-1.1-.35-2.3-.5-3.5-.5c-1.7 0-4.15.65-5.5 1.5V8c1.35-.85 3.8-1.5 5.5-1.5c1.2 0 2.4.15 3.5.5z"/><path fill="currentColor" d="M17.5 10.5c.88 0 1.73.09 2.5.26V9.24c-.79-.15-1.64-.24-2.5-.24c-1.7 0-3.24.29-4.5.83v1.66c1.13-.64 2.7-.99 4.5-.99M13 12.49v1.66c1.13-.64 2.7-.99 4.5-.99c.88 0 1.73.09 2.5.26V11.9c-.79-.15-1.64-.24-2.5-.24c-1.7 0-3.24.3-4.5.83m4.5 1.84c-1.7 0-3.24.29-4.5.83v1.66c1.13-.64 2.7-.99 4.5-.99c.88 0 1.73.09 2.5.26v-1.52c-.79-.16-1.64-.24-2.5-.24"/></svg>
           <strong>Training</strong>
         </div>

      Step-by-step tutorials on how to use IMAS-ParaView.

Gallery
=======

Explore what others have created with IMAS-ParaView.

.. jinja:: gallery_ctx

   .. grid:: 1 2 3 3
      :gutter: 2

   {% for entry in entries[:3] %}
      .. grid-item-card::
         :img-top: /gallery/examples/{{ entry.dir_name }}/{{ entry.image_name }}
         :link: {{ entry.name }}
         :link-type: ref

         {{ entry.title }}
   {% endfor %}

:ref:`See all gallery examples <gallery>`

.. toctree::
   :maxdepth: 2
   :hidden:

   installing
   usage
   cli
   training/training
   gallery/gallery
   api
   code_style
   ci_config
   dev_guide
