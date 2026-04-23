.. _`training`:

Training Examples
=================

This section presents a step-by-step walkthrough of the examples from the IMAS-ParaView training held on June 17–18, 2025.

These steps assume that you have IMAS-ParaView up and running, for details on this, see :ref:`installing`. Furthermore, it is assumed that you have access to SDCC, as the data used for the examples is loaded from SDCC. It is advised to have some experience with IMAS-ParaView before starting this, so consider having a look at the instructions for :ref:`plugins`. It is recommended to stick to the order as presented here, as the later training sections introduce more advanced concepts.

At the start of each section, a ParaView ``.pvsm`` state file is also provided for your convenience. This can be loaded into ParaView and will automatically load the entire example setup for you. However, if it is your first time going through the examples, it is recommended to follow the step-by-step walkthrough manually, to become familiar with the workflow of IMAS-ParaView.

.. note::
   The ParaView state files have been tested with IMAS-ParaView version ``2.4.0``, 
   but may not work for different versions. If you get any errors when loading the state file, 
   it is recommended to revert your IMAS-ParaView version to ``2.4.0``. Please also 
   `report an issue <https://github.com/iterorganization/IMAS-ParaView/issues>`_ so the
   state files can be updated.


.. grid:: 1 2 3 3
   :gutter: 2

   .. grid-item-card::
      :img-top: /images/training/solps_electron_pressure.png
      :class-img-top: fixed-150height
      :link: training_solps
      :link-type: doc

      **1. SOLPS-ITER Case**

   .. grid-item-card::
      :img-top: /images/training/jorek.webp
      :class-img-top: fixed-150height
      :link: training_jorek
      :link-type: doc

      **2. JOREK Case**

   .. grid-item-card::
      :img-top: /images/training/machine_description.webp
      :class-img-top: fixed-150height
      :link: training_md
      :link-type: doc

      **3. Machine Description Case**

   .. grid-item-card::
      :img-top: /images/training/jintrac_mapper.png
      :class-img-top: fixed-150height
      :link: training_jintrac
      :link-type: doc

      **4. JINTRAC Case**

   .. grid-item-card::
      :img-top: /images/training/solps_batch1.png
      :class-img-top: fixed-150height
      :link: training_batch
      :link-type: doc

      **5. Batch Processing**

.. toctree::
   :maxdepth: 2
   :hidden:

   training_solps
   training_jorek
   training_md
   training_jintrac
   training_batch
