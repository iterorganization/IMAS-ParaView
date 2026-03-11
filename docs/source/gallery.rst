.. _`gallery`:

=======
Gallery
=======

The gallery showcases the capabilities of what you can do with the IMAS-ParaView plugins.
You can click on the images to enlarge them and to display their descriptions.

.. grid:: 1 2 3 3
   :gutter: 2

   .. grid-item-card::

      .. thumbnail:: images/gallery/jorek.gif

         Animation of the electron temperature and wall currents in a JOREK simulation.
         A step-by-step tutorial on how to re-create this animation using the IMAS-ParaView
         tools has been provided in the :ref:`JOREK tutorial <training_jorek>`.
         Data provided by J. Artola, using the following URI:

         ``imas:hdf5?user=public;pulse=112111;run=2;database=ITER_DISRUPTIONS;version=4``

Contributing to the Gallery
---------------------------

Have a cool visualisation made with IMAS-ParaView? We'd love to feature it here!
Contributions are welcome via a pull request on the
`IMAS-ParaView GitHub repository <https://github.com/iterorganization/IMAS-ParaView>`_.

Please follow the steps below to contribute:

#. Create a `fork <https://github.com/iterorganization/IMAS-ParaView/fork>`_ of the IMAS-ParaView repository.
#. Add your image or animation to the ``docs/source/images/gallery/`` directory.
#. Edit ``docs/source/gallery.rst`` and add a ``.. grid-item-card::`` entry containing a 
   ``.. thumbnail::`` pointing to your image, along with a short description.
   Ensure you have permission to use the image and properly credit the data owners.
#. Open a pull request from your fork back to the main repository on the ``develop`` branch.
   A maintainer will review and merge it.
