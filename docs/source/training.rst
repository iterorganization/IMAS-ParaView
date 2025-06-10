.. _`training`:

Training
========

This page explains the examples shown in the IMAS-ParaView training that took place 17th and 18th of June 2025. 
These steps assume that you have IMAS-ParaView up and running, for this see :ref:`installing`. It is advised to have some experience with IMAS-ParaView before starting this, so consider reading the instructions for :ref:`usage`.

SOLPS-ITER Case
---------------

TODO

JOREK Case
----------
For this example, we will visualize a JOREK disruption case available from the confluence page https://confluence.iter.org/display/IMP/The+JOREK+disruption+cases. 
We will visualize the electron temperature from the `plasma_profiles` IDS and corresponding current magnitude in the surrounding wall. We will create a video of this behaviour over time.

We following URI will be used:


.. code-block:: bash

   imas:hdf5?user=public;pulse=112111;run=2;database=ITER_DISRUPTIONS;version=4


.. |ico1| image:: images/rotate_axis.png

Loading the Electron Temperature
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In this subsection, we will load the JOREK grid and visualize the electron temperatures on this grid.

#. Navigate to Sources > IMAS Tools and select the JOREK Reader.
#. Enter the URI shown above in the ``Enter URI`` field of the reader plugin.
#. Select the ``plasma_profiles/1`` IDS in the IDS/Occurrence dropdown menu.
#. Select ``Apply`` to load the plasma profiles GGD grid.
#. After the GGD grid is loaded, bring the grid into view by aligning the viewpoint in the positive Y direction using the following button: |ico1|.
#. Select the ``Electrons Temperature`` from the attribute array selection window.
#. Select ``Apply`` to load the electron temperature values on the grid.
#. Select ``Electrons Temperature [eV]`` in the coloring dropdown to visualize the electron temperature.
#. Enable log scale coloring by editing the color map, and selecting ``Use Log Scale When Mapping Data To Colors``
#. Set the ``N plane`` to 3 and the ``Phi range`` from 0 to 180 degrees in the Bezier interpolation settings.

Loading the Wall Current
^^^^^^^^^^^^^^^^^^^^^^^^

In this subsection, we will load the wall currents of the simulation using the GGD Reader, and apply a clip mask.

#. Navigate to Sources > IMAS Tools and select the GGD Reader.
#. Enter the URI shown above in the ``Enter URI`` field of the reader plugin.
#. Select the ``wall`` IDS in the IDS/Occurrence dropdown menu.
#. Select ``Apply`` to load the wall grid.
#. Select the ``J_totel`` from the attribute array selection window.
#. Select ``Apply`` to load the current on the wall GGD grid.
#. Select ``Description_ggd J_total [A.m^-2]`` in the coloring dropdown and select the ``Magnitude`` to visualize the total wall current.
#. As the wall shows an enclosed surface, which is a bit hard to see, we will apply a clip filter to the wall grid. To do this, select the clip [TODO icon] filter.
#. Set the normal vector to ``0, -1, 0`` and select ``Apply`` to apply the filter
#. In order to distinguish between the wall currents and the electron temperature grid, we can change the color map wall currents. Edit the color map and select ``Select a color map from default presets`` and choose a different color map.

Temporal Interpolation and Animation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this subsection, we will create an animation of the loaded electron temperature and wall currents. To this end, we will interpolate the time basis of both plugins to a linearly spaced basis.

#. To visualize the current time in the video, we will add a time value in the corner of the viewer. This can be added using Sources > Annotation > Annotate Time. Press ``Apply`` to apply the time annotation source.
#. Select the JOREKReader and apply a ``Temporal Interpolator`` filter. This can be found under Filters > Temporal > Temporal Interpolator.
#. Set the ``Discrete Time Step Interval`` to 0.01, and select ``Apply`` to apply the temporal interpolation.
#. Repeat the previous two steps for the GGDReader containing the Wall Currents.
#. The clip must be applied to the Temporal Interpolator filter instead, so we can right-click the clip filter, select ``Copy Pipeline``, select the temporal interpolator and select ``Paste Pipeline``. Ensure to select the wall currents again in the coloring section.
#. Verify that both temporal interpolators are working, this can be done by opening View > Time Manager and checking if the two temporal interpolators have the same number of time steps, and that the time steps are of equal size.
#. We will now create an animation of the JOREK electron temperature and the wall currents over time. Place the objects in the viewpoints in the orientation that you want for the video. To create a video, go to File > Save Animation, provide a directory and a name for the video and select ``OK``.
#. In the pop-up window, the video settings, such as image resolution, compression and more can be changed. You can change these to whichever you like, but in this example we will only change the frame rate to 24. Press ``OK`` to start generating the animation, this may take a while.

Machine Description case
------------------------

TODO

JINTRAC case
------------

TODO
