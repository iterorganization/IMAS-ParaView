IMAS-ParaView development guide
===============================

This page is for developers that want to improve or extend IMAS-ParaView.

Project folder structure
------------------------

- ``.github`` contains GitHub specific files, such as the GitHub Actions configuration
- ``benchmarks`` contains performance benchmarks for the conversion logic. We use the
  Airspeed Velocity (``asv``) as benchmark framework. More information on this framework
  can be found on their website: <https://asv.readthedocs.io/en/stable/index.html>
- ``data`` contains datasets that are used in the CI tests.
- ``docs`` contains all source files for the documentation (which you are currently
  reading). The documentation is based on `Sphinx
  <https://www.sphinx-doc.org/en/master/index.html>`__ and uses the `reStructuredText`
  (.rst) format. More information on this format can be found on the Sphinx website.
- ``imas_paraview`` contains all source code for the plugins and conversion tool:

  - ``io`` contains modules to read GGD data from IDSs into VTK datasets
  - ``paraview_support`` contains tools to more easily work with the ParaView plugin
    library.
  - ``plugins`` contains all ParaView plugins
  - ``tests`` contains all unit and integration tests

More details on the Python methods can be found in the :ref:`api reference`, or as
docstrings in the files.


High-level overview of the plugins
----------------------------------

Several plugins exist in this package. All source plugins have logic in common, which is
why they all derive from the same base class. The class diagram, at the time of writing,
is as follows:

- :py:class:`~imas_paraview.plugins.base_class.GGDVTKPluginBase`. This is the base
  class, which defines all 

  - :py:class:`~imas_paraview.plugins.ggd_base_reader.GGDBaseReader`. This is the base
    class for plugins that read GGD data:

    - :py:class:`~imas_paraview.plugins.vtkggdreader.GGDReader`. This is the GGD
      Reader plugin.
    - :py:class:`~imas_paraview.plugins.jorek.JorekGGDReader`. This is the JOREK Reader
      plugin.

  - :py:class:`~imas_paraview.plugins.beam.BeamReader`. This is the Beam Reader.
  - :py:class:`~imas_paraview.plugins.line_of_sight.LineOfSightReader`. This is the Line
    of Sight Reader.
  - :py:class:`~imas_paraview.plugins.position.PositionReader`. This is the Position
    Reader.
  - :py:class:`~imas_paraview.plugins.profiles_1d.Profiles1DReader`. This is the 1D
    Profiles Reader.
  - :py:class:`~imas_paraview.plugins.profiles_2d.Profiles2DReader`. This is the 2D
    Profiles Reader.
  - :py:class:`~imas_paraview.plugins.wall_limiter.WallLimiterReader`. This is the Wall
    Limiter Reader.
- :py:class:`~imas_paraview.plugins.profiles_1d_mapper.Profiles1DMapper`. This is the
  1D Profiles Mapper filter.


The plugin base class
---------------------

As mentioned in the previous section,
:py:class:`~imas_paraview.plugins.base_class.GGDVTKPluginBase` is the base class for all
source plugins. The base class is responsible for all ParaView properties (the UI that
can be adjusted in the ParaView Properties Panel):

- URI selection.
- IDS selection. Each implementing class indicates which IDSs are supported (through the
  ``supported_ids`` argument of the
  :py:meth:`~imas_paraview.plugins.base_class.GGDVTKPluginBase.__init__` method. The
  base class will only advertise supported IDSs that exist in the selected URI.
- Sub index selection, called ``ParentIndex`` in the code. This is only enabled in the
  GGD Reader (which sets ``_show_parent_indices_dropdown`` to ``True``).
- The attribute selection (Select attribute Arrays selector in the UI). Subclasses
  fill the ``_selectable`` attribute after loading the IDS to tell the base class which
  items to display in this list.
- Bézier interpolation settings. These settings are only used by the JOREK Reader.

Note that all Properties / UI elements are declard on the base class. The reason is that
this is the only way to keep their ordering consistent in ParaView. The ParaView Python
Plugin API does some interesting things with subclasses...

The base class also handles timing logic. Plugins that work with time-dependent data
need to indicate this when subclassing, see for example `the Profiles1DReader
<https://github.com/iterorganization/IMAS-ParaView/blob/develop/imas_paraview/plugins/profiles_1d.py#L22>`__.


The base reader has some methods which can be implemented by subclasses to
provide the actual functionality of the plugin:

- :py:meth:`~imas_paraview.plugins.base_class.GGDVTKPluginBase.get_heterogeneous_time_array`
  if the reader wants to support IDSs using heterogeneous time.
- :py:meth:`~imas_paraview.plugins.base_class.GGDVTKPluginBase.request_information` for
  any logic that needs to execute during the RequestInformation stage in ParaView.
- :py:meth:`~imas_paraview.plugins.base_class.GGDVTKPluginBase.setup_ids` is called
  after loading a (new) IDS. This is mostly used for indicating which data is available
  in the IDS.


How to extend the list of supported IDSs in a plugin
----------------------------------------------------

Most plugins can handle multiple input IDSs and will use the Data Dictionary definitions
provided by `IMAS-Python <https://github.com/iterorganization/IMAS-Python>`__ to
dynamically search for data that can be displayed.

For example, the GGD Reader has a list of ``SUPPORTED_IDS_NAMES`` at the top of the
Python file. New IDSs can be added to this list to allow the GGD Reader to read those as
well.
