"""IMASPy version of the paraview plugin classes.
"""

import logging

from paraview.util.vtkAlgorithm import smhint, smproxy

from vtkggdtools.plugins.base_class import GGDVTKPluginBase

logger = logging.getLogger("vtkggdtools")


@smproxy.source(label="IMASPy GGDReader")
@smhint.xml("""<ShowInMenu category="VTKGGDTools" />""")
class IMASPyGGDReader(GGDVTKPluginBase):
    pass
