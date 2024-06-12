"""IMASPy version of the paraview plugin classes.
"""

import imaspy
from paraview import logger
from paraview.util.vtkAlgorithm import smdomain, smhint, smproperty, smproxy
from vtkmodules.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtkmodules.vtkCommonCore import vtkPoints
from vtkmodules.vtkCommonDataModel import (
    vtkCompositeDataSet,
    vtkDataAssembly,
    vtkDataObject,
    vtkPartitionedDataSetCollection,
)

from vtkggdtools.io import read_bezier, read_geom, read_ps
from vtkggdtools.util import get_first_grid


@smproxy.source(label="IMASPy GGDReader")
@smhint.xml("""<ShowInMenu category="VTKGGDTools" />""")
class IMASPyGGDReader(VTKPythonAlgorithmBase):
    """GGD Reader based on IMASPy"""

    def __init__(self):
        super().__init__(
            nInputPorts=0,
            nOutputPorts=1,
            outputType="vtkPartitionedDataSetCollection",
        )
        # Properties
        self._uri = ""
        self._idsname = ""
        self._occurrence = 0
        self._n_plane = 0
        self._phi_start = 0
        self._phi_end = 0

        # Data caches
        self._dbentry = None
        self._ids = None

    @smproperty.stringvector(
        name="IMAS URI",
        default_values="imas:hdf5?path=/home/maarten/projects/iter-python/ggd-vtk/data/ggd-vtk-testdb/",
    )
    def SetURI(self, uri):
        if uri != self._uri:
            self._uri = uri
            if self._dbentry is not None:
                self._dbentry.close()
                self._ids = self._dbentry = None
            self.Modified()

    @smproperty.stringvector(name="IDS", default_values="mhd")
    def SetIDS(self, idsname):
        if idsname != self._idsname:
            self._idsname = idsname
            self._ids = None
            self.Modified()

    @smproperty.intvector(name="Occurrence", default_values=0)
    def SetOccurrence(self, occurrence):
        if occurrence != self._occurrence:
            self._occurrence = occurrence
            self._ids = None
            self.Modified()

    @smproperty.intvector(name="N plane", default_values=0)
    def SetNPlane(self, val):
        if val != self._n_plane:
            self._n_plane = val
            self.Modified()

    @smproperty.doublevector(name="Phi range", default_values=[0, 0])
    @smdomain.doublerange(min=0, max=360.0)
    def SetPhiRange(self, val, val2):
        if val != self._phi_start or val2 != self._phi_end:
            self._phi_start = val
            self._phi_end = val2
            self.Modified()

    def FillOutputPortInformation(self, port, info):
        info.Set(vtkDataObject.DATA_TYPE_NAME(), "vtkPartitionedDataSetCollection")
        return 1

    def RequestDataObject(self, request, inInfo, outInfo):
        output = vtkPartitionedDataSetCollection()
        outInfo.GetInformationObject(0).Set(vtkDataObject.DATA_OBJECT(), output)
        outInfo.GetInformationObject(0).Set(
            vtkDataObject.DATA_EXTENT_TYPE(), output.GetExtentType()
        )
        return 1

    def _ensure_ids(self):
        if self._ids is None:
            if self._dbentry is None:
                self._dbentry = imaspy.DBEntry(self._uri, "r")

            # TODO: enable lazy loading
            logger.info("Loading IDS %s/%d ...", self._idsname, self._occurrence)
            self._ids = self._dbentry.get(
                self._idsname,
                self._occurrence,
                autoconvert=False,
            )
            logger.info("Done loading IDS.")

    def RequestData(self, request, inInfo, outInfo):
        self._ensure_ids()

        # TODO: allow selecting other grids
        _aos_index_values = FauxIndexMap()
        grid_ggd = get_first_grid(self._ids)

        # We now have the grid_ggd
        # Check if we have anything to read:
        if len(grid_ggd.grid_subset) < 1 and len(grid_ggd.space) < 1:
            logger.info("No points to read from grid_ggd")
            return 0

        output = vtkPartitionedDataSetCollection.GetData(outInfo)
        num_subsets = len(grid_ggd.grid_subset)
        points = vtkPoints()
        space_idx = 0

        read_geom.fill_vtk_points(grid_ggd, space_idx, points, self._idsname)
        assembly = vtkDataAssembly()
        output.SetDataAssembly(assembly)

        # Interpolate JOREK Fourier space
        # TODO: figure out if we can put this functionality in a post-processing step?
        if self._n_plane != 0:
            number_of_spaces = len(grid_ggd.space)
            if number_of_spaces > 1 and len(grid_ggd.space[0].coordinates_type) == 2:
                n_period = grid_ggd.space[1].geometry_type.index
                if n_period > 0:  # Fourier space with periodicity (JOREK)
                    logger.info(
                        f"Reading Bezier mesh with Fourier periodiciy {n_period}"
                    )
                    ugrid = read_bezier.convert_grid_subset_to_unstructured_grid(
                        self._idsname,
                        self._ids,
                        _aos_index_values,
                        self._n_plane,
                        self._phi_start,
                        self._phi_end,
                    )
                    output.SetPartition(0, 0, ugrid)
                    child = assembly.AddNode(self._idsname, 0)
                    assembly.AddDataSetIndex(child, 0)
                    output.GetMetaData(0).Set(vtkCompositeDataSet.NAME(), self._idsname)
                else:
                    logger.error(
                        f"Number of planes {self._n_plane} invalid for this "
                        f"{number_of_spaces} number of spaces"
                    )
            else:
                logger.error(
                    f"Number of planes {self._n_plane} invalid for this IDS type."
                    " Try using N = 0"
                )
            return 1

        def fill_grid_and_plasma_state(subset_idx, partition):
            subset = None if subset_idx < 0 else grid_ggd.grid_subset[subset_idx]
            ugrid = read_geom.convert_grid_subset_geometry_to_unstructured_grid(
                grid_ggd, subset_idx, points
            )
            read_ps.read_plasma_state(
                self._idsname, self._ids, _aos_index_values, subset_idx, ugrid
            )
            output.SetPartition(partition, 0, ugrid)
            label = str(subset.identifier.name) if subset else self._idsname
            child = assembly.AddNode(label.replace(" ", "_"), 0)
            assembly.AddDataSetIndex(child, partition)
            output.GetMetaData(partition).Set(vtkCompositeDataSet.NAME(), label)

        # Regular grid reading
        if num_subsets <= 1:
            logger.info("No subsets to read from grid_ggd")
            output.SetNumberOfPartitionedDataSets(1)
            fill_grid_and_plasma_state(-1, 0)
        elif self._idsname == "wall":
            # FIXME: what if num_subsets is 2 or 3?
            output.SetNumberOfPartitionedDataSets(num_subsets - 3)
            fill_grid_and_plasma_state(-1, 0)
            for subset_idx in range(4, num_subsets):
                fill_grid_and_plasma_state(subset_idx, subset_idx - 3)
                self.UpdateProgress(self.GetProgress() + 1 / num_subsets)
        else:
            output.SetNumberOfPartitionedDataSets(num_subsets)
            for subset_idx in range(num_subsets):
                fill_grid_and_plasma_state(subset_idx, subset_idx)
                self.UpdateProgress(self.GetProgress() + 1 / num_subsets)

        return 1


class FauxIndexMap:
    def __getitem__(self, item):
        return 0

    def get(self, name, default=None):
        return 0
