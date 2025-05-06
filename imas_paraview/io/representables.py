from collections import OrderedDict, defaultdict
from dataclasses import dataclass, field

from vtkmodules.vtkCommonDataModel import vtkUnstructuredGrid


@dataclass
class GridGGDRepresentable:
    """
    A representation of the index transforms required to go from point and cell IDS
    in the VTK unstructured grid to the object indices in a grid_ggd IDS node.
    The indices are stored in C/python index convention (0, ... n)
    Prior to assigning these indices to a grid_ggd/.../attribute, increment by 1.
    """

    ugrid: vtkUnstructuredGrid = None
    edges: OrderedDict = field(
        default_factory=OrderedDict
    )  #: Go from edge (p1,p2) to 1d grid_ggd object index.
    faces: OrderedDict = field(
        default_factory=OrderedDict
    )  #: Go from face (p1,..) to 2d grid_ggd object index.
    volumes: OrderedDict = field(
        default_factory=OrderedDict
    )  #: Go from a 3D cell's id to 3d grid_ggd object index.


@dataclass
class GridSubsetRepresentable:
    ugrid: vtkUnstructuredGrid = None
    num_subsets: int = 0
    element_list: defaultdict = field(default_factory=lambda: defaultdict(list))
    subset_cell_ids: defaultdict = field(default_factory=lambda: defaultdict(list))
    subset_cell_types: defaultdict = field(default_factory=lambda: defaultdict(list))
