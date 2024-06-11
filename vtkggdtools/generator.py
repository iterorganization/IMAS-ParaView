import os
import sys
import xml.etree.ElementTree as ET
from collections import defaultdict
from dataclasses import dataclass, field
from datetime import datetime
from xml.etree.ElementTree import Element

from identifiers.ggd_identifier import ggd_identifier
from identifiers.ggd_space_identifier import ggd_space_identifier
from jinja2 import Environment, FileSystemLoader, select_autoescape

dd_version = os.getenv("IMAS_VERSION")
ual_version = os.getenv("UAL_VERSION")
imas_root = os.getenv("IMAS_PREFIX")
idsdefxml = os.path.join(imas_root, "include/IDSDef.xml")

base_dir = os.path.dirname(__file__)


@dataclass
class IDSTemplate:
    name: str = ""
    label: str = ""
    description: str = ""
    grid_ggd_path: list = field(default_factory=list)
    aos_indices: dict = field(default_factory=dict)


@dataclass
class IMASDefaults:
    backend: int = 12
    imas_backends: dict = field(default_factory=dict)
    shot: int = 122408
    run: int = 3
    username: str = "public"
    tokamak: str = "ITER"
    major_ver: str = "3"
    occurrence: int = 0
    delimiter: str = "/"
    ggd_type: int = 4
    ggd_types: dict = field(default_factory=dict)
    ggd_space_type: int = 1
    ggd_space_types: dict = field(default_factory=dict)
    n_plane: int = 0  # number of planes for 3D interpolation
    phi_start: float = 0.0  # start phi coordinate in degrees
    phi_end: float = 90.0  # end phi coordinate in degrees

    def __post_init__(self):
        self.ggd_types.update(ggd_identifier)
        self.ggd_space_types.update(ggd_space_identifier)
        self.imas_backends.update({"ASCII": 11, "MDSPLUS": 12, "HDF5": 13, "UDA": 15})


def delete_obsoleted(nodes: list):
    for node in nodes.copy():
        if node.attrib.get("lifecycle_status") == "obsolescent":
            nodes.remove(node)


def to_camel_case(inp_string: str, delimiter: str = "_"):
    return "".join([text_item.capitalize() for text_item in inp_string.split("_")])


def get_ancestor_with_tag(node: Element, lut: dict, tag: str = ""):
    child = node
    parent = lut[child]

    while parent.tag != tag:
        child = parent
        parent = lut[child]
    return parent


def get_ids_ancestor(node: Element, lut: dict):
    return get_ancestor_with_tag(node, lut, "IDS")


def get_ancestors_upto(node: Element, target: Element, lut: dict):
    child = node
    parent = lut[child]
    result = []

    while parent != target:
        result.append(parent)
        child = parent
        parent = lut[child]
    if len(result):
        result.reverse()
    return result


def is_aos_node(node: Element):
    return node.attrib.get("data_type") == "struct_array"


def get_aos_var_name(node: Element):
    is_dynamic = node.attrib.get("type") == "dynamic"

    if not is_aos_node(node):
        return ""

    if is_dynamic and "itime" in node.attrib["path_doc"]:
        return "TimeIdx"
    else:
        name = node.attrib["name"]
        return to_camel_case(name) + "Idx"


def prepare_blueprints(tree: ET):
    parent_map = {c: p for p in tree.iter() for c in p}

    structure_references = ["generic_grid_dynamic", "generic_grid_aos3_root"]
    paths = defaultdict(set)
    blueprints = dict()

    for struct_ref in structure_references:

        nodes = tree.getroot().findall(
            f".//field/..[@structure_reference='{struct_ref}']"
        )
        delete_obsoleted(nodes)

        for node in nodes:
            try:
                ids_ancestor = get_ids_ancestor(node, parent_map)
                other_ancestors = get_ancestors_upto(node, ids_ancestor, parent_map)
            except KeyError:
                print(
                    f"Warning: can't find ancestor of {node.attrib['path']}."
                    " Continuing."
                )
                continue
            ids_name = ids_ancestor.attrib["name"]

            bp = blueprints.get(ids_name)
            if not bp:
                first_time_init = True
            else:
                first_time_init = False

            if first_time_init:
                bp = IDSTemplate()
                blueprints.update({ids_name: bp})
                # used in templating.
                bp.name = ids_name
                bp.description = ids_ancestor.attrib["documentation"]
                # Convert to vtk class name convention. edge_profiles <-> EdgeProfiles
                bp.label = to_camel_case(ids_name)

            # full_path = ids_name + IMASDefaults.delimiter + node.attrib["path"]
            full_path_doc = ids_name + IMASDefaults.delimiter + node.attrib["path_doc"]
            paths[ids_name].add(full_path_doc)

            for _node in [*other_ancestors, node]:

                name = _node.attrib["name"]

                if first_time_init:
                    bp.grid_ggd_path.append(name)

                if is_aos_node(_node):
                    bp.aos_indices.update({name: get_aos_var_name(_node)})

    for path in paths.items():
        print(path)
    return blueprints


def generate(input_dir: str, output_dir: str):
    env = Environment(
        loader=FileSystemLoader(input_dir),
        autoescape=select_autoescape(),
        trim_blocks=True,
        lstrip_blocks=True,
    )

    tree = ET.parse(idsdefxml)
    blueprints = prepare_blueprints(tree)
    template = env.get_template("VTKGGDTools.py.jinja")
    template.globals["now"] = datetime.now
    output = template.render(
        ids_blueprints=blueprints,
        default=IMASDefaults(),
        dd_ver=dd_version,
        ual_ver=ual_version,
        warning_msg="WARNING: Please refrain from editing this file.",
    )

    with open(os.path.join(output_dir, "VTKGGDTools.py"), "w+") as f:
        f.write(output)


if __name__ == "__main__":

    if len(sys.argv) >= 2:
        generate(sys.argv[1], sys.argv[2])
    else:
        generate(
            os.path.join(base_dir, "plugins/templates"),
            os.path.join(base_dir, "plugins"),
        )
