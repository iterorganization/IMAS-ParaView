import logging

import numpy as np
from vtkmodules.numpy_interface import dataset_adapter as dsa
from vtkmodules.vtkCommonCore import vtkDoubleArray
from vtkmodules.vtkCommonDataModel import vtkCellData, vtkPointData, vtkUnstructuredGrid

# We'll need these below when we create some units manually:
import vtkggdtools.imashelper

# For the units:
from vtkggdtools.imashelper import get_units

logger = logging.getLogger(__name__)

u_pre = vtkggdtools.imashelper.u_pre
u_post = vtkggdtools.imashelper.u_post


SUPPORTED_IDS_NAMES = [
    "distribution_sources",
    "distributions",
    "edge_profiles",
    "edge_sources",
    "edge_transport",
    "equilibrium",
    "mhd",
    "radiation",
    "tf",
    "transport_solver_numerics",
    "wall",
    "waves",
]


def read_plasma_state(
    ids_name: str,
    ids_obj,
    aos_index_values: dict,
    subset_idx: int,
    ugrid: vtkUnstructuredGrid,
) -> None:
    """
    Reads plasma state data arrays from the ggd node. These arrays are added as point
    data or cell data to the unstructured grid.
    This is just a dispatching routine, the real work is done in the called
    funtions read_<IDS_name>(ids_obj, aos_index_values, subset_idx, ugrid)
    :param ids_name: name of the top level IDS.
    :param ids_obj: an ids_obj of type ids_name
    :param aos_index_values: the values that shall be used to navigate the AoS and reach
        the scalar arrays.
    :param subset_idx: an index into grid_ggd/grid_subset AoS
    :param ugrid: the unstructured grid instance.
    :return: None
    """
    if ids_name == "distribution_sources":
        read_distribution_sources(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == "distributions":
        read_distributions(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == "edge_profiles":
        read_edge_profiles(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == "edge_sources":
        read_edge_sources(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == "edge_transport":
        read_edge_transport(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == "equilibrium":
        read_equilibrium(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == "mhd":
        read_mhd(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == "radiation":
        read_radiation(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == "tf":
        read_tf(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == "transport_solver_numerics":
        read_transport_solver_numerics(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == "wall":
        read_wall(ids_obj, aos_index_values, subset_idx, ugrid)

    elif ids_name == "waves":
        read_waves(ids_obj, aos_index_values, subset_idx, ugrid)

    else:
        logger.warn(f"Reading plasma state from IDS {ids_name} not implemented.")


def read_distribution_sources(
    ids_obj, aos_index_values: dict, subset_idx: int, ugrid: vtkUnstructuredGrid
) -> None:

    # ggd is at /distribution_sources/source(i1)/ggd(itime)
    source_idx = aos_index_values.get("SourceIdx")
    time_idx = aos_index_values.get("TimeIdx")
    try:
        ggd = ids_obj.source[source_idx].ggd[time_idx]
    except IndexError:
        return

    name = "Particle Density"
    _add_aos_scalar_array_to_vtk_field_data(ggd.particles, subset_idx, name, ugrid)


def read_distributions(
    ids_obj, aos_index_values: dict, subset_idx: int, ugrid: vtkUnstructuredGrid
) -> None:

    # ggd is at /distributions/distribution(i1)/ggd(itime)
    distribution_idx = aos_index_values.get("DistributionIdx")
    time_idx = aos_index_values.get("TimeIdx")
    try:
        ggd = ids_obj.distribution[distribution_idx].ggd[time_idx]
    except IndexError:
        return

    # - expansion(i2)
    num_expansions = len(ggd.expansion)
    expansion_components = []

    for i in range(num_expansions):
        expansion = ggd.expansion[i].grid_subset
        if subset_idx >= len(expansion):
            continue
        expansion_components.append(expansion[subset_idx].values)
        name = f"Distribution Function Order: {i}"
        _add_aos_scalar_array_to_vtk_field_data(expansion, subset_idx, name, ugrid)

    if len(expansion_components):
        total_df = np.sum(expansion_components, axis=0)
        name = "Distribution Function"
        _add_scalar_array_to_vtk_field_data(total_df, name, ugrid)


def read_edge_profiles(
    ids_obj, aos_index_values: dict, subset_idx: int, ugrid: vtkUnstructuredGrid
) -> None:

    # ggd is at /edge_profiles/ggd(itime)
    time_idx = aos_index_values.get("TimeIdx")
    try:
        ggd = ids_obj.ggd[time_idx]
    except IndexError:
        return

    # for finding units:
    ids_name = "edge_profiles"
    ggd_path = "ggd"

    # electrons
    units_path = ggd_path + "/electrons"
    # scalar quantities, each value is a tuple with ('description', units?):
    quantities = {
        # scalars
        "temperature": ("Electron Temperature", True),
        "density": ("Electron Density", True),
        "density_fast": ("Electron Density (fast)", True),
        "pressure": ("Electron Pressure", True),
        "pressure_fast_perpendicular": ("Electron Pressure (fast perpendicular)", True),
        "pressure_fast_parallel": ("Electron Pressure (fast parallel)", True),
        "distribution_function": ("Electron Distribution Function", False),  # no units
    }
    for q_name in quantities:
        (name, use_units) = quantities[q_name]
        if use_units:
            units = get_units(ids_name, f"{units_path}/{q_name}")
            name += f" {units}"
        _add_aos_scalar_array_to_vtk_field_data(
            getattr(ggd.electrons, q_name), subset_idx, name, ugrid
        )
    # vector array:
    units = get_units(ids_name, units_path + "/velocity")
    name = f"Electron Velocity {units}"
    _add_aos_vector_array_to_vtk_field_data(
        ggd.electrons.velocity, subset_idx, name, ugrid
    )

    # ions, neutrals = heavy particles = hp
    # Add ' (all ions)' to the ion name to distinguish from state below:
    ion_names = [" (all ions)"] * len(ggd.ion) + [""] * len(ggd.neutral)
    for hp, i_name in list(zip([*ggd.ion, *ggd.neutral], ion_names)):
        hp_name = hp.label + i_name  # ' add (all ions)' to ions or '' to neutrals
        if i_name == "":
            units_hp_path = ggd_path + "/neutral"
        else:
            units_hp_path = ggd_path + "/ion"
        # scalar arrays:
        quantities = {
            "temperature": f"{hp_name} Temperature",
            "density": f"{hp_name} Density",
            "density_fast": f"{hp_name} Density (fast)",
            "pressure": f"{hp_name} Pressure",
            "pressure_fast_perpendicular": f"{hp_name} Pressure (fast perpendicular)",
            "pressure_fast_parallel": f"{hp_name} Pressure (fast parallel)",
            "energy_density_kinetic": f"{hp_name} Kinetic Energy Density",
        }
        for q_name in quantities:
            units = get_units(ids_name, f"{units_hp_path}/{q_name}")
            _add_aos_scalar_array_to_vtk_field_data(
                getattr(hp, q_name), subset_idx, f"{quantities[q_name]} {units}", ugrid
            )
        # vector array:
        units = get_units(ids_name, units_hp_path + "/velocity")
        name = f"{hp_name} Velocity {units}"
        _add_aos_vector_array_to_vtk_field_data(hp.velocity, subset_idx, name, ugrid)
        # heavy particles: state
        units_hp_path += "/state"
        for state in hp.state:
            # scalars, each value is a tuple with ('description', units?):
            quantities = {
                "temperature": (f"{state.label} Temperature", True),
                "density": (f"{state.label} Density", True),
                "density_fast": (f"{state.label} Density (fast)", True),
                "pressure": (f"{state.label} Pressure", True),
                "pressure_fast_perpendicular": (
                    f"{state.label} Pressure (fast perpendicular)",
                    True,
                ),
                "pressure_fast_parallel": (
                    f"{state.label} Pressure (fast parallel)",
                    True,
                ),
                "energy_density_kinetic": (
                    f"{state.label} Kinetic Energy Density",
                    True,
                ),
                "distribution_function": (
                    f"{state.label} Distribution Function",
                    False,
                ),  # no units
            }
            if (
                i_name
            ):  # only for ions, units manually written here as 'e' (e- unit charge)
                quantities["z_average"] = (f"{state.label} <Z> {u_pre}e{u_post}", False)
                quantities["z_square_average"] = (
                    f"{state.label} <Z^2> {u_pre}e{u_post}",
                    False,
                )
                quantities["ionisation_potential"] = (
                    f"{state.label} Ionisation Potential {u_pre}e{u_post}",
                    False,
                )
            for q_name in quantities:
                (name, use_units) = quantities[q_name]
                if use_units:
                    units = get_units(ids_name, f"{units_hp_path}/{q_name}")
                    name += f" {units}"
                _add_aos_scalar_array_to_vtk_field_data(
                    getattr(state, q_name), subset_idx, name, ugrid
                )
            # vectors:
            quantities = {
                "velocity": f"{state.label} Velocity",
                "velocity_diamagnetic": f"{state.label} Velocity (diamagnetic)",
                "velocity_exb": f"{state.label} Velocity (ExB)",
            }
            for q_name in quantities:
                units = get_units(ids_name, f"{units_hp_path}/{q_name}")
                _add_aos_vector_array_to_vtk_field_data(
                    getattr(state, q_name),
                    subset_idx,
                    f"{quantities[q_name]} {units}",
                    ugrid,
                )

    # Other scalar quantities:
    units_path = ggd_path
    quantities = {
        "t_i_average": "Average t_i",
        "n_i_total_over_n_e": "Total n_i over n_e",
        "zeff": "Zeff",
        "pressure_thermal": "Pressure (thermal)",
        "pressure_perpendicular": "Pressure (perpendicular)",
        "pressure_parallel": "Pressure (parallel)",
        "j_parallel": "Current (parallel)",
        "phi_potential": "Potential Phi",
    }
    for q_name in quantities:
        units = get_units(ids_name, f"{units_path}/{q_name}")
        _add_aos_scalar_array_to_vtk_field_data(
            getattr(ggd, q_name), subset_idx, f"{quantities[q_name]} {units}", ugrid
        )

    # other vector quantities
    quantities = {
        "j_total": "Current (total)",
        "j_anomalous": "Current (anomalous)",
        "j_inertial": "Current (inertial)",
        "j_ion_neutral_friction": "Current (ion+neutral friction)",
        "j_parallel_viscosity": "Current (parallel viscosity)",
        "j_perpendicular_viscosity": "Current (perpendicular viscosity)",
        "j_heat_viscosity": "Current (heat viscosity)",
        "j_pfirsch_schlueter": "Current (Pfirsch-Schlueter effects)",
        "j_diamagnetic": "Current (diamagnetic)",
        "e_field": "E",
    }
    for q_name in quantities:
        units = get_units(ids_name, f"{units_path}/{q_name}")
        _add_aos_vector_array_to_vtk_field_data(
            getattr(ggd, q_name), subset_idx, f"{quantities[q_name]} {units}", ugrid
        )

    # TODO: ggd_fast...


def read_edge_profiles_lecad(
    ids_obj, aos_index_values: dict, subset_idx: int, ugrid: vtkUnstructuredGrid
) -> None:
    # ggd is at /edge_profiles/ggd(itime)
    time_idx = aos_index_values.get("TimeIdx")
    try:
        ggd = ids_obj.ggd[time_idx]
    except IndexError:
        return

    # for finding units:
    ids_name = "edge_profiles"
    ggd_path = "ggd"

    # quantities for electrons, ions and neutrals
    # tuple shape: (object discription, if values are scalars or not)
    # velocity values are vectors, other are scalars
    quantities = {
        "density": ("Density", True),
        "temperature": ("Temperature", True),
        "density_fast": ("Density Fast", True),
        "pressure": ("Pressure", True),
        "pressure_fast_perpendicular": ("Pressure Fast Perpendicular", True),
        "pressure_fast_parallel": ("Pressure Fast Parallel", True),
        "velocity": ("Velocity", False),
    }
    for q_name in quantities:
        (name, scalar) = quantities[q_name]
        # electrons
        units_path = ggd_path + "/electrons"
        units = get_units(ids_name, f"{units_path}/{q_name}")
        if scalar:
            _add_aos_scalar_array_to_vtk_field_data(
                getattr(ggd.electrons, q_name),
                subset_idx,
                f"Electron {name} {units}",
                ugrid,
            )
        else:
            _add_aos_vector_array_to_vtk_field_data(
                getattr(ggd.electrons, q_name),
                subset_idx,
                "Electron {name} {units}",
                ugrid,
            )

        # ions
        for j in range(len(ggd.ion)):
            # find ion label
            if ggd.ion[j].label != "":
                ionLabel = ggd.ion[j].label
            else:
                ionLabel = ggd.ion[j].state[0].label
            units_path = ggd_path + "/ion"
            for k in range(len(getattr(ggd.ion[j], q_name))):
                # Check if there is actually an array
                if len(getattr(ggd.ion[j], q_name)[k].values) < 1:
                    break
                units = get_units(ids_name, f"{units_path}/{q_name}")
                name = f"Ion {ionLabel} {q_name} {units}"
                if scalar:
                    _add_aos_scalar_array_to_vtk_field_data(
                        getattr(ggd.ion[j], q_name), subset_idx, name, ugrid
                    )
                else:
                    _add_aos_vector_array_to_vtk_field_data(
                        getattr(ggd.ion[j], q_name), subset_idx, name, ugrid
                    )
        # neutrals
        for j in range(len(ggd.neutral)):
            # find neutral label
            if ggd.neutral[j].label != "":
                neutralLabel = ggd.neutral[j].label
            else:
                neutralLabel = ggd.neutral[j].state[0].label
            units_path = ggd_path + "/neutral"
            try:
                for k in range(len(getattr(ggd.neutral[j], q_name))):
                    # Check if there is actually an array
                    if len(getattr(ggd.neutral[j], q_name)[k].values) < 1:
                        break
                    units = get_units(ids_name, f"{units_path}/{q_name}")
                    name = f"Neutral {neutralLabel} {q_name} {units}"

                    # FIXME: path was undefined, but this will definitely raise an error
                    # inside the getattr calls below!
                    path = None
                    if scalar:
                        _add_aos_scalar_array_to_vtk_field_data(
                            getattr(ggd.neutral[j], path), subset_idx, name, ugrid
                        )
                    else:
                        _add_aos_vector_array_to_vtk_field_data(
                            getattr(ggd.neutral[j], path), subset_idx, name, ugrid
                        )
            except Exception:
                break

    # TODO ggd_fast


def read_edge_sources(
    ids_obj, aos_index_values: dict, subset_idx: int, ugrid: vtkUnstructuredGrid
) -> None:

    # ggd is at /edge_sources/source/ggd(itime)
    time_idx = aos_index_values.get("TimeIdx")
    # model_idx = aos_index_values.get("ModelIdx")
    # For units reading:
    ids_name = "edge_sources"
    ggd_path = "source/ggd"

    for source in ids_obj.source:
        try:
            ggd = source.ggd[time_idx]
        except IndexError:
            return
        s_name = source.identifier.name

        # electrons
        units_path = ggd_path + "/electrons"
        units = get_units(ids_name, units_path + "/particles")
        name = f"Electron Density ({s_name}) {units}"
        _add_aos_scalar_array_to_vtk_field_data(
            ggd.electrons.particles, subset_idx, name, ugrid
        )
        units = get_units(ids_name, units_path + "/energy")
        name = f"Electron Energy ({s_name}) {units}"
        _add_aos_scalar_array_to_vtk_field_data(
            ggd.electrons.energy, subset_idx, name, ugrid
        )

        # ions, neutrals = heavy particles = hp
        # Add ' (all ions)' to the ion name to distinguish from state below:
        ion_names = [" (all ions)"] * len(ggd.ion) + [""] * len(ggd.neutral)
        for hp, i_name in list(zip([*ggd.ion, *ggd.neutral], ion_names)):
            hp_name = hp.label + i_name  # add ' (all ions)' to ions or '' to neutrals
            if i_name == "":
                units_hp_path = ggd_path + "/neutral"
            else:
                units_hp_path = ggd_path + "/ion"
            units = get_units(ids_name, units_hp_path + "/particles")
            name = f"{hp_name} Density ({s_name}) {units}"
            _add_aos_scalar_array_to_vtk_field_data(
                hp.particles, subset_idx, name, ugrid
            )
            units = get_units(ids_name, units_hp_path + "/energy")
            name = f"{hp_name} Energy ({s_name}) {units}"
            _add_aos_scalar_array_to_vtk_field_data(hp.energy, subset_idx, name, ugrid)
            units = get_units(ids_name, units_hp_path + "/momentum")
            name = f"{hp_name} Momentum ({s_name}) {units}"
            _add_aos_vector_array_to_vtk_field_data(
                hp.momentum, subset_idx, name, ugrid
            )
            # heavy particles: state
            for state in hp.state:
                # sometimes state.label is not filled:
                try:
                    st_name = state.label
                except Exception:
                    st_name = "? " + hp_name
                units_path = units_hp_path + "/state"
                units = get_units(ids_name, units_path + "/particles")
                name = f"{st_name} Density ({s_name}) {units}"
                _add_aos_scalar_array_to_vtk_field_data(
                    state.particles, subset_idx, name, ugrid
                )
                units = get_units(ids_name, units_path + "/energy")
                name = f"{st_name} Energy ({s_name}) {units}"
                _add_aos_scalar_array_to_vtk_field_data(
                    state.energy, subset_idx, name, ugrid
                )
                units = get_units(ids_name, units_path + "/momentum")
                name = f"{st_name} Momentum ({s_name}) {units}"
                _add_aos_vector_array_to_vtk_field_data(
                    state.momentum, subset_idx, name, ugrid
                )

        # total_ion_energy
        units = get_units(ids_name, ggd_path + "/total_ion_energy")
        name = f"Total Ion Energy ({s_name}) {units}"
        _add_aos_scalar_array_to_vtk_field_data(
            ggd.total_ion_energy, subset_idx, name, ugrid
        )

        # momentum
        units = get_units(ids_name, ggd_path + "/momentum")
        name = f"Momentum ({s_name}) {units}"
        _add_aos_vector_array_to_vtk_field_data(ggd.momentum, subset_idx, name, ugrid)

        # current
        units = get_units(ids_name, ggd_path + "/current")
        name = f"Current ({s_name}) {units}"
        _add_aos_scalar_array_to_vtk_field_data(ggd.current, subset_idx, name, ugrid)

        # TODO: ggd_fast


def read_edge_transport(
    ids_obj, aos_index_values: dict, subset_idx: int, ugrid: vtkUnstructuredGrid
) -> None:

    # ggd is at /edge_transport/model(i1)/ggd(itime)
    # TODO: ModelIdx is not detected by generator.py and hence it is not present.
    # model_idx = aos_index_values.get('ModelIdx') or 0
    # time_idx = aos_index_values.get('TimeIdx')
    # try:
    #     ggd = ids_obj.model[model_idx].ggd[time_idx]
    # except IndexError:
    #     return

    time_idx = aos_index_values.get("TimeIdx")
    ids_name = "edge_transport"
    ggd_path = "model/ggd"
    for model in ids_obj.model:
        try:
            ggd = model.ggd[time_idx]
        except IndexError:
            continue

        # Use model name but only if we have more than one model:
        if len(ids_obj.model) > 1:
            m_name = f"({model.identifier.name})"
        else:
            m_name = ""

        # conductivity
        units_path = ggd_path + "/conductivity"
        units = get_units(ids_name, units_path)
        name = f"Conductivity {m_name} {units}"
        _add_aos_vector_array_to_vtk_field_data(
            ggd.conductivity, subset_idx, name, ugrid
        )

        # electrons
        electrons = ggd.electrons
        for q_name in ["particles", "energy"]:
            units_path = ggd_path + "/electrons/" + q_name
            quantity = getattr(electrons, q_name)
            units = get_units(ids_name, units_path + "/d")
            name = f"Electron Diffusivity ({q_name}) {units} {m_name}"
            _add_aos_scalar_array_to_vtk_field_data(quantity.d, subset_idx, name, ugrid)
            units = get_units(ids_name, units_path + "/v")
            name = f"Electron Convection ({q_name}) {units} {m_name}"
            _add_aos_scalar_array_to_vtk_field_data(quantity.v, subset_idx, name, ugrid)
            units = get_units(ids_name, units_path + "/flux")
            name = f"Electron Flux ({q_name}) {units} {m_name}"
            _add_aos_scalar_array_to_vtk_field_data(
                quantity.flux, subset_idx, name, ugrid
            )
            name = f"Electron Flux Limiter Coefficient ({q_name}) {m_name}"
            _add_aos_scalar_array_to_vtk_field_data(
                quantity.flux_limiter, subset_idx, name, ugrid
            )

        # total_ion_energy
        units_path = ggd_path + "/total_ion_energy"
        quantity = ggd.total_ion_energy
        q_name = "Total Ion Energy"
        units = get_units(ids_name, units_path + "/d")
        name = f"{q_name} Diffusivity {units} {m_name}"
        _add_aos_scalar_array_to_vtk_field_data(quantity.d, subset_idx, name, ugrid)
        units = get_units(ids_name, units_path + "/v")
        name = f"{q_name} Convection {units} {m_name}"
        _add_aos_scalar_array_to_vtk_field_data(quantity.v, subset_idx, name, ugrid)
        units = get_units(ids_name, units_path + "/flux")
        name = f"{q_name} Flux {units} {m_name}"
        _add_aos_scalar_array_to_vtk_field_data(quantity.flux, subset_idx, name, ugrid)
        name = f"{q_name} Flux Limiter Coefficient {m_name}"
        _add_aos_scalar_array_to_vtk_field_data(
            quantity.flux_limiter, subset_idx, name, ugrid
        )

        # momentum
        units_path = ggd_path + "/momentum"
        quantity = ggd.momentum
        q_name = "Momentum"
        units = get_units(ids_name, units_path + "/d")
        name = f"{q_name} Diffusivity {units} {m_name}"
        _add_aos_vector_array_to_vtk_field_data(quantity.d, subset_idx, name, ugrid)
        units = get_units(ids_name, units_path + "/v")
        name = f"{q_name} Convection {units} {m_name}"
        _add_aos_vector_array_to_vtk_field_data(quantity.v, subset_idx, name, ugrid)
        units = get_units(ids_name, units_path + "/flux")
        name = f"{q_name} Flux {units} {m_name}"
        _add_aos_vector_array_to_vtk_field_data(quantity.flux, subset_idx, name, ugrid)
        name = f"{q_name} Flux Limiter Coefficient {m_name}"
        _add_aos_vector_array_to_vtk_field_data(
            quantity.flux_limiter, subset_idx, name, ugrid
        )

        # ions, neutrals = heavy particles = hp
        # Add ' (all ions)' to the ion name to distinguish from state below:
        ion_names = [" (all ions)"] * len(ggd.ion) + [""] * len(ggd.neutral)
        for hp, i_name in list(zip([*ggd.ion, *ggd.neutral], ion_names)):
            hp_name = hp.label + i_name  # add ' (all ions)' to ions or '' to neutrals
            if i_name == "":
                units_hp_path = ggd_path + "/neutral"
            else:
                units_hp_path = ggd_path + "/ion"
            # heavy particles: particles, energy
            for q_name in ["particles", "energy"]:
                quantity = getattr(hp, q_name)
                units_path = units_hp_path + "/" + q_name
                units = get_units(ids_name, units_path + "/d")
                name = f"{hp_name} Diffusivity ({q_name}) {units} {m_name}"
                _add_aos_scalar_array_to_vtk_field_data(
                    quantity.d, subset_idx, name, ugrid
                )
                units = get_units(ids_name, units_path + "/v")
                name = f"{hp_name} Convection ({q_name}) {units} {m_name}"
                _add_aos_scalar_array_to_vtk_field_data(
                    quantity.v, subset_idx, name, ugrid
                )
                units = get_units(ids_name, units_path + "/flux")
                name = f"{hp_name} Flux ({q_name}) {units} {m_name}"
                _add_aos_scalar_array_to_vtk_field_data(
                    quantity.flux, subset_idx, name, ugrid
                )
                name = f"{hp_name} Flux Limiter Coefficient ({q_name}) {m_name}"
                _add_aos_scalar_array_to_vtk_field_data(
                    quantity.flux_limiter, subset_idx, name, ugrid
                )
            # heavy particles: momentum
            quantity = hp.momentum
            q_name = "Momentum"
            units_path = units_hp_path + "/momentum"
            units = get_units(ids_name, units_path + "/d")
            name = f"{hp_name} Diffusivity ({q_name}) {units} {m_name}"
            _add_aos_vector_array_to_vtk_field_data(quantity.d, subset_idx, name, ugrid)
            units = get_units(ids_name, units_path + "/v")
            name = f"{hp_name} Convection ({q_name}) {units} {m_name}"
            _add_aos_vector_array_to_vtk_field_data(quantity.v, subset_idx, name, ugrid)
            units = get_units(ids_name, units_path + "/flux")
            name = f"{hp_name} Flux ({q_name}) {units} {m_name}"
            _add_aos_vector_array_to_vtk_field_data(
                quantity.flux, subset_idx, name, ugrid
            )
            name = f"{hp_name} Flux Limiter Coefficient ({q_name}) {m_name}"
            _add_aos_vector_array_to_vtk_field_data(
                quantity.flux_limiter, subset_idx, name, ugrid
            )
            # heavy particles: state
            for state in hp.state:
                for q_name in ["particles", "energy"]:
                    units_path = units_hp_path + "/state/" + q_name
                    quantity = getattr(state, q_name)
                    units = get_units(ids_name, units_path + "/d")
                    name = f"{state.label} Diffusivity ({q_name}) {units} {m_name}"
                    _add_aos_scalar_array_to_vtk_field_data(
                        quantity.d, subset_idx, name, ugrid
                    )
                    units = get_units(ids_name, units_path + "/v")
                    name = f"{state.label} Convection ({q_name}) {units} {m_name}"
                    _add_aos_scalar_array_to_vtk_field_data(
                        quantity.v, subset_idx, name, ugrid
                    )
                    units = get_units(ids_name, units_path + "/flux")
                    name = f"{state.label} Flux ({q_name}) {units} {m_name}"
                    _add_aos_scalar_array_to_vtk_field_data(
                        quantity.flux, subset_idx, name, ugrid
                    )
                    name = f"{state.label} Flux Limiter Coefficient ({q_name}) {m_name}"
                    _add_aos_scalar_array_to_vtk_field_data(
                        quantity.flux_limiter, subset_idx, name, ugrid
                    )

                q_name = "momentum"
                units_path = units_hp_path + "/state/" + q_name
                quantity = state.momentum
                units = get_units(ids_name, units_path + "/d")
                name = f"{state.label} Diffusivity ({q_name}) {units} {m_name}"
                _add_aos_vector_array_to_vtk_field_data(
                    quantity.d, subset_idx, name, ugrid
                )
                units = get_units(ids_name, units_path + "/v")
                name = f"{state.label} Convection ({q_name}) {units} {m_name}"
                _add_aos_vector_array_to_vtk_field_data(
                    quantity.v, subset_idx, name, ugrid
                )
                units = get_units(ids_name, units_path + "/flux")
                name = f"{state.label} Flux ({q_name}) {units} {m_name}"
                _add_aos_vector_array_to_vtk_field_data(
                    quantity.flux, subset_idx, name, ugrid
                )
                name = f"{state.label} Flux Limiter Coefficient ({q_name}) {m_name}"
                _add_aos_vector_array_to_vtk_field_data(
                    quantity.flux_limiter, subset_idx, name, ugrid
                )
    # TODO: ggd_fast


def read_equilibrium(
    ids_obj, aos_index_values: dict, subset_idx: int, ugrid: vtkUnstructuredGrid
) -> None:

    # ggd is at /equilibrium/time_slice(itime)/ggd(i1)
    time_idx = aos_index_values.get("TimeIdx")
    ggd_idx = aos_index_values.get("GridIdx")
    try:
        ggd = ids_obj.time_slice[time_idx].ggd[ggd_idx]
    except IndexError:
        return

    # for finding units:
    ids_name = "equilibrium"
    ggd_path = "time_slice/ggd"
    units_path = ggd_path

    quantities = {
        "r": "Major Radius",
        "z": "Height",
        "psi": "Poloidal Flux",
        "phi": "Toroidal Flux",
        "theta": "Poloidal Angle",
        "j_tor": "Toroidal Plasma Current Density",
        "j_parallel": "Parallel Plasma Current Density",
        "b_field_r": "Magnetic Field Br",
        "b_field_z": "Magnetic Field Bz",
        "b_field_tor": "Magnetic Field Btor",
    }
    for q_name in quantities:
        units = get_units(ids_name, f"{units_path}/{q_name}")
        _add_aos_vector_array_to_vtk_field_data(
            getattr(ggd, q_name), subset_idx, f"{quantities[q_name]} {units}", ugrid
        )


def read_mhd(
    ids_obj, aos_index_values: dict, subset_idx: int, ugrid: vtkUnstructuredGrid
) -> None:

    # ggd is at /mhd/ggd(itime)
    time_idx = aos_index_values.get("TimeIdx")
    try:
        ggd = ids_obj.ggd[time_idx]
    except IndexError:
        return

    # for finding units:
    ids_name = "mhd"
    ggd_path = "ggd"
    units_path = ggd_path

    # - electrons.temperature
    q_name = "electrons/temperature"
    units = get_units(ids_name, f"{units_path}/{q_name}")
    name = f"Electron Temperature {units}"
    _add_aos_scalar_array_to_vtk_field_data(
        ggd.electrons.temperature, subset_idx, name, ugrid
    )
    # all the others...
    quantities = {
        "t_i_average": "Ion Temperature (average)",
        "n_i_total": "Ion Density (total)",
        "zeff": "Z effective",
        "b_field_r": "Magnetic Field Br",
        "b_field_z": "Magnetic Field Bz",
        "b_field_tor": "Magnetic Field Btor",
        "a_field_r": "Magnetic Potential Ar",
        "a_field_z": "Magnetic Potential Az",
        "a_field_tor": "Magnetic Potential Ator",
        "psi": "Poloidal Flux",
        "velocity_r": "Plasma Velocity Vr",
        "velocity_z": "Plasma Velocity Vz",
        "velocity_tor": "Plasma Velocity Vtor",
        "velocity_parallel": "Plasma Velocity Vparallel",
        "velocity_parallel_over_b_field": "Vparallel / |B|",
        "phi_potential": "Electric Potential",
        "vorticity": "Vorticity",
        "vorticity_over_r": "Vorticity / Major Radius",
        "j_r": "Current Density Jr",
        "j_z": "Current Density Jz",
        "j_tor": "Current Density Jtor",
        "j_tor_r": "Jtor x Major Radius",
        "mass_density": "Mass Density",
    }
    for q_name in quantities:
        units = get_units(ids_name, f"{units_path}/{q_name}")
        _add_aos_scalar_array_to_vtk_field_data(
            getattr(ggd, q_name), subset_idx, f"{quantities[q_name]} {units}", ugrid
        )

    # NOT TESTED! Can't find a data entry to test.


def read_radiation(
    ids_obj, aos_index_values: dict, subset_idx: int, ugrid: vtkUnstructuredGrid
) -> None:

    # ggd is at /radiation/process(i1)/ggd(itime)
    # TODO: ProcessIdx is not detected by generator.py and hence it is not present.
    # process_idx = aos_index_values.get('ProcessIdx') or 0
    # time_idx = aos_index_values.get('TimeIdx')
    # try:
    #     ggd = ids_obj.process[process_idx].ggd[time_idx]
    # except IndexError:
    #     return

    # for finding units:
    ids_name = "radiation"
    ggd_path = "process/ggd"
    units_path = ggd_path

    time_idx = aos_index_values.get("TimeIdx")
    for process in ids_obj.process:
        try:
            ggd = process.ggd[time_idx]
        except IndexError:
            continue

        # Use process name in brackets but only if we have more than one process:
        if len(ids_obj.process) > 1:
            p_name = f"({process.identifier.name})"
        else:
            p_name = ""

        # electrons.emissivity
        q_name = "electrons/emissivity"
        units = get_units(ids_name, f"{units_path}/{q_name}")
        name = f"Electron Emissivity {p_name} {units}"
        _add_aos_scalar_array_to_vtk_field_data(
            ggd.electrons.emissivity, subset_idx, name, ugrid
        )

        # ions, neutrals = heavy particles = hp
        for hp in [*ggd.ion, *ggd.neutral]:

            # heavy particles: emissivity
            quantity = hp.emissivity
            # for units, no need to distinguish between ions and neutrals, use ions all
            # the time:
            units_path = f"{ggd_path}/ion"
            q_name = "emissivity"
            units = get_units(ids_name, f"{units_path}/{q_name}")
            name = f"{hp.label} Emissivity {p_name} {units}"
            _add_aos_scalar_array_to_vtk_field_data(quantity, subset_idx, name, ugrid)

            # heavy particles: state
            for state in hp.state:
                quantity = state.emissivity
                # for units, no need to distinguish between ions and neutrals, use ions
                # all the time:
                units_path = f"{ggd_path}/ion/state"
                q_name = "emissivity"
                units = get_units(ids_name, f"{units_path}/{q_name}")
                name = f"{state.label} Emissivity {p_name} {units}"
                _add_aos_scalar_array_to_vtk_field_data(
                    quantity, subset_idx, name, ugrid
                )


def read_tf(
    ids_obj, aos_index_values: dict, subset_idx: int, ugrid: vtkUnstructuredGrid
) -> None:
    # A 'ggd' node is absent,
    # but there are scalar arrays.
    # /tf/field_map(itime)/
    #  - b_field_r(i1)
    #  - b_field_z(i1)
    #  - b_field_tor(i1)
    #  - a_field_r(i1)
    #  - a_field_z(i1)
    #  - a_field_tor(i1)
    time_idx = aos_index_values.get("TimeIdx")

    # for finding units:
    ids_name = "tf"
    units_path = "field_map"

    for attr_name in [
        "b_field_r",
        "b_field_z",
        "b_field_tor",
        "a_field_r",
        "a_field_z",
        "a_field_tor",
    ]:
        try:
            aos_scalar_node = getattr(ids_obj.field_map[time_idx], attr_name)
        except IndexError:
            continue
        component_name = attr_name.split("_")[-1]
        units = get_units(ids_name, f"{units_path}/{attr_name}")
        name = f"Magnetic Field B{component_name} {units}"
        _add_aos_scalar_array_to_vtk_field_data(
            aos_scalar_node, subset_idx, name, ugrid
        )


def read_transport_solver_numerics(
    ids_obj, aos_index_values: dict, subset_idx: int, ugrid: vtkUnstructuredGrid
) -> None:

    # ggd is at /transport_solver_numerics/boundary_conditions_ggd(itime)
    time_idx = aos_index_values.get("TimeIdx")
    try:
        ggd = ids_obj.boundary_conditions_ggd[time_idx]
    except IndexError:
        return

    # for finding units:
    ids_name = "transport_solver_numerics"
    ggd_path = "boundary_conditions_ggd"

    # - current(i1)
    name = f"Current Boundary Condition - {ggd.current.identifier.name.capitalize()}"
    _add_aos_scalar_array_to_vtk_field_data(ggd.current, subset_idx, name, ugrid)
    # - electrons/particles(i1)
    try:
        units_path = ggd_path + "/electrons/particles"
        units = get_units(ids_name, units_path)
        name = (
            "Electron Density Boundary Condition - "
            f"{ggd.electrons.particles[subset_idx].identifier.name} {units}"
        )
        _add_aos_scalar_array_to_vtk_field_data(
            ggd.electrons.particles, subset_idx, name, ugrid
        )
    except IndexError:
        # can raise from name = "...", when subset_idx >= ggd.electrons.particles
        pass
    # - electrons/energy(i1)
    try:
        units_path = ggd_path + "/electrons/energy"
        units = get_units(ids_name, units_path)
        name = (
            "Electron Energy Boundary Condition - "
            f"{ggd.electrons.energy[subset_idx].identifier.name} {units}"
        )
        _add_aos_scalar_array_to_vtk_field_data(
            ggd.electrons.energy, subset_idx, name, ugrid
        )
    except IndexError:
        # can raise from name = "...", when subset_idx >= ggd.electrons.energy
        pass

    # ions
    for ion in ggd.ion:
        # - ion/particles
        quantity = ion.particles
        units_path = ggd_path + "/ion/particles"
        units = get_units(ids_name, units_path)
        name = (
            f"{ion.label} Density Boundary Condition ({quantity.identifier.name})"
            f" {units}"
        )
        _add_aos_scalar_array_to_vtk_field_data(quantity, subset_idx, name, ugrid)
        # - ion/energy
        quantity = ion.energy
        units_path = ggd_path + "/ion/energy"
        units = get_units(ids_name, units_path)
        name = (
            f"{ion.label} Energy Boundary Condition ({quantity.identifier.name})"
            f" {units}"
        )
        _add_aos_scalar_array_to_vtk_field_data(quantity, subset_idx, name, ugrid)

        for state in ion.state:
            # - ion/state/particles
            quantity = state.particles
            units_path = ggd_path + "/ion/state/particles"
            units = get_units(ids_name, units_path)
            name = (
                f"{state.label} Density Boundary Condition "
                f"({quantity.identifier.name}) {units}"
            )
            _add_aos_scalar_array_to_vtk_field_data(quantity, subset_idx, name, ugrid)
            # - ion/state/energy
            quantity = state.energy
            units_path = ggd_path + "/ion/state/energy"
            units = get_units(ids_name, units_path)
            name = (
                f"{state.label} Energy Boundary Condition "
                f"({quantity.identifier.name}) {units}"
            )
            _add_aos_scalar_array_to_vtk_field_data(quantity, subset_idx, name, ugrid)


def read_wall(
    ids_obj, aos_index_values: dict, subset_idx: int, ugrid: vtkUnstructuredGrid
) -> None:

    # ggd is at /wall/description_ggd(i1)/ggd(itime)
    description_ggd_idx = aos_index_values.get("DescriptionGgdIdx")
    time_idx = aos_index_values.get("TimeIdx")

    # for finding units:
    ids_name = "wall"
    ggd_path = "description_ggd/ggd"
    units_path = ggd_path

    try:
        if description_ggd_idx is not None:
            ggd = ids_obj.description_ggd[description_ggd_idx].ggd[time_idx]
        else:
            ggd = ids_obj.description_ggd[0].ggd[time_idx]
    except IndexError:
        return

    # Scalar quantities:
    quantities = {
        "power_density": "Power Density",
        "temperature": "Temperature",
        "v_biasing": "External Electric Potencial",
        "psi": "Poloidal Flux",
        "phi_potential": "Electric Potential",
        "resistivity": "Resistivity",
    }
    for q_name in quantities:
        try:
            units = get_units(ids_name, f"{units_path}/{q_name}")
            _add_aos_scalar_array_to_vtk_field_data(
                getattr(ggd, q_name), subset_idx, f"{quantities[q_name]} {units}", ugrid
            )
        except AttributeError:
            logger.info(
                f"No quantity {q_name} found, "
                "perhaps mismatched data dictionary versions?"
            )
        except IndexError:
            logger.info(f"No subset {subset_idx} found for quantity{q_name}.")

    # Vector quantities:
    quantities = {
        "j_total": "Current Density (total)",
        "e_field": "Electric Field",
        "a_field": "Magnetic Potential",
    }
    for q_name in quantities:
        try:
            units = get_units(ids_name, f"{units_path}/{q_name}")
            _add_aos_vector_array_to_vtk_field_data(
                getattr(ggd, q_name), subset_idx, f"{quantities[q_name]} {units}", ugrid
            )
        except AttributeError:
            logger.info(
                f"No quantity {q_name} found, "
                "perhaps missmatched data dictionary versions?"
            )
        except IndexError:
            logger.info(f"No subset {subset_idx} found for quantity{q_name}.")


def read_waves(
    ids_obj, aos_index_values: dict, subset_idx: int, ugrid: vtkUnstructuredGrid
) -> None:

    # A 'ggd' node is absent,
    # but there are scalar and vector arrays.
    # /waves/coherent_wave(i1)/full_wave(itime)/
    #  - e_field/plus(i2)
    #  - e_field/minus(i2)
    #  - e_field/parallel(i2)
    #  - e_field/normal(i2)
    #  - e_field/bi_normal(i2)
    #  - b_field/plus(i2)
    #  - b_field/minus(i2)
    #  - b_field/parallel(i2)
    #  - b_field/normal(i2)
    #  - b_field/bi_normal(i2)
    #  - k_perpendicular(i2)
    coherent_wave_idx = aos_index_values.get("CoherentWaveIdx")
    time_idx = aos_index_values.get("TimeIdx")

    # for finding units:
    ids_name = "waves"
    ggd_path = "coherent_wave/full_wave"

    try:
        e_field = ids_obj.coherent_wave[coherent_wave_idx].full_wave[time_idx].e_field
        units_path = ggd_path + "/e_field"
        quantities = {
            "plus": "Electric Field - LH polarized",
            "minus": "Electric Field - RH polarized",
            "parallel": "Electric Field - Parallel to B",
            "normal": "Electric Field - Normal to flux surface",
            "bi_normal": "Electric Field - Tangential to flux surface",
        }
        for q_name in quantities:
            units = get_units(ids_name, f"{units_path}/{q_name}")
            _add_aos_scalar_array_to_vtk_field_data(
                getattr(e_field, q_name),
                subset_idx,
                f"{quantities[q_name]} {units}",
                ugrid,
            )
    except IndexError:
        pass

    try:
        b_field = ids_obj.coherent_wave[coherent_wave_idx].full_wave[time_idx].b_field
        units_path = ggd_path + "/b_field"
        quantities = {
            "parallel": "Magnetic Field - Parallel to B",
            "normal": "Magnetic Field - Normal to flux surface",
            "bi_normal": "Magnetic Field - Tangential to flux surface",
        }
        for q_name in quantities:
            units = get_units(ids_name, f"{units_path}/{q_name}")
            _add_aos_scalar_array_to_vtk_field_data(
                getattr(b_field, q_name),
                subset_idx,
                f"{quantities[q_name]} {units}",
                ugrid,
            )
    except IndexError:
        pass

    try:
        k_perp = (
            ids_obj.coherent_wave[coherent_wave_idx].full_wave[time_idx].k_perpendicular
        )
        units_name = ggd_path + "/k_perpendicular"
        units = get_units(ids_name, f"{units_name}")
        name = f"Perpendicular Wave Vector {units}"
        _add_aos_scalar_array_to_vtk_field_data(k_perp, subset_idx, name, ugrid)
    except IndexError:
        pass


def _add_scalar_array_to_vtk_field_data(
    array: np.ndarray, name: str, ugrid: vtkUnstructuredGrid
) -> None:
    """
    Add a named array as scalars to a vtkUnstructuredGrid instance
    :param array: the numpy array.
    :param name: the name string.
    :param ugrid: an instance of vtkUnstructuredGrid
    :return: None
    """
    logger.debug(f"           {name}...")
    point_data: vtkPointData = ugrid.GetPointData()
    num_points = ugrid.GetNumberOfPoints()
    cell_data: vtkCellData = ugrid.GetCellData()
    num_cells = ugrid.GetNumberOfCells()

    # To see interpolation of point data on cells, just point data is necessary.
    vtk_arr = dsa.numpyTovtkDataArray(array, name)
    if len(array) == num_points:
        point_data.AddArray(vtk_arr)
    if len(array) == num_cells:
        cell_data.AddArray(vtk_arr)


def _add_aos_scalar_array_to_vtk_field_data(
    aos_scalar_node, subset_idx: int, name: str, ugrid: vtkUnstructuredGrid
):
    """
    Add the array under the aos_scalar_node to the unstructured grid.
    :param aos_scalar_node: A node with scalar values for each grid subset.
    :param subset_idx: an index into aos_scalar_node
    :param name: this becomes the array name in VTK
    :param ugrid: an unstructured grid instance
    :return: None
    """
    logger.debug(f"           {name}...")
    if subset_idx >= len(aos_scalar_node):
        return

    # For wall IDS nodes, edges, cells, volumes in one partition.
    if subset_idx == -1:
        for i in range(4):
            try:
                if hasattr(aos_scalar_node[i], "values") and len(
                    aos_scalar_node[i].values
                ):
                    _add_scalar_array_to_vtk_field_data(
                        aos_scalar_node[i].values, name, ugrid
                    )
            except IndexError:
                logger.info(f"           no index {i} for subset {subset_idx}...")
            except AttributeError:
                logger.info(f"           no index {i} for subset {subset_idx}...")
    else:
        if hasattr(aos_scalar_node[subset_idx], "values") and len(
            aos_scalar_node[subset_idx].values
        ):
            _add_scalar_array_to_vtk_field_data(
                aos_scalar_node[subset_idx].values, name, ugrid
            )


def _add_aos_vector_array_to_vtk_field_data(
    aos_vector_node, subset_idx: int, name: str, ugrid: vtkUnstructuredGrid
):
    """
    Add the array under the aos_vector_node to the unstructured grid.
    :param aos_vector_node: A node with component vectors for each grid subset.
    :param subset_idx: an index into aos_vector_node
    :param name: this becomes the array name in VTK
    :param ugrid: an unstructured grid instance
    :return: None
    """
    logger.debug(f"           {name}...")
    if subset_idx >= len(aos_vector_node):
        return

    point_data: vtkPointData = ugrid.GetPointData()
    num_points = ugrid.GetNumberOfPoints()
    cell_data: vtkCellData = ugrid.GetCellData()
    num_cells = ugrid.GetNumberOfCells()

    # Only add the components that have data:
    components = dict()  # name and values
    for component_name in [
        "radial",
        "diamagnetic",
        "parallel",
        "poloidal",
        "toroidal",
        "r",
        "z",
    ]:
        try:
            values = getattr(aos_vector_node[subset_idx], component_name)
        except Exception:
            continue
        if len(values):
            components[component_name] = values

    vtk_arr = vtkDoubleArray()
    vtk_arr.SetName(name)
    vtk_arr.SetNumberOfComponents(len(components))
    num_tuples = 0

    for i, component_name in enumerate(components):
        vtk_arr.SetComponentName(i, component_name.capitalize())
        scalar_arr = dsa.numpyTovtkDataArray(
            components[component_name], name + "-" + component_name.capitalize()
        )
        if num_tuples == 0:
            num_tuples = scalar_arr.GetNumberOfTuples()
            vtk_arr.SetNumberOfTuples(num_tuples)

        vtk_arr.CopyComponent(i, scalar_arr, 0)

    if num_tuples == num_points:
        point_data.AddArray(vtk_arr)
    if num_tuples == num_cells:
        cell_data.AddArray(vtk_arr)
