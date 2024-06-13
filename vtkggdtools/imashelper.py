"""
Set of functions to help using the Access Layer API
"""

# Reading units from an IDS path:
from functools import lru_cache

from imaspy import IDSFactory
from paraview import logger as pvlog

# Units pre- and post- formatting:
u_pre = "["
u_post = "]"


# For MathML:
# u_pre = '[$' # start MathText
# u_post = '$]' # end MathText
@lru_cache
def get_units(ids_name: str, path: str, pre=u_pre, post=u_post) -> str:
    """
    Read the units from an IDS path.
    : param ids_name: name of the IDS to read the units from (eg, "equilibrium")
    : type ids_name: str
    : param path: the path of the entry to get the units
        (eg, "time_slice(:)/profiles_2d(0)/psi")
    : type path: str
    : return: a string with the units, wrapped in u_pre and u_post
    """
    try:
        # FIXME: using default IDS Factory as temporary workaround
        # A more permanent solution will use ``node.metadata.units`` directly instead of
        # this method
        units = f"{pre}{IDSFactory().new(ids_name).metadata[path].units}{post}"
    except Exception:
        pvlog.warn(f"Can't read units for {ids_name}/{path}.")
        units = ""
    return units
