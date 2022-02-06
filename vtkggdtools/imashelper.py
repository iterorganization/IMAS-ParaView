"""
Exception reinforced imas open/create
"""

from vtkggdtools.errors import InvalidIDSIOError


def imas_env_call(ids_func, *ids_args):
    """
    Open/Create an imas ids.
    :param ids_func: an io function (imas.ids.open_env or im)
    :param shot: a shot number
    :type shot: int
    :param run: a run number
    :type run: int
    :param user: username
    :type user: str
    :param tokamak: a tokamak populated with idss
    :type tokamak: str
    :param imasver: major version of IMAS
    :type imasver: str
    :param backend: backend type
    :type backend: int
    :return: None
    :raises InvalidIDSIOError
    """
    # imas.ids.open_env returns a tuple (status, idx)
    # The return codes are detailed in a comment in this jira issue.
    # [https://jira.iter.org/browse/IMAS-2859?focusedCommentId=668912&page=com.atlassian.jira.plugin.system.issuetabpanels:comment-tabpanel#comment-668912]
    # The codes:
    # status = 0, success
    # status = -2 context exception
    # status = -3 backend exception
    # status = -4 Abstract Low Layer exception
    # status = -1 unknown error
    # idx = 0 or positive, success
    # idx = -2 context exception
    # idx = -3 backend exception
    # idx = -4 Abstract Low Layer exception
    # idx = -1 unknown error
    status, idx = ids_func(*ids_args)
    if (status != 0) or (idx < 0):
        message = "Unknown exception was raised."
        if status == -2:
            message = "Context exception was raised."
        elif status == -3:
            message = "Backend exception was raised."
        elif status == -4:
            message = "Abstract Low Layer exception was raised."
        elif status == -1:
            message = "Unknown exception was raised."
        if idx == -2:
            message = "Context exception was raised."
        elif idx == -3:
            message = "Backend exception was raised."
        elif idx == -4:
            message = "Abstract Low Layer exception was raised."
        elif idx == -1:
            message = "Unknown exception was raised."
        e = InvalidIDSIOError(message, ids_func, *ids_args)
        raise e
