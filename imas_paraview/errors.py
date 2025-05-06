"""
Error handling for imas open/create or ids related get/put
"""


class InvalidIDSIOError(Exception):
    def __init__(self, msg, ids_func, *ids_args):
        self.message = (
            f"Error message: {msg} | function {ids_func}, args: {list(ids_args)}"
        )

    def __str__(self):
        return self.message


class InvalidIDSParametersError(Exception):
    def __init__(self, **params):
        self.message = (
            f"Invalid IDS with parameters {params if len(params) else 'None'}"
        )

    def __str__(self):
        return self.message


class InputPointsError(Exception):
    def __init__(self, *args, **kwargs):
        pass
