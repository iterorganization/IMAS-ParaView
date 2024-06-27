def profile(fname):
    """Decorator to cProfile the code inside the function.

    Args:
        fname: Filename to dump the cProfile statistics to.

    Example:
        .. code-block:: python

            class Test(VTKPythonAlgorithmBase):
                @profile("RequestData.stats")
                def RequestData(self, request, inInfo, outInfo):
                    ...
    """
    import cProfile
    import functools

    def inner(func):
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            prof = cProfile.Profile()
            try:
                prof.enable()
                return func(*args, **kwargs)
            finally:
                prof.disable()
                prof.dump_stats(fname)
                print("Statistics written to", fname)

        return wrapper

    return inner
