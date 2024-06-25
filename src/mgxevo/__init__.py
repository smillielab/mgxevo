# from importlib.metadata import version
# __init__.py from jsonschema 3.2.0
try:
    from importlib.metadata import version
except ImportError: # for Python<3.8
    from importlib_metadata import version

from . import ut, tl, pl, data

__all__ = ["ut", "tl", "pl", "data"]

# __version__ = version("mgxevo")


