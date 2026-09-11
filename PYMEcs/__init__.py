# __version__ = '26.09.02'
import importlib.metadata

try:
    # figure out version from installed metadata
    __version__ = importlib.metadata.version("PYME-extra")
except importlib.metadata.PackageNotFoundError:
    # Fallback if the package is not installed (e.g., during local development)
    try:
        from ._version import __version__
    except ImportError:
        __version__ = "0.0.0.dev0+unknown"
