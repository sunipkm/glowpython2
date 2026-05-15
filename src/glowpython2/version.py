from importlib.metadata import version

try:
    __version__ = version(__package__ or __name__)
except Exception:
    __version__ = "unknown"
