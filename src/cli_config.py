"""Select the configuration file from `--config PATH` before `params` is imported.

params.py reads the file named by the TIDAL_NS_CONFIG environment variable at import time. The
command-line tools call `select_from_argv()` first thing, so that `--config` works for every
module that imports params.
"""
import os
import sys


def select_from_argv(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    for i, a in enumerate(argv):
        if a == "--config" and i + 1 < len(argv):
            os.environ["TIDAL_NS_CONFIG"] = argv[i + 1]
        elif a.startswith("--config="):
            os.environ["TIDAL_NS_CONFIG"] = a.split("=", 1)[1]
