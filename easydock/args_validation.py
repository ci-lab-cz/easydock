import logging
import os
import argparse
from multiprocessing import cpu_count
from easydock.auxiliary import expand_path as _expand_path

protonation_programs = ['chemaxon', 'pkasolver', 'molgpka', 'molgpka_fix']

def protonation_type(value: str):
    if value.endswith(".sif") and os.path.isfile(_expand_path(value)):
        return _expand_path(value)
    return value


def cpu_type(x):
    return max(1, min(int(x), cpu_count()))


def filepath_type(x):
    if x:
        return os.path.abspath(_expand_path(x))
    else:
        return x


def str_lower_type(x):
    if x:
        return x.lower()
    else:
        return x


def log_level_type(value):
    level = logging.getLevelName(value.upper())
    if not isinstance(level, int):
        raise argparse.ArgumentTypeError(
            f"Invalid log level: {value!r}. Choose from: NOTSET, DEBUG, INFO, WARNING, ERROR, CRITICAL."
        )
    return level
