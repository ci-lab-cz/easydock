import os
import argparse
from multiprocessing import cpu_count

protonation_programs = ['chemaxon', 'pkasolver', 'molgpka']

def protonation_type(value: str):
    allowed = protonation_programs
    if value in allowed:
        return value
    elif value.endswith(".sif") and os.path.isfile(value):
        return value
    else:
        raise argparse.ArgumentTypeError(
            f"Invalid --protonation argument '{value}'. "
            f"Must be one of {sorted(allowed)} or an existing .sif file."
        )


def cpu_type(x):
    return max(1, min(int(x), cpu_count()))


def filepath_type(x):
    if x:
        return os.path.abspath(x)
    else:
        return x


def str_lower_type(x):
    if x:
        return x.lower()
    else:
        return x
