import os
import argparse

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

