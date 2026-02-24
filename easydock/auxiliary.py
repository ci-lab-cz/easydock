import gzip
import os
import pickle
from itertools import islice
from math import ceil
import signal
import functools
import platform
import threading
from typing import List, Tuple, Optional


def take(n, iterable):
    return list(islice(iterable, n))


def mol_name_split(mol_name):
    return tuple(mol_name.rsplit('_', 1))


def empty_func(*args, **kwargs):
    pass


def empty_generator(*args, **kwargs):
    return
    yield


#https://stackoverflow.com/questions/75928586/how-to-stop-the-execution-of-a-function-in-python-after-a-certain-time/75928879#75928879
def timeout(seconds=5, default=None):

    def decorator(func):

        @functools.wraps(func)
        def wrapper(*args, **kwargs):

            if platform.system() == 'Windows':
                # Windows fallback: use threading
                result = {}

                def target():
                    try:
                        result["value"] = func(*args, **kwargs)
                    except Exception as e:
                        result["error"] = e

                t = threading.Thread(target=target)
                t.daemon = True
                t.start()
                t.join(seconds)

                if t.is_alive():
                    return None
                    # raise TimeoutError()
                # if "error" in result:
                #     raise result["error"]
                return result.get("value", None)

            else:

                # Unix: use SIGALRM
                def handle_timeout(signum, frame):
                    raise TimeoutError()

                signal.signal(signal.SIGALRM, handle_timeout)
                signal.alarm(seconds)

                result = func(*args, **kwargs)

                signal.alarm(0)

                return result

        return wrapper

    return decorator


def split_generator_to_chunks(generator, chunk_size):
    """Yield successive chunks from a generator"""
    chunk = []

    for item in generator:
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = [item]
        else:
            chunk.append(item)

    if chunk:
        yield chunk


def chunk_into_n(smi_l: List[str], n: int):
    smi_size = ceil(len(smi_l) / n)
    return list(map(lambda x: smi_l[x * smi_size:x * smi_size + smi_size], list(range(n))))


def chunk_into_size_n(smi_l: List[str], n: int):
    return [smi_l[i:i+n] for i in range(0,len(smi_l),n)]


def expand_path(path):
    return os.path.expanduser(os.path.expandvars(path))


def count_input_structures(input_fname: str) -> Optional[int]:
    """
    Count all records in input file, including structures that may fail RDKit parsing.
    """
    def _count_sdf_structures(handle) -> int:
        count = 0
        has_data = False
        for line in handle:
            if line.strip():
                has_data = True
            if line.strip() == '$$$$':
                count += 1
                has_data = False
        if has_data:
            count += 1
        return count

    if not input_fname:
        return None

    fname = os.path.abspath(input_fname)
    lower_fname = fname.lower()

    try:
        if lower_fname.endswith('.smi') or lower_fname.endswith('.smiles'):
            with open(fname) as f:
                return sum(1 for line in f if line.strip())

        if lower_fname.endswith('.sdf'):
            with open(fname) as f:
                return _count_sdf_structures(f)

        if lower_fname.endswith('.sdf.gz'):
            with gzip.open(fname, 'rt', encoding='utf-8', errors='ignore') as f:
                return _count_sdf_structures(f)

        if lower_fname.endswith('.pkl'):
            count = 0
            with open(fname, 'rb') as f:
                while True:
                    try:
                        pickle.load(f)
                        count += 1
                    except EOFError:
                        break
            return count
    except OSError:
        return None

    return None
