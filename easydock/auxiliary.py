from itertools import islice
from math import ceil

def take(n, iterable):
    return list(islice(iterable, n))


def mol_name_split(mol_name):
    return mol_name.rsplit('_', 1)


def empty_func(*args, **kwargs):
    pass


def empty_generator(*args, **kwargs):
    return
    yield

def chunk_into_n(smi_l: list[str], n: int):
    smi_size = ceil(len(smi_l) / n)
    return list(map(lambda x: smi_l[x * smi_size:x * smi_size + smi_size], list(range(n))))

def chunk_into_size_n(smi_l: list[str], n: int):
    return [smi_l[i:i+n] for i in range(0,len(smi_l),n)]