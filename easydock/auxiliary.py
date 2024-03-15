from itertools import islice


def take(n, iterable):
    return list(islice(iterable, n))


def mol_name_split(mol_name):
    return mol_name.rsplit('_', 1)


def empty_func(*args, **kwargs):
    pass


def empty_generator(*args, **kwargs):
    return
    yield
