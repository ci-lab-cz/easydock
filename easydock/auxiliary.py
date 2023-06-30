from itertools import islice


def take(n, iterable):
    return list(islice(iterable, n))
