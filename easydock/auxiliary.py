from itertools import islice
from math import ceil
import signal
import functools

def take(n, iterable):
    return list(islice(iterable, n))


def mol_name_split(mol_name):
    return mol_name.rsplit('_', 1)


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


def chunk_into_n(smi_l: list[str], n: int):
    smi_size = ceil(len(smi_l) / n)
    return list(map(lambda x: smi_l[x * smi_size:x * smi_size + smi_size], list(range(n))))


def chunk_into_size_n(smi_l: list[str], n: int):
    return [smi_l[i:i+n] for i in range(0,len(smi_l),n)]