from itertools import islice
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