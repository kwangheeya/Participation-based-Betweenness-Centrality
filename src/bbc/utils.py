from collections import Iterable
from contextlib import contextmanager
import traceback


def print_status(*obj, print_op=True):
    if print_op & isinstance(obj, Iterable):
        for i in obj:
            print(i, end=' ')
        print()


def get_function_name():
    return traceback.extract_stack(None, 2)[0][2]


@contextmanager
def fopen(filepath, mode='r'):
    try:
        f = open(filepath, mode)
        print('+ Open file [{0}]'.format(filepath))
        yield f
    except IOError as err:
        print('- Error: can\'t find [{0}] or read data ({1})'.format(filepath, err))
    finally:
        f.close()

