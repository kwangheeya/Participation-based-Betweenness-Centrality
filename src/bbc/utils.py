from collections import Iterable
import traceback


def print_status(*obj, print_op=True):
    if print_op & isinstance(obj, Iterable):
        for i in obj:
            print(i, end=' ')
        print()

def get_function_name():
    return traceback.extract_stack(None, 2)[0][2]

