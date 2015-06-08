from schema import Schema, SchemaError
from itertools import izip, repeat
from functools import partial
import inspect
from collections import namedtuple
functype = 'FUNCTYPE'
def fzip(funcs, args): #Embedded file name: func.py
    for func, arg in izip(funcs, args):
        yield func(arg)

def _and(*x):
    #Embedded file name: func.py
    return reduce(_and2, x)


def _and2(f, g):
    #Embedded file name: func.py
    return lambda x: f(x) and g(x)


def _id(x):
    #Embedded file name: func.py
    return x


def _not(f):
    #Embedded file name: func.py
    return lambda x: not f(x)


def _or(f, g):
    #Embedded file name: func.py
    return lambda x: f(x) or g(x)


def apply_each(funcs, arg):
    #Embedded file name: func.py
    return fzip(funcs, repeat(arg))


def apply_key_func(k, v, funcdict):
    #Embedded file name: func.py
    return funcdict.get(k, _id)(v)




def _bnot(x):
    #Embedded file name: func.py
    return ~x


#def <lambda>(x):
#    #Embedded file name: func.py
#    if x:
#        return 1
#    return 0

message = "Function {0} argument {1} expected type {2} but was type {3}.".format
def check_type(func, _type, obj):
    #Embedded file name: func.py
    assert isinstance(obj, _type), message(func.__name__, obj, _type, type(obj))


def cmp2(f, g):
    #Embedded file name: func.py
    return lambda a, b: (f(a), g(b))


def compose2(outer, inner):
    #Embedded file name: func.py
    def newfunc(*args, **kwargs):
        return outer(inner(*args, **kwargs))
    return newfunc

def compose_all(*funcs):
    return reduce(compose, funcs)

compose = compose_all

def dict_intersect(d1, d2):
    #Embedded file name: func.py
    return {x:d1[x] for x in d1 if x in d2}


def dictmap(func, _dict):
    #Embedded file name: func.py
    return dict(((key, func(val)) for key, val in _dict.items()))


def flatten_list(a):
    #Embedded file name: func.py
    if type(a) != list:
        return a
    return a[0]




def get_funcs():
    #Embedded file name: func.py
    return inspect.getmembers(sys.modules[__name__], predicate=lambda f: inspect.isfunction(f) and f.__module__ == __name__)


def handle_check(etype, arg, kwargs, func):
    #Embedded file name: func.py
    if type(etype) is dict:
        scheme, matched_args = dict_intersect(etype, kwargs), dict_intersect(kwargs, etype)
        try:
            Schema(scheme).validate(matched_args)
        except SchemaError as e:
            print 'scheme {0} did not fit kwargs {1}'.format(scheme, kwargs)
            raise e

    else:
        if arg is None or etype is None:
            return
        if etype == functype:
            assert hasattr(arg, '__call__')
        else:
            check_type(func, etype, arg)


def ilen(iterable):
    return sum((1 for _ in iterable))


def is_local(object):
    return isinstance(object, types.FunctionType) and object.__module__ == __name__


def iter_until_stop(f, *args, **kwargs):
    #Embedded file name: func.py
    while True:
        try:
            yield f(*args, **kwargs)
        except StopIteration:
            break


def kstarcompose(*funcs):
    #Embedded file name: func.py
    return reduce(kstarcompose2, funcs)


def kstarcompose2(f, g):
    #Embedded file name: func.py
    return lambda x: f(**g(x))


def merge_dicts(*dict_args):
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def nameddict(Name, _dict):
    names, values = unzip(_dict.items())
    return namedtuple(Name, names)(*values)


def partial2(method, param):
    def t(x):
        return method(x, param)
    return t


def partition(pred, iterable):
    t1, t2 = it.tee(iterable)
    return (it.ifilterfalse(pred, t1), it.ifilter(pred, t2))


def pjoin(x):
    #Embedded file name: func.py
    return partial(string.join, sep=x)


def psplit(x):
    #Embedded file name: func.py
    return partial(string.split, sep=x)


def pstrip(x):
    #Embedded file name: func.py
    return partial(string.strip, chars=x)


def reverse(collection):
    return collection[::-1]


def rpartial(func, *args):
    return lambda *a: func(*(a + args))


def split_list(A, idx):
    #Embedded file name: func.py
    return (A[:idx], A[idx:])


def starcompose(*funcs):
    #Embedded file name: func.py
    return reduce(starcompose2, funcs)


def starcompose2(f, g):
    #Embedded file name: func.py
    return lambda x: f(*g(x))


def typecheck(*types):
    def decorator(func):

        @wraps(func)
        def typechecked(*args, **kwargs):
            check = partial(handle_check, func=func, kwargs=kwargs)
            map(check, types, args)
            return func(*args, **kwargs)

        return typechecked


    return decorator

