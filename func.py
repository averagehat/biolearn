from functools import partial, wraps
import itertools as it
import string
import sys
from collections import namedtuple
from operator import itemgetter, attrgetter as attr
from schema import Schema, SchemaError

PY3 = sys.version[0] == '3'
imap, ifilter, izip = (map, filter, zip) if PY3 else (it.imap, it.ifilter, it.izip)

#notin = compose(_not, operator.methodcaller('__contains__'))
#notin = compose(_not, attr('__contains__'))
#mismatches = pfilter(notin('M='))
def merge_dicts(*dict_args):
    '''
    from http://stackoverflow.com/a/26853961
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    '''
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

def _not(x):
    return not x

def partial2(method, param):
      def t(x):
              return method(x, param)
      return t

def _id(x): return x

def apply_key_func(k, v, funcdict):
    return funcdict.get(k, _id)(v)
def compose_all(*funcs):
    return reduce(compose, funcs)

#def k_compose(outer, **okwargs):
#    ''' compose(f, g)(x) == f(g(x)) '''
#    def newfunc(*args, **ikwargs):
#        _kwargs = dict( (k, apply_key_func(k, v, okwargs)) for k, v in ikwargs.items())
#        return outer(*args, **_kwargs)
#    return newfunc

def compose(outer, inner):
    ''' compose(f, g)(x) == f(g(x)) '''
    def newfunc(*args, **kwargs):
        return outer(inner(*args, **kwargs))
    return newfunc


starcompose2 = lambda f, g: lambda x: f(*g(x))

#starcompose = partial(reduce, starcompose2) #need splat
def starcompose(*funcs):
    return reduce(starcompose2, funcs)

def compose(outer, inner):
    ''' compose(f, g)(x) == f(g(x)) '''
    def newfunc(*args, **kwargs):
        return outer(inner(*args, **kwargs))
    return newfunc
def fzip(funcs, args):
    for func, arg in izip(funcs, args):
        yield func(arg)

def dictmap(func, _dict):
    return dict( (key, func(val)) for key, val in _dict.items())


def ilen(iterable):
    return sum(1 for _ in iterable)

def reverse(collection): return collection[::-1]
pifilter = partial(partial, ifilter)
compose_list = partial(reduce, compose)
#   compose_all = compose(compose_list, lambda *a: a)
pmap = partial(partial, map)
pfilter = partial(partial, filter)
#TODO: could use partial2 instead
pstrip = lambda x: partial(string.strip, chars=x)
psplit = lambda x: partial(string.split, sep=x)
pjoin = lambda x: partial(string.join, sep=x)
boolint = lambda x: 1 if x else 0
dictzip = compose(dict, zip)
cmp2=lambda f, g: lambda a, b: (f(a), g(b))
#ilen = compose(sum, pmap(boolint))
#Given a list of functons and names, return the result of those functions dictzipped witht the names.

#TODO:
''' dictfilter '''
def apply_each(funcs, arg):
    return fzip(funcs, it.repeat(arg))

import inspect
import types
def is_local(object):
   return isinstance(object, types.FunctionType) and object.__module__ == __name__
    #use inspect.isfunction
def get_funcs():
    return inspect.getmembers(sys.modules[__name__], \
                             predicate = lambda f: inspect.isfunction(f) and f.__module__ == __name__)
    #return inspect.getmembers(sys.modules[__name__], predicate=is_local)
    #return dict( ((name, func)) for name, func in locals().items() if is_local(name))
#      for key, value in locals().items():
#          if callable(value) and value.__module__ == __name__:
#              l.append(key)


'''
compose columns + object + getters => dict
- unzip

have fqframe return a dict of functions, excluding get_row, openframe; instead passing it to a function which
arranges the getters, applies them to a get_object function, and creates an intermediate dictionary.
This function allows for optional extras, like samframe & fastq (rather than fasta
'''


unzip = starcompose(zip, _id)

def nameddict(Name, _dict):
    ''' dict to named tuple '''
    names, values = unzip(_dict.items())
    return namedtuple(Name, names)(*values)

ppartial = partial(partial)
apply_to_object = compose(apply, ppartial)

kstarcompose2 = lambda f, g: lambda x: f(**g(x))
def kstarcompose(*funcs):
    return reduce(kstarcompose2, funcs)

#kstarcompose = partial(reduce, kstarcompose2)

#use str.endswith( (tuple, of, vals)
extension = compose(itemgetter(-1), psplit('.'))
fileext = compose(extension, attr('filename'))


def iter_until_stop(f, *args, **kwargs):
    while True:
        try:
            yield f(*args, **kwargs)
        except StopIteration:
            break

flatten_list = lambda a: a if type(a) != list else a[0]

def split_list(A, idx):
    return A[:idx], A[idx:]


'''http://docs.python.org/3.4/library/itertools.html#itertools-recipes'''


def partition(pred, iterable):
    """Use a predicate to partition entries into false entries and true entries
    partition(is_odd, range(10)) --> 0 2 4 6 8 and 1 3 5 7 9
    http://docs.python.org/3.4/library/itertools.html#itertools-recipes
    """
    t1, t2 = it.tee(iterable)
    return it.ifilterfalse(pred, t1), it.ifilter(pred, t2)


_and2 = lambda f, g: lambda x: f(x) and g(x)
_and = lambda *x: reduce(_and2, x)
_not = lambda f: lambda x: not f(x)
_or = lambda f, g: lambda x: f(x) or g(x)
bnot = lambda x: ~x


#def pool_map(func, *args, **kwargs):
#    pool = multiprocessing.Pool()
#    pool
#
def rpartial(func, *args):
    return lambda *a: func(*(a + args))



message = "in function {0}, arg {1} is not an instance of {2}, but is type {3}.".format
def check_type(func, _type, obj):
    assert isinstance(obj, _type), message(func.__name__, obj, _type, type(obj))

#get_type = partial(getattr, __module__)
def dict_intersect(d1, d2):
    #return dict(d1.viewitems() & d2.viewitems())
    return {x:d1[x] for x in d1 if x in d2}
functype = 'FUNCTYPE'

def handle_check(etype, arg, kwargs, func):
    if type(etype) is dict:
        scheme, matched_args = dict_intersect(etype, kwargs), dict_intersect(kwargs, etype)
        try:
            Schema(scheme).validate(matched_args)#, error=message(func.__name__, arg, etype))
        except SchemaError as e:
            print "scheme {0} did not fit kwargs {1}".format(scheme, kwargs)
            raise e
    elif arg is None or etype is None:
        #TODO: log?
        #print "Warning " + message(func.__name__, arg, etype)
        return
    elif etype == functype:
        assert hasattr(arg, '__call__')
    else:
        check_type(func, etype, arg)
    #TODO:
    '''
    from string import splitfields
    if type(etype) is tuple:
       handle *args (how?)

    if etype.startswith('has'):
        attrs = etype.split(' ')[1:]
        for _attr in attrs:
            assert hasattr(arg, _attr), "arg {0} did not have attribute {1} in function {2}".format(arg, _attr, func.__name__)
'''

''' could be: *types, **kwtypes '''
def typecheck(*types):
    def decorator(func):
        @wraps(func)
        def typechecked(*args, **kwargs):
            check = partial(handle_check, func=func, kwargs=kwargs)
            map(check, types, args)
            return func(*args, **kwargs)
        return typechecked
    return decorator

