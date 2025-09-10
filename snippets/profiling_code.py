import itertools
import timeit
from random import random

bad_dict = {'[':None, ']':None}

def uz(lelist, lestring):
    for x in lelist:
        if lestring.count(x):
            return 'Yep. "%s" contains characters from "%s" item.' % (lestring, x)

def ab(lelist, lestring):
    return [e for e in lelist if e in lestring]


def ab_2(lelist):
    # return [e in lestring for e in lelist if e in lestring]
    # return len([e for e in lelist if e in lestring])
    return lelist.translate(str.maketrans(bad_dict))

def jamie(lelist, lestring):
    return next(itertools.chain((e for e in lelist if e in lestring), (None,))) is not None


lestring = "Text123"
lelist = ["Text", "foo", "bar"]
line_name = '[O3]_5007A'

t_ab = timeit.Timer("ab(lelist, lestring)", setup="from __main__ import lelist, lestring, ab")
t_uz = timeit.Timer("uz(lelist, lestring)", setup="from __main__ import lelist, lestring, uz")
t_ab2 = timeit.Timer("ab_2(line_name)", setup="from __main__ import line_name, bad_dict, ab_2")


print(t_ab.timeit(100000))
print(t_uz.timeit(100000))
print(t_ab2.timeit(100000))
