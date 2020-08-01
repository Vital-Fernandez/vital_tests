import timeit

testcode1 = '''
def myCoso1(a, b, c):
    return a + b + c

a, b, c = 2, 3, 4

myCoso1(2,3,4)    
'''

testcode2 = '''
def myCoso1(a, b, **kwargs):
    return a + b + kwargs['c']

a, b, c = 2, 3, 4

myCoso1(2,3,c=c)    
'''

testcode3 = '''
def myCoso1(a, b, **kwargs):
    return a + b + kwargs['c']+ kwargs['d']

a, b, c, d = 2, 3, 4, 5

myCoso1(2,3,c=c,d=d)   
'''


def myCoso2(a, b, c):
    return a + b + c


# myCoso1(2, 3, c=4)
# myCoso2(2, 3, 4)

print(timeit.timeit(stmt=testcode1, number=1000000))
print(timeit.timeit(stmt=testcode2, number=1000000))
print(timeit.timeit(stmt=testcode3, number=1000000))
