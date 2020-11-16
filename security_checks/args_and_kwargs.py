
# concatenate_keys.py
def myFun(arg1, **kwargs):
    if 'Bicho' in kwargs:
        print('No bicho')
    if 'first' in kwargs:
        print('first si')
    for key, value in kwargs.items():
        print ("%s == %s" %(key, value))

myFun("Hi", first ='Geeks', mid ='for', last='Geeks')
