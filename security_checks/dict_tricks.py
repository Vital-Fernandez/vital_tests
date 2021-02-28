import numpy as np

AMP_CONF_orig = dict(value=None, min=0, max=np.inf, stderr=None, vary=True, expr=None, brute_Step=None)
MU_CONF_orig = dict(value=None, min=-np.inf, max=np.inf, stderr=None, vary=True, expr=None, brute_Step=None)
SIG_CONF_orig = dict(value=None, min=0, max=np.inf, stderr=None, vary=True, expr=None, brute_Step=None)

AMP_CONF_DEFAULT = dict(value=None, min=0, max=np.inf, stderr=None, vary=True, expr=None, brute_Step=None)
MU_CONF_DEFAULT = dict(value=None, min=-np.inf, max=np.inf, stderr=None, vary=True, expr=None, brute_Step=None)
SIG_CONF_DEFAULT = dict(value=None, min=0, max=np.inf, stderr=None, vary=True, expr=None, brute_Step=None)

myConf = dict(value=3.0, vary=False, mierdaNewville=True)

print('\nOption 1')
print(AMP_CONF_DEFAULT == AMP_CONF_orig)
print('myConf', myConf)
newConf = {**AMP_CONF_DEFAULT, **myConf}
print(newConf)
print(AMP_CONF_DEFAULT == AMP_CONF_orig)
print(newConf==AMP_CONF_DEFAULT, newConf==AMP_CONF_DEFAULT)

print('\nOption 2')
print(AMP_CONF_DEFAULT == AMP_CONF_orig)
print('myConf', myConf)
print({**myConf, **AMP_CONF_DEFAULT})
print(AMP_CONF_DEFAULT == AMP_CONF_orig)

print('\nOption 3')
print(AMP_CONF_DEFAULT == AMP_CONF_orig)
emptyConf = {}
print('emptyConf', emptyConf)
print({**AMP_CONF_DEFAULT, **emptyConf})
print(AMP_CONF_DEFAULT == AMP_CONF_orig)


for i, item in enumerate(AMP_CONF_orig.items()):
   key, value = item
   print(i, key, value)

print()
print(myConf)
print(list(myConf))
a = np.array(myConf)
print(a)