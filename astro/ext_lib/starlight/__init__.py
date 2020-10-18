'''
Created on Mar 7, 2012

@author: Andre Luiz de Amorim

'''

from distutils.version import LooseVersion
from atpy import __version__ as atpy_version

if LooseVersion(atpy_version) >= LooseVersion('0.9.6'):
    from atpy.registry import register_set_reader #@UnresolvedImport
    from atpy.registry import register_reader #@UnresolvedImport
else:
    from atpy import register_set_reader #@UnresolvedImport @Reimport
    from atpy import register_reader #@UnresolvedImport @Reimport

import astro.ext_lib.starlight.starlighttable
import astro.ext_lib.starlight.starlighttablev4

register_set_reader('starlight', starlighttable.read_set)
register_set_reader('starlightv4', starlighttablev4.read_set)


