"""
desitarget.targetmask
=====================

This looks more like a script than an actual module.
"""
import os.path
import yaml
from pathlib import Path
from pkg_resources import resource_filename



# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
================
desiutil.bitmask
================

Mask bits for the spectro pipeline.

Individual packages will define their own mask bits and use this as a utility
access wrapper.  Typical users will get their bitmasks pre-made from those
packages, not from here.

Stephen Bailey, Lawrence Berkeley National Lab
Fall 2015

Examples
--------

desispec_ could create a ccdmask like this:

>>> from desiutil.bitmask import BitMask
>>> import yaml
>>> _bitdefs = yaml.safe_load('''
... ccdmask:
...     - [BAD,       0, "Pre-determined bad pixel (any reason)"]
...     - [HOT,       1, "Hot pixel"]
...     - [DEAD,      2, "Dead pixel"]
...     - [SATURATED, 3, "Saturated pixel from object"]
...     - [COSMIC,    4, "Cosmic ray"]
... ''')
...
>>> ccdmask = BitMask('ccdmask', _bitdefs)

Users would then access this mask with:

>>> from desispec.bitmasks import ccdmask
>>> ccdmask.COSMIC | ccdmask.SATURATED  #- 2**4 + 2**3
24
>>> ccdmask.mask('COSMIC')     # 2**4, same as ccdmask.COSMIC
16
>>> ccdmask.mask(4)            # 2**4, same as ccdmask.COSMIC
16
>>> ccdmask.COSMIC             # 2**4, same as ccdmask.mask('COSMIC')
16
>>> ccdmask.bitnum('COSMIC')
4
>>> ccdmask.bitname(4)
'COSMIC'
>>> ccdmask.names()
['BAD', 'HOT', 'DEAD', 'SATURATED', 'COSMIC']
>>> ccdmask.names(3)
['BAD', 'HOT']
>>> ccdmask.comment(0)
'Pre-determined bad pixel (any reason)'
>>> ccdmask.comment('COSMIC')
'Cosmic ray'


.. _desispec: http://desispec.readthedocs.io
"""


class _MaskBit(int):
    """A single mask bit.

    Subclasses :class:`int` to act like an :class:`int`, but allows the
    ability to extend with blat.name, blat.comment, blat.mask, blat.bitnum.

    Attributes
    ----------
    name : :class:`str`
        The name of the bit.
    bitnum : :class:`int`
        The number of the bit.  The value of the bit is ``2**bitnum``.
    mask : :class:`int`
        The value of the bit, ``2**bitnum``.
    comment : :class:`str`
        A comment explaining the meaning of the bit.
    """
    def __new__(cls, name, bitnum, comment, extra=dict()):
        self = super(_MaskBit, cls).__new__(cls, 2**bitnum)
        self.name = name
        self.bitnum = bitnum
        self.mask = 2**bitnum
        self.comment = comment
        self._extra = extra
        for key, value in extra.items():
            if hasattr(self, key):
                raise AttributeError(
                    "Bit {0} extra key '{1}' is already in use by int objects.".format(name, key))
            self.__dict__[key] = value
        return self

    def __str__(self):
        return ('{0.name:16s} bit {0.bitnum} mask 0x{0.mask:X} - ' +
                '{0.comment}').format(self)

    # def __repr__(self):
    #     return "_MaskBit(name='{0.name}', bitnum={0.bitnum:d}, comment='{0.comment}')".format(self)


#  Class to provide mask bit utility functions
class BitMask(object):
    """BitMask object to represent bit names, masks, and comments.

    Typical users are not expected to create BitMask objects directly;
    other packages like desispec and desitarget will have used this
    to pre-create the bitmasks for them using definition files in those
    packages.

    Parameters
    ----------
    name : :class:`str`
        Name of this mask, must be key in `bitdefs`.
    bitdefs : :class:`dict`
        Dictionary of different mask bit definitions;
        each value is a list of ``[bitname, bitnum, comment]``.
        A 4th entry is optional, which must be a dictionary.
    """

    def __init__(self, name, bitdefs):
        """Init.
        """
        self._bits = dict()
        self._name = name
        for x in bitdefs[name]:
            bitname, bitnum, comment = x[0:3]
            if len(x) == 4:
                extra = x[3]
                if not isinstance(extra, dict):
                    raise ValueError(
                        '{} extra values should be a dict'.format(bitname))
            else:
                extra = dict()
            self._bits[bitname] = _MaskBit(bitname, bitnum, comment, extra)
            self._bits[bitnum] = self._bits[bitname]

    def __getitem__(self, bitname):
        """Return mask for individual bitname.
        """
        return self._bits[bitname]

    def bitnum(self, bitname):
        """Return bit number (int) for this `bitname` (string).

        Parameters
        ----------
        bitname : :class:`str`
            The bit name.

        Returns
        -------
        :class:`int`
            The bit value.
        """
        return self._bits[bitname].bitnum

    def bitname(self, bitnum):
        """Return bit name (string) for this `bitnum` (integer).

        Parameters
        ----------
        bitnum : :class:`int`
            The number of the bit.

        Returns
        -------
        :class:`str`
            The name of the bit.
        """
        return self._bits[bitnum].name

    def comment(self, bitname_or_num):
        """Return comment for this bit name or bit number.

        Parameters
        ----------
        bitname_or_num : :class:`int` or :class:`str`
            Name of number of the mask.

        Returns
        -------
        :class:`str`
            The comment string.
        """
        return self._bits[bitname_or_num].comment

    def mask(self, name_or_num):
        """Return mask value.

        Parameters
        ----------
        name_or_num : :class:`int` or :class:`str`
            Name of number of the mask.

        Returns
        -------
        :class:`int`
            The value of the mask.

        Examples
        --------

        >>> bitmask.mask(3)         #  2**3
        8
        >>> bitmask.mask('BLAT')
        >>> bitmask.mask('BLAT|FOO')
        """
        if isinstance(name_or_num, int):
            return self._bits[name_or_num].mask
        else:
            mask = 0
            for name in name_or_num.split('|'):
                mask |= self._bits[name].mask
            return mask

    def names(self, mask=None):
        """Return list of names of masked bits.

        Parameters
        ----------
        mask : :class:`int`, optional
            The mask integer to convert to names. If not supplied,
            return names of all known bits.

        Returns
        -------
        :class:`list`
            The list of names contained in the mask.
        """
        names = list()
        if mask is None:
            # return names in sorted order of bitnum
            bitnums = [x for x in self._bits.keys() if isinstance(x, int)]
            for bitnum in sorted(bitnums):
                names.append(self._bits[bitnum].name)
        else:
            mask = int(mask)  # workaround numpy issue #2955 for uint64
            bitnum = 0
            while 2**bitnum <= mask:
                if (2**bitnum & mask):
                    if bitnum in self._bits.keys():
                        names.append(self._bits[bitnum].name)
                    else:
                        names.append('UNKNOWN' + str(bitnum))
                bitnum += 1

        return names

    def __getattr__(self, name):
        """Enable ``mask.BITNAME`` equivalent to ``mask['BITNAME']``.
        """
        if name in self._bits:
            return self._bits[name]
        else:
            raise AttributeError('Unknown mask bit name ' + name)

    def __repr__(self):
        '''Return yaml representation defining the bits of this bitmask.
        '''
        result = list()
        result.append(self._name + ':')
        # return names in sorted order of bitnum
        bitnums = [x for x in self._bits.keys() if isinstance(x, int)]
        for bitnum in sorted(bitnums):
            bit = self._bits[bitnum]
            # format the line for single bit, with or without extra keys
            line = '  - [{:16s} {:2d}, "{}"'.format(
                bit.name+',', bit.bitnum, bit.comment)
            if len(bit._extra) > 0:
                line = line + ', '+str(bit._extra)+']'
            else:
                line = line + ']'
            result.append(line)

        return "\n".join(result)


def load_mask_bits(prefix=""):
    """Load bit definitions from yaml file.
    """
    us = ""
    if len(prefix) > 0:
        us = '_'
    prename = prefix+us
    fn = os.path.join(prefix, "data", "{}targetmask.yaml".format(prename))
    # _filepath = resource_filename('desitarget', fn)
    _filepath = Path(fn)
    with open(_filepath) as fx:
        bitdefs = yaml.safe_load(fx)
        try:
            bitdefs = _load_mask_priorities(bitdefs, handle="priorities", prename=prename)
        except TypeError:
            pass
        try:
            bitdefs = _load_mask_priorities(bitdefs, handle="numobs", prename=prename)
        except TypeError:
            pass
    return bitdefs


def _load_mask_priorities(bitdefs, handle="priorities", prename=""):
    """Priorities and NUMOBS are defined in the yaml file, but they aren't
    a bitmask and so require some extra processing.
    """
    for maskname, priorities in bitdefs[handle].items():
        for bitname in priorities:
            # -"SAME_AS_XXX" enables one bit to inherit priorities from another
            if isinstance(priorities[bitname], str) and priorities[bitname].startswith('SAME_AS_'):
                other = priorities[bitname][8:]
                priorities[bitname] = priorities[other]

            # -fill in default "more" priority to be same as "unobs"
            # ADM specifically applies to dictionary of priorities
            if handle == 'priorities':
                if isinstance(priorities[bitname], dict):
                    if 'MORE_ZWARN' not in priorities[bitname]:
                        priorities[bitname]['MORE_ZWARN'] = priorities[bitname]['UNOBS']
                    if 'MORE_ZGOOD' not in priorities[bitname]:
                        priorities[bitname]['MORE_ZGOOD'] = priorities[bitname]['UNOBS']

                    # - fill in other states as priority=1
                    for state, blat, foo in bitdefs[prename+'obsmask']:
                        if state not in priorities[bitname]:
                            priorities[bitname][state] = 1
                else:
                    priorities[bitname] = dict()

        # - add to the extra info dictionary for this target mask
        for bitdef in bitdefs[maskname]:
            bitname = bitdef[0]
            bitdef[3][handle] = priorities[bitname]
    return bitdefs


# -convert to BitMask objects
# if bitdefs is None:
#     load_bits()
_bitdefs = load_mask_bits()
try:
    desi_mask = BitMask('desi_mask', _bitdefs)
    mws_mask = BitMask('mws_mask', _bitdefs)
    bgs_mask = BitMask('bgs_mask', _bitdefs)
    scnd_mask = BitMask('scnd_mask', _bitdefs)
    obsconditions = BitMask('obsconditions', _bitdefs)
    obsmask = BitMask('obsmask', _bitdefs)
    targetid_mask = BitMask('targetid_mask', _bitdefs)
    zwarn_mask = BitMask('zwarn_mask', _bitdefs)
except TypeError:
    desi_mask = object()
    mws_mask = object()
    bgs_mask = object()
    scnd_mask = object()
    obsconditions = object()
    obsmask = object()
    targetid_mask = object()
    zwarn_mask = object()


