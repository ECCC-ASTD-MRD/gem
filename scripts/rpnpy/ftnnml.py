#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author: Stephane Chamberland <stephane.chamberland@canada.ca>
# Copyright: LGPL 2.1
"""
Module ftnnml contains the classes used to manipulate fortran namelist files,
extract/get values, get list of namelists and vars, beautify,
set/rename vat values, create new namelist, delete values or namelists, ...
"""

## Profiling: python -m cProfile -- ./bin/getnml2 -l -f phydict.nml

## TODO: option to switch derived type form to list form (vice versa)
## TODO: trim too long list...  with repeated values?
##       take into account fixed length string list (no commas or quote when write)
## TODO: catch value format like: kfctrig4=2*0.0000000E+00  , 2*5.0000001E-02
## TODO: Make Logicals format consistent (F/false/.f./.false.)

import re
import sys

_islisttype   = lambda x: isinstance(x, (list, tuple))
_isstringtype = lambda x: isinstance(x, str)
_cleanName    = lambda x: x.lower().replace('\n', ' ').strip()

class FtnNmlObj(object):
    """
    Fortran Namlist container base class (need to be subclassed)
    """
    allowedSubClass = None
    (splitPattern, matchStartPattern, matchEndPattern)  = (None, None, None)

    @classmethod
    def encode(cls, data):
        """Encode a strng before parsing"""
        return data

    @classmethod
    def decode(cls, data):
        """Decode what have been encoded with encode"""
        return data

    @classmethod
    def _parseSubContent(cls, mystr):
        """Deffer parsing to allowedSubClass"""
        mystr2 = cls.decode(mystr)
        return cls.allowedSubClass.parseToList(mystr2 \
            if cls.allowedSubClass \
            else mystr2)

    @classmethod
    def parseToList(cls, mystr):
        """Parse string"""
        mystr2 = cls.encode(mystr)
        myitemList = [mystr2]
        if cls.splitPattern:
            myitemList = [item.replace("\x01", '\n')
                          for item in re.split(cls.splitPattern, mystr2)]
        (myitem, mySubClassObj, listdata) = ('', None, [])
        for myline in myitemList:
            (m0, m1) = (None, None)
            if cls.matchStartPattern:
                m0 = re.match(cls.matchStartPattern, myline)
            if m0:
                if mySubClassObj:
                    mySubClassObj.set(cls._parseSubContent(myitem))
                    listdata.append(mySubClassObj)
                elif len(myitem) > 0:
                    listdata.append(cls.decode(myitem))
                (myitem, mySubClassObj) = ('', cls(''))
                for mykey in m0.groupdict().keys():
                    if mykey in ('strStart', 'name'):
                        mySubClassObj.rename(m0.group(mykey))
                    else:
                        mySubClassObj.prop[mykey] = m0.group(mykey)
            else:
                if cls.matchEndPattern:
                    m1 = re.match(cls.matchEndPattern, myline)
                if m1 and mySubClassObj:
                    junk = ''
                    for mykey in m1.groupdict().keys():
                        if    mykey == 'data': myitem += m1.group(mykey)
                        elif  mykey == 'junk': junk   += m1.group(mykey)
                        else: mySubClassObj.prop[mykey] = m1.group(mykey)
                    mySubClassObj.set(cls._parseSubContent(myitem))
                    listdata.append(mySubClassObj)
                    if junk: listdata.append(cls.decode(junk))
                    (myitem, mySubClassObj) = ('', None)
                else:
                    myitem += myline
        if mySubClassObj:
            mySubClassObj.set(cls._parseSubContent(myitem))
            listdata.append(mySubClassObj)
        elif len(myitem) > 0:
            listdata.append(cls.decode(myitem))
        return listdata

    @classmethod
    def prepStr(cls, mystr, clean=False, uplowcase=None):
        changeCase = lambda x: x
        if uplowcase:
            changeCase = lambda x: x.lower()
            if _isstringtype(uplowcase) and uplowcase[0].lower() == 'u':
                changeCase = lambda x: x.upper()
        return changeCase(mystr.lstrip() if clean else mystr)

    def __init__(self, name):
        (self.name, self.data) = (_cleanName(name), [])
        self.prop = { #Start SepS Sep1 data End SepE
            'strStart' : name,
            'strSepS'  : '',
            'strSep1'  : '',
            'strEnd'   : '',
            'strSepE'  : '',
           }

    def rename(self, name):
        """Properly set name of object (please avoid obj.name = 'name')"""
        (self.name, self.prop['strStart']) = (_cleanName(name), name)

    def get(self, name=None):
        """Return sub object with given name
           return the list of contained objects otherwise"""
        if not name: return self.data
        try:
            return self.data[self.keyIndex(name)]
        except:
            ## sys.stderr.write('Known Keys:'+repr(self.keys()))
            raise KeyError(" ({0}) Oops! get, Key not found: {1}"\
                           .format(self.__class__.__name__, name))

    def set(self, namedata, data=None):
        """Set sub object data with given name
           replace the list of contained objects otherwise"""
        name = (namedata if data else None)
        if not data: data = namedata
        if name: self.get(name).set(data)
        else:    self.data = (data if _islisttype(data) else [data])

    def add(self, data=None):
        if not isinstance(data, self.allowedSubClass):
            raise TypeError(" ({0}) Oops! add, provided data is of wrong type: {1} (accepting: {2})"\
                            .format(self.__class__.__name__, str(type(data)),
                                    str(self.allowedSubClass)))
        if _cleanName(data.name) in self.keys():
            raise KeyError(" ({0}) Oops! add, Key already exists: {1}"\
                           .format(self.__class__.__name__, data.name))
        self.data.append(data)

    def rm(self, name=None):
        """Delete sub object with given name
           delete all contained objects otherwise"""
        if not name:
            self.data = []
        else:
            try:
                del self.data[self.keyIndex(name)]
            except:
                raise KeyError(" ({0}) Oops! rm, Key not found: {1}"\
                               .format(self.__class__.__name__, name))

    def keyIndex(self, name):
        (name2, myindex) = (_cleanName(name), -1)
        for item in self.data:
            myindex += 1
            if isinstance(item, FtnNmlObj) and item.name == name2:
                return myindex
        raise KeyError(" ({0}) Oops! keyIndex, Key not found: {1}"\
                       .format(self.__class__.__name__, name))

    def keys(self):
        """Return list of keys (contained objects names)"""
        return [item.name for item in self.data
                if isinstance(item, FtnNmlObj) and item.name]

    def __repr__(self):
        return "{0}({1},{2},{3},{4},d={5},{6},{7})"\
                .format(self.__class__.__name__, repr(self.name),
                        repr(self.prop['strStart']),
                        repr(self.prop['strSepS']),
                        repr(self.prop['strSep1']),
                        repr(self.data),
                        repr(self.prop['strEnd']),
                        repr(self.prop['strSepE'])
                    )

    def __str__(self):
        return self.toStr()

    def startStr(self, clean=False, uplowcase=None):
        return self.prepStr(self.prop['strStart']+self.prop['strSepS'],
                            clean, uplowcase)
    def sepStr(self, clean=False, uplowcase=None):
        return self.prepStr(self.prop['strSep1'], clean, uplowcase)
    def endStr(self, clean=False, uplowcase=None):
        return self.prepStr(self.prop['strEnd']+self.prop['strSepE'],
                            clean, uplowcase)

    def _myToStr(self, data, clean=False, uplowcase=None, updnsort=None):
        if isinstance(data, FtnNmlObj):
            return data.toStr(clean, uplowcase, updnsort)
        else:
            return self.prepStr(str(data), clean, uplowcase)

    def toStr(self, clean=False, uplowcase=None, updnsort=None):
        """Return String representation of the FtnNml Object, recursively"""
        if updnsort: clean=True
        data = (self.data if _islisttype(self.data) else [self.data])
        if clean: data = [s for s in data if type(s) != type(' ')]
        if updnsort:
            datakeys = []
            ii = 0
            for item in data:
                datakeys.append((item.name.lower(), ii))
                ii += 1
            datakeys.sort()
            data2 = []
            for (nn, ii) in datakeys:
                data2.append(data[ii])
            data = data2
        return self.startStr(clean, uplowcase) \
            + self.sepStr(clean, uplowcase) \
            + ''.join([self._myToStr(s, clean, uplowcase, updnsort)
                       for s in data]) \
            + self.endStr(clean, uplowcase)\


class FtnNmlVal(FtnNmlObj):
    """
    Fortran Namlist value container
    """

    @classmethod
    def parseToList(cls, mystr):
        if not mystr: return ['']
        m = re.match(r'^([\s\t\n]*)([\w\W]+?)([\s\t\n,]*)$', mystr, re.I)
        if m:
            return [s for s in [str(m.group(1)), cls(str(m.group(2))),
                                str(m.group(3))] if s]
        else:
            return [cls(str(mystr))]

    def __init__(self, value):
        FtnNmlObj.__init__(self, 'v')
        self.data = value

    def __repr__(self):
        return "{0}(d={1})".format(self.__class__.__name__, repr(self.data))

    def toStr(self, clean=False, uplowcase=None, updnsort=None):
        data = (self.data[0] if _islisttype(self.data) else self.data)
        if clean: return str(data).replace('\n', '')
        else:     return str(data)

    def rename(self, name):
        pass


class FtnNmlKeyVal(FtnNmlObj):
    """
    Fortran Namlist key=value container
    """

    allowedSubClass   = FtnNmlVal
    # ftnVarPattern     = r'\w+'
    # ftnVarPattern     = r'(?:\w+|\w+%\w+)'
    ftnVarPattern     = r'(?:\w+|\w+%\w+|\w+\([0-9]+\)|\w+\([0-9]+\)%\w+)'
    splitPattern      = re.compile(r'([\s\t]*'+ftnVarPattern+r'[\s\t]*=[\s\t]*)')
    matchStartPattern = re.compile(r'^(?P<strStart>[\s\t]*'+ftnVarPattern+r'[\s\t]*)(?P<strSep1>=[\s\t]*)$')
    matchEndPattern   = re.compile(r'^(?P<data>[\w\W]+)(?P<strEnd>[,\n])(?P<junk>[,\n\t\s]*)$')

    @classmethod
    def encode(cls, data):
        return re.sub(r'("[^"]+?"|\'[^\']+?\')',
                      lambda m: m.group(0).replace("=", "\x00"),
                      data.replace('\n', "\x01"))

    @classmethod
    def decode(cls, data):
        return data.replace("\x00", "=").replace("\x01", '\n')

    def __init__(self, name, data=None):
        FtnNmlObj.__init__(self, name)
        self.prop['strSep1'] = '='
        self.prop['strEnd']  = ',\n'
        if data: self.set(data)

    def startStr(self, clean=False, uplowcase=None):
        if clean: return self.prepStr(self.name, clean, uplowcase)
        else:     return FtnNmlObj.startStr(self, clean, uplowcase)

    def sepStr(self, clean=False, uplowcase=None):
        if clean: return '='
        else:     return FtnNmlObj.sepStr(self, clean, uplowcase)

    def endStr(self, clean=False, uplowcase=None):
        if clean: return '\n'
        else:     return FtnNmlObj.endStr(self, clean, uplowcase)


class FtnNmlSection(FtnNmlObj):
    """
    Fortran Namlist 'namelist' container
    """

    allowedSubClass = FtnNmlKeyVal
    splitPattern      = re.compile(r'([^\n]*\n)')
    matchStartPattern = re.compile(r'^(?P<strStart>[\s\t]*&[^\s\t]+[\s\t]*)(?P<strSepS>\n)$')
    matchEndPattern   = re.compile(r'^(?P<strEnd>[\s\t]*/[\s\t]*)(?P<strSepE>\n)$')

    def __init__(self, name):
        FtnNmlObj.__init__(self, name)
        self.prop['strStart'] = '&'+name
        self.prop['strSepS']  = '\n'
        self.prop['strEnd']   = '/'
        self.prop['strSepE']  = '\n'

    def rename(self, name):
        ## if re.match(self.matchStartPattern,name): #TODO
        if re.match('(&|[^&\w]+&)', name):
            self.prop['strStart'] = name
        else:
            self.prop['strStart'] = '&'+name.lstrip().replace('\n', '')
        self.name = _cleanName(self.prop['strStart'].replace('&', ''))


class FtnNmlFile(FtnNmlObj):
    """
    Fortran Namlist file container
    """

    allowedSubClass = FtnNmlSection

    @classmethod
    def parseToList(cls, mystr):
        return cls._parseSubContent(mystr)

    def __init__(self, name, fromFile=True):
        FtnNmlObj.__init__(self, name)
        self.prop['strStart'] = ''
        if fromFile: self.read(name)

    def parse(self, mystr):
        self.data = self.__class__.parseToList(mystr)

    def read(self, filename):
        """Read and parse file"""
        rawdata = ""
        try:
            fd = open(filename, "rb")
            try:
                if sys.stdin.encoding:
                        rawdata = "".join([x.decode(sys.stdin.encoding) for x in fd.readlines()])
                else:
                        rawdata = "".join(fd.readlines())
            finally:
                fd.close()
        except IOError:
            raise IOError(" Oops! File does not exist or is not readable: {0}".
                          format(filename))
        self.parse(rawdata.replace("\x00", ""))

    def write(self, filename, clean=False, uplowcase=None, updnsort=None):
        """Write nml to file"""
        try:
            fd = open(filename, "wb")
            try:
                fd.write(self.toStr(clean, uplowcase, updnsort).encode('ascii'))
            except IOError:
                raise IOError(" Oops! Cannot wrtie to file: {0}".
                              format(filename))
            finally:
                fd.close()
        except IOError:
            raise IOError(" Oops! Cannot open file: {0}".format(filename))



if __name__ == '__main__':
    import doctest
    doctest.testmod()
    #TODO: base class on list or dict
    #TODO: itf consistenncy (name,data) in set,get,rename,rm,add,__init__
    #TODO: data should be converted to list if not already
    #TODO: implement sort
    #TODO: doctests

    ## import pprint
    ## pp = pprint.PrettyPrinter(indent=4)
    ## filename = 'gem_settings.nml'
    ## b = FtnNmlFile(filename)
    ## print b
    ## print repr(b)

    ## #T: get
    ## print '---- List of nml'
    ## mynmls = b.keys()
    ## print filename+': ',', '.join(mynmls)
    ## print '---- List of var per nml'
    ## for nmlkey in mynmls:
    ##     nml = b.get(nmlkey)
    ##     mykeys = nml.keys()
    ##     print '&'+nml.name+': '+', '.join(mykeys)
    ## print '---- List of values'
    ## for nmlkey in mynmls:
    ##     nml = b.get(nmlkey)
    ##     mykeys = nml.keys()
    ##     for varkey in mykeys:
    ##         kv = nml.get(varkey)
    ##         val = kv.get('v')
    ##         #valstr = str(val) #equivalent to: val.toStr()
    ##         valstr = val.toStr(clean=True)
    ##         print '&'+nml.name+'/'+varkey+'='+valstr

    ## #T: set
    ## print '---- Change value'
    ## val = b.get('gem_cfgs').get('lctl_debug_l').get('v')
    ## print '&gem_cfgs/lctl_debug_l= '+str(val)
    ## val.set('.T.')
    ## print '&gem_cfgs/lctl_debug_l= '+str(val)
    ## print '---- Change var name'
    ## kv = b.get('gem_cfgs').get('hyb')
    ## print [kv.name,kv.prop['strStart']]
    ## kv.rename('levels')
    ## print [kv.name,kv.prop['strStart']]
    ## print '---- Change nml name'
    ## nml = b.get('gem_cfgs')
    ## print [nml.name,nml.prop['strStart']]
    ## nml.rename('sps_cfgs')
    ## print [nml.name,nml.prop['strStart']]

    ## #T: del
    ## print '---- Delete nml var'
    ## nml = b.get('sps_cfgs')
    ## print '&'+nml.name+': '+', '.join(nml.keys())
    ## nml.rm('sol_type2_s')
    ## nml.rm('etiket')
    ## print '&'+nml.name+': '+', '.join(nml.keys())
    ## print '---- Delete nml'
    ## print filename+': ',', '.join(b.keys())
    ## b.rm('grid_gu')
    ## print filename+': ',', '.join(b.keys())

    ## #T: add
    ## print repr(b)
    ## print '---- Add nml var'
    ## nml = b.get('sps_cfgs')
    ## print '&'+nml.name+': '+', '.join(nml.keys())
    ## nml.add(FtnNmlKeyVal('newvar',FtnNmlVal(1)))
    ## nml.add(FtnNmlKeyVal('n2',FtnNmlVal(4)))
    ## print '&'+nml.name+': '+', '.join(nml.keys())
    ## print '---- Add nml'
    ## print filename+': ',', '.join(b.keys())
    ## b.add(FtnNmlSection('mytoto'))
    ## print filename+': ',', '.join(b.keys())
    ## nml = b.get('mytoto')
    ## nml.add(FtnNmlKeyVal('totavar',FtnNmlVal('w')))

    ## print '----'
    ## print b
    ## print b.toStr(clean=True)
    ## print repr(b)

    ## ## ## import doctest
    ## ## ## doctest.testmod()
    ## ## ## verbose = 0
    ## ## ## a = FtnNmlFile('gem_settings.nml')


# -*- Mode: C; tab-width: 4; indent-tabs-mode: nil -*-
# vim: set expandtab ts=4 sw=4:
# kate: space-indent on; indent-mode cstyle; indent-width 4; mixedindent off;
