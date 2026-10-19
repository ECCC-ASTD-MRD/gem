#!/usr/bin/env python
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
    Fortran Namelist container base class (need to be subclassed)
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
        ## print(cls.__name__," tolist:",mystr2,":") 
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
        (self._name, self._data) = (_cleanName(name), [])
        self.prop = { #Start SepS Sep1 data End SepE
            'strStart' : name,
            'strSepS'  : '',
            'strSep1'  : '',
            'strEnd'   : '',
            'strSepE'  : '',
           }

    def rename(self, name):
        """Properly set name of object (please avoid obj._name = 'name')"""
        (self._name, self.prop['strStart']) = (_cleanName(name), name)

    def get(self, name=None):
        """Return sub object with given name
           return the list of contained objects otherwise"""
        if not name: return self._data
        try:
            return self._data[self.keyIndex(name)]
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
        else:    self._data = (data if _islisttype(data) else [data])

    def add(self, data=None):
        if not isinstance(data, self.allowedSubClass):
            raise TypeError(" ({0}) Oops! add, provided data is of wrong type: {1} (accepting: {2})"\
                            .format(self.__class__.__name__, str(type(data)),
                                    str(self.allowedSubClass)))
        if _cleanName(data._name) in self.keys():
            raise KeyError(" ({0}) Oops! add, Key already exists: {1}"\
                           .format(self.__class__.__name__, data._name))
        self._data.append(data)

    def rm(self, name=None):
        """Delete sub object with given name
           delete all contained objects otherwise"""
        if not name:
            self._data = []
        else:
            try:
                del self._data[self.keyIndex(name)]
            except:
                raise KeyError(" ({0}) Oops! rm, Key not found: {1}"\
                               .format(self.__class__.__name__, name))

    def keyIndex(self, name):
        (name2, myindex) = (_cleanName(name), -1)
        for item in self._data:
            myindex += 1
            if isinstance(item, FtnNmlObj) and item._name == name2:
                return myindex
        raise KeyError(" ({0}) Oops! keyIndex, Key not found: {1}"\
                       .format(self.__class__.__name__, name))

    def keys(self):
        """Return list of keys (contained objects names)"""
        return [item._name for item in self._data
                if isinstance(item, FtnNmlObj) and item._name]

    def __repr__(self):
        return "{0}({1},{2},{3},{4},d={5},{6},{7})"\
                .format(self.__class__.__name__, repr(self._name),
                        repr(self.prop['strStart']),
                        repr(self.prop['strSepS']),
                        repr(self.prop['strSep1']),
                        repr(self._data),
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
        # sys.stderr.write("{}._myToStr({}, {})\n".format(self.__class__, repr(data), clean))
        if isinstance(data, FtnNmlObj):
            return data.toStr(clean, uplowcase, updnsort)
        else:
            return self.prepStr(str(data), clean, uplowcase)

    def toStr(self, clean=False, uplowcase=None, updnsort=None):
        """Return String representation of the FtnNml Object, recursively"""
        if updnsort: clean=True
        data = (self._data if _islisttype(self._data) else [self._data])
        if clean: data = [s for s in data if type(s) != type(' ')]
        if updnsort:
            datakeys = []
            ii = 0
            for item in data:
                datakeys.append((item._name.lower(), ii))
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
    Fortran Namelist value container
    """
            
    _re_encode = re.compile(r'("[^"]+?"|\'[^\']+?\'|""|\'\')')
    _re_split1 = re.compile(r'([\s\t]*![^\n]*\n|[\s\t]*![^\n]*$)')
    _re_split2 = re.compile(r'([\s\t\n]*,[\s\t\n]*)')

    _SPLIT1  = lambda cls, s: re.split(cls._re_split1, cls.encode(s))  #Note: minimize calling this as it is very heavy
    _FILTER1 = lambda cls, s: ("comment", cls.decode(s)) if (len(s.lstrip()) > 0 and s.lstrip()[0] == '!') else s
    _SPLIT2  = lambda cls, s: re.split(cls._re_split2, s)
    _FILTER2 = lambda cls, s: ("sep", cls.decode(s)) if (len(s.lstrip()) > 0 and s.lstrip()[0] == ',') else cls.decode(s)
    _FILTER3 = lambda cls, s: [cls._FILTER2(cls, s) for s in cls._SPLIT2(cls, s) if len(s) > 0] if isinstance(s, str) else [s]

    @classmethod
    def encode(cls, data):
        return re.sub(cls._re_encode,
                      lambda m: m.group(0).replace("!", "\x00").replace(",", "\x01"),
                      data)

    @classmethod
    def decode(cls, data):
        return data.replace("\x00", "!").replace("\x01", ',')

    @classmethod
    def parseToList(cls, mystr):
        if not mystr:
            return ['']
        mystr = cls.encode(mystr)
        if '!' in mystr:
            mylist = [cls._FILTER1(cls, s) for s in cls._SPLIT1(cls, mystr)]
        else:
            mylist = (cls._FILTER1(cls, mystr), )
        listdata = [x for y in mylist for x in cls._FILTER3(cls, y)]
        return cls(listdata)

    def get(self, name=None, unquote=False):
        """Return sub object with given name
           return the list of contained objects otherwise"""
        #if not name: return self._data
        ## sys.stderr.write("{}.get({})\n".format(self.__class__, name))
        return self.toStr(clean=True, unquote=unquote)

    def __init__(self, value):
        FtnNmlObj.__init__(self, 'v')
        self._data = value

    def __repr__(self):
        return "{0}(d={1})".format(self.__class__.__name__, repr(self._data))

    def toStr(self, clean=False, uplowcase=None, updnsort=None, unquote=False):
        ## sys.stderr.write("{}.toStr({},{},{})\n".format(self.__class__, clean,uplowcase,updnsort))
        # if updnsort: clean=False
        # print(repr(self._data))
        data = (self._data if _islisttype(self._data) else [self._data])
        data = [s if _islisttype(s) else ('',s) for s in data]
        ## print(repr(data))
        if clean:
            data = [s for s in data if s[0] != 'comment']
            if data[-1][0] == 'sep':
                del data[-1]
            ## sys.stderr.write("FtnNmlVal.toStr(): data={}\n".format(repr(data)))
            # print(1,repr(data))
            data = [s[1] if isinstance(s, (list, tuple)) else s for s in data]
            # print(2,repr(data))
            data = [s.strip() for s in data]
            # print(3,repr(data))
            if unquote:
                data = [s.strip(s[0]) if s[0] in ("'", '"') else s for s in data]
        else:
            #Note: doesn't make sens to "unquote" if not clean. should we warn?
            data = [s[1] for s in data]
        try:
           return ''.join(data)
        except:
            print(repr(data))
            raise
 
    def rename(self, name):
        pass


class FtnNmlKeyVal(FtnNmlObj):
    """
    Fortran Namelist key=value container
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
        s = re.sub(r'("[^"]+?"|\'[^\']+?\'|""|\'\')',
                   lambda m: m.group(0).replace("=", "\x00").replace("!", "\x02"),
                   data.replace('\n', "\x01"))
        s = re.sub(r'([\s\t]*![^\x01]*\x01)',
                   lambda m: m.group(0).replace("=", "\x00"),
                   s)
        return s

    @classmethod
    def decode(cls, data):
        return data.replace("\x00", "=").replace("\x02", '!').replace("\x01", '\n')

    def __init__(self, name, data=None):
        FtnNmlObj.__init__(self, name)
        self.prop['strSep1'] = '='
        self.prop['strEnd']  = ',\n'
        if data: self.set(data)

    def startStr(self, clean=False, uplowcase=None):
        if clean: return self.prepStr(self._name, clean, uplowcase)
        else:     return FtnNmlObj.startStr(self, clean, uplowcase)

    def sepStr(self, clean=False, uplowcase=None):
        if clean: return '='
        else:     return FtnNmlObj.sepStr(self, clean, uplowcase)

    def endStr(self, clean=False, uplowcase=None):
        if clean: return '\n'
        else:     return FtnNmlObj.endStr(self, clean, uplowcase)


class FtnNmlSection(FtnNmlObj):
    """
    Fortran Namelist 'namelist' container
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
        self._name = _cleanName(self.prop['strStart'].replace('&', ''))


class FtnNmlFile(FtnNmlObj):
    """
    Fortran Namelist file container
    """

    allowedSubClass = FtnNmlSection

    @classmethod
    def parseToList(cls, mystr):
        return cls._parseSubContent(mystr)

    def __init__(self, name, fromFile=True):
        FtnNmlObj.__init__(self, name)
        self.prop['strStart'] = ''
        if fromFile: self.read(name)

    def printtmp(self, myobj):
        for item in myobj:
            if isinstance(item, FtnNmlObj):
                print("1 ", repr(item))
                if isinstance(item._data, (list, FtnNmlObj)):
                    self.printtmp(item._data)
                ## else:
                ##     # print("2 ",repr(item._data),repr(item))
                ##     print("2 ", repr(item))
            ## else:
            ##     print("3 ", repr(item))
        
    def parse(self, mystr):
        self._data = self.__class__.parseToList(mystr)
        # self.printtmp(self._data)
 
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
    ##     print '&'+nml._name+': '+', '.join(mykeys)
    ## print '---- List of values'
    ## for nmlkey in mynmls:
    ##     nml = b.get(nmlkey)
    ##     mykeys = nml.keys()
    ##     for varkey in mykeys:
    ##         kv = nml.get(varkey)
    ##         val = kv.get('v')
    ##         #valstr = str(val) #equivalent to: val.toStr()
    ##         valstr = val.toStr(clean=True)
    ##         print '&'+nml._name+'/'+varkey+'='+valstr

    ## #T: set
    ## print '---- Change value'
    ## val = b.get('gem_cfgs').get('lctl_debug_l').get('v')
    ## print '&gem_cfgs/lctl_debug_l= '+str(val)
    ## val.set('.T.')
    ## print '&gem_cfgs/lctl_debug_l= '+str(val)
    ## print '---- Change var name'
    ## kv = b.get('gem_cfgs').get('hyb')
    ## print [kv._name,kv.prop['strStart']]
    ## kv.rename('levels')
    ## print [kv._name,kv.prop['strStart']]
    ## print '---- Change nml name'
    ## nml = b.get('gem_cfgs')
    ## print [nml._name,nml.prop['strStart']]
    ## nml.rename('sps_cfgs')
    ## print [nml._name,nml.prop['strStart']]

    ## #T: del
    ## print '---- Delete nml var'
    ## nml = b.get('sps_cfgs')
    ## print '&'+nml._name+': '+', '.join(nml.keys())
    ## nml.rm('sol_type2_s')
    ## nml.rm('etiket')
    ## print '&'+nml._name+': '+', '.join(nml.keys())
    ## print '---- Delete nml'
    ## print filename+': ',', '.join(b.keys())
    ## b.rm('grid_gu')
    ## print filename+': ',', '.join(b.keys())

    ## #T: add
    ## print repr(b)
    ## print '---- Add nml var'
    ## nml = b.get('sps_cfgs')
    ## print '&'+nml._name+': '+', '.join(nml.keys())
    ## nml.add(FtnNmlKeyVal('newvar',FtnNmlVal(1)))
    ## nml.add(FtnNmlKeyVal('n2',FtnNmlVal(4)))
    ## print '&'+nml._name+': '+', '.join(nml.keys())
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
