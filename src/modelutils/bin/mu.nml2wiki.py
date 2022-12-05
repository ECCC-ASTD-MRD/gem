#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Author: Stephane Chamberland <stephane.chamberland@canada.ca>
# Copyright: LGPL 2.1
"""
"""

import re
import sys

def nml2wiki(nmlfile,nmlname,wikifile):
    print("Processing nml {0} in file {1} to {2}".format(nmlname,nmlfile,wikifile))
    try:
        fd = open(nmlfile,"rb")
        try:     rawdata = "".join([x.strip()+'\n' for x in fd.readlines()])
        finally: fd.close()
    except IOError:
        raise IOError(" Oops! File does not exist or is not readable: {0}".format(nmlfile))
    try:
        fd = open(wikifile,"wb")
    except IOError:
        raise IOError(" Oops! Unable to open file for writing: {0}".format(wikifile))
    data = [x for x in re.split('\n\n+',rawdata) if re.search(r'\bnamelist[\s\t]*[/][\s\t]*'+nmlname+'[\s\t]*[/]',x,re.I)]
    fd.write("""
=== {0} Namelist ===

{{| class="wikitable"
|-
! style="width: 10em;" | Name 
! style="width: 40em;" | Description 
! style="width: 10em;" | Default Value 
! Type
""".format(nmlname.upper()))
    for x in data:
        a = re.search(r'\bnamelist[\s\t]*[/][\s\t]*'+nmlname+'[\s\t]*[/][\s\t]*([a-zA-Z][a-zA-Z0-9_]+)',x,re.I)
        if not a:
            sys.stderr.write("Problem1 with: \n"+repr(x)+'\n\n')
            continue
        varname = a.groups()[0]
        #a = re.search(r'([^\s\n]+)[\s\t]*::[\s\t]*'+varname+'[\s\t]*=[\s\t]*([^\n]+)',x,re.I)
        a = re.search(r'([^\n]+)[\s\t]*::[\s\t]*'+varname+'[\s\t]*=[\s\t]*([^\n]+)',x,re.I)
        if not a:
            #a = re.search(r'([^\s\n]+)[\s\t]*::[\s\t]*'+varname+'[(0-9)]+[\s\t]*=[\s\t]*([^\n]+)',x,re.I)
            a = re.search(r'([^\n]+)[\s\t]*::[\s\t]*'+varname+'[(0-9)]+[\s\t]*=[\s\t]*([^\n]+)',x,re.I)
        if not a:
            #a = re.search(r'([^\s\n]+)[\s\t]*::[\s\t]*'+varname+'[(][a-zA-Z0-9_+-=]+[)][\s\t]*=[\s\t]*([^\n]+)',x,re.I)
            a = re.search(r'([^\n]+)[\s\t]*::[\s\t]*'+varname+'[(][a-zA-Z0-9_+-=]+[)][\s\t]*=[\s\t]*([^\n]+)',x,re.I)
        if not a:
            sys.stderr.write("Problem2 with: \n"+repr(x)+'\n\n')
            continue
        vartype = a.groups()[0]
        varval  = a.groups()[1]
        vardesc = "\n".join([re.sub(r'^[\s\t]*[!][#][\s\t]?','',x1) for x1 in x.split("\n") if re.match(r'^[\s\t]*[!][#]',x1)])

        fd.write("""
|-
| {0} ||
{1}
| {2} || {3}\n""".format(varname,vardesc,varval,vartype))
    fd.write('|}')
    fd.close()
    
if __name__ == '__main__':
    a = (
        ('grid_options.F90', 'grid', 'grid.wiki.txt'),
        ('grdc_options.F90', 'grdc', 'grdc.wiki.txt'),
        ('step_options.F90', 'step', 'step.wiki.txt'),
        ('adw_options.F90', 'adw_cfgs', 'adw.wiki.txt'),
        ('adv_options.F90', 'adv_cfgs', 'adv.wiki.txt'),
        ('ens_options.F90', 'ensembles', 'ens.wiki.txt'),
        ('gem_options.F90', 'gem_cfgs', 'gem.wiki.txt'),
        ('phy_options.F90', 'physics_cfgs', 'physics_cfgs.wiki.txt'),
        ('cnv_options.F90', 'convection_cfgs', 'convection_cfgs.wiki.txt'),
        ('sfc_options.F90', 'surface_cfgs', 'surface_cfgs.wiki.txt'),
    ##     ('./src/series/series.cdk', 'series', 'series.wiki.txt'),
        )
    for (nf,nm,wf) in a:
        nml2wiki(nf,nm,wf)
        
# -*- Mode: C; tab-width: 4; indent-tabs-mode: nil -*-
# vim: set expandtab ts=4 sw=4:
# kate: space-indent on; indent-mode cstyle; indent-width 4; mixedindent off;
