#!/usr/bin/env python3

#/* Part of the Maestro sequencer software package.
# * Copyright (C) 2011-2015  Canadian Meteorological Centre
# *                          Environment Canada
# *
# * Maestro is free software; you can redistribute it and/or
# * modify it under the terms of the GNU Lesser General Public
# * License as published by the Free Software Foundation,
# * version 2.1 of the License.
# *
# * Maestro is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# * Lesser General Public License for more details.
# *
# * You should have received a copy of the GNU Lesser General Public
# * License along with this library; if not, write to the
# * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
# * Boston, MA 02111-1307, USA.
# */

#-------------------------------------------------------------------
# task_setup.py
#
# Module / executable to perform task setup operations
#-------------------------------------------------------------------

"""Create and fill a task runtime directory

INTRODUCTION
This task setup utility makes use of a pair of internal classes
('Section' and 'Config') to convert information contained in a
configuration file into the complete layout of a task directory.
The primary user access to this utility is through the 'Config'
class.

CONFIG CLASS

  METHODS
    cfg = Config('/path/to/config/file.cfg') - Class constructor.  This
         method reads and caches the contents of the named input
         configuration file.  The output from this method is a instance
         of the Config class.
    cfg.getSections() - Parse the sections of the configuration file
         to generate a list of 'links', 'targets' and 'options'.  These
         values are stored internally by the instance.
    cfg.link() - Generate the subdirectories and links to the files
         identified in the config file.
    cfg.setOption(option,value) - Set the named option ('delimiter_exec',
         'verbosity','cleanup','force') to the value specified in the argument.
         This method should be called before the 'getSections' method
         to ensure that keywords are properly resolved.

  CLASS VARIABLES
    configData - Cached copy of the data read from the configuration file.
    configFile - Name of the configuration file.
    taskdir    - User-specified (not real-path'd) path to task directory.
    basepath   - Real path to working directory below task level.
    taskname   - Task name.
    subdir_sectionMap - Name mapping from configuration file sections to
         task subdirectories.
    verbosity  - Integer to control verbosity level.
    cleanup    - Boolean to clean task directory before setup.
    force      - Boolean to force actions despite warnings.
    error      - Error code for return.
    ok         - Successful completion code for return.

SECTION CLASS

  METHODS
    s = Section(section) - Class constructor.  This method returns an
         instance of the 'Section' class for the particular section
         identified in the argument list.
    s.add(line) - Add the contents of a configuration file line to
         this this instance of the Section class.  This information
         will be appended to any existing additions made to the instance.

  CLASS VARIABLES
    delimiter_exec   - Delimiter for embedded commands (default '`')
    delimiter_target - Delimiter for multiple targets on the RHS (default \s or \n)
    verbosity        - Integer to control verbosity level.
    cleanup          - Boolean to clean task directory before setup.
    force            - Force action despite warnings.
"""

__version__ = "0.16.0"
__author__  = "Ron McTaggart-Cowan (ron.mctaggart-cowan@ec.gc.ca)"

#---------
# Imports
#---------
import os
import sys
import shutil
import re
import optparse
import tempfile
import types
import shlex
import copy
from time import time

class Store(object):
    """Space for saving values using a callable object"""
    def __init__(self):
        """Class constructor"""
        self.saved = []
    def __call__(self,value):
        """Add an entry to the saved space"""
        self.saved.append(value)

def mkdir_p(path):
    import os,sys,errno
    try:
        os.makedirs(path)
    except OSError:
        value = sys.exc_info()[1][0]
        if value == errno.EEXIST:
            pass
        else:
            sys.stderr.write('task_setup.py::os.makedirs() returned the following error information on an attempt to create ' \
                             +path+': '+str(sys.exc_info())+"\n")
            raise

def which(name,path=None,verbose=True):
    """Duplicates the functionality of UNIX 'which' command"""
    if re.search('/',name):
        return(name)
    bin_path = path and path or os.environ['PATH']
    for dir in re.split(':',bin_path):
        fullname=os.path.join(dir,name)
        try:
            if os.path.isfile(fullname):
                if os.access(fullname,os.X_OK): return(fullname)
        except:
            continue
    if (verbose): print("Warning: unable to find "+name+" in path:\n"+bin_path)
    return('')

def path2host(machine,path):
    """Convert a machine/abspath pair to the heirarchical part of a URI"""
    if machine:
        return(machine+':'+path)
    else:
        return(path)

def resolveKeywords(entry,delim_exec='',set=None,verbose=False,internals={}):
    """Resolve special keywords in the entry (no procesing for keywords in embedded commands)"""
    delim_start='\$\{'
    delim_end='}'
    delim = re.compile(delim_start+'(.*?)'+delim_end)
    dollar = re.compile('\$')
    elements = delim_exec and re.split(delim_exec+'(.*?)'+delim_exec,entry) or [entry]
    found_internal = False
    for i in range(0,len(elements)):
        element_orig=elements[i]
        if i%2:
            # This is an embedded command.  Add delimiters and do not substitute keywords.
            elements[i] = delim_exec+elements[i]+delim_exec
        else:
            # This is a standard string.  Attempt to replace all keywords.
            for keyword in delim.findall(elements[i]):
                vartype = 'unknown'
                if not keyword: continue
                if keyword in internals:
                    this_keyword = internals[keyword]
                    found_internal = True
                    vartype = 'internal'
                if vartype != 'internal':
                    try:
                        this_keyword = os.environ[keyword]
                        vartype = 'found'
                    except KeyError:
                        vartype = 'environment'
                    if set:
                        try:
                            this_keyword = set[keyword]
                            vartype='found'
                        except KeyError:
                            vartype = vartype+'/set'
                    if vartype != 'found':
                        if not keyword in undef_list.saved:
                            warnline = "Warning: "+vartype+" variable "+keyword+" undefined ... empty substitution performed"
                            sys.stderr.write(warnline+'\n')
                            if (verbose): print(warnline)
                            undef_list(keyword)
                        this_keyword = ''
                elements[i] = re.sub(delim_start+keyword+delim_end,this_keyword,elements[i])
            # Final substitution attempt to support deep internal indexing
            for keyword in delim.findall(elements[i]):
                this_keyword = ''
                if not keyword: continue
                if keyword in internals:
                    this_keyword = internals[keyword]
                    found_internal = True
                elements[i] = re.sub(delim_start+keyword+delim_end,this_keyword,elements[i])
            # Check for leftover $ symbols and generate error message
            if dollar.search(elements[i]):
                warnline="Error: found a $ character after resolution of "+element_orig+" to "+elements[i]+ \
                          "\n  The result of external keyword resolution cannot contain un-expanded shell variables.  Evaluate the\n"+\
                          "  string or remove extra quoting / escape characters before the task_setup call to avoid this problem.  "
                sys.stderr.write(warnline+'\n')
                if (verbose): print(warnline)
    updated = ''.join(elements)
    return({'string':updated,'contains_internal':found_internal})

def getRealPath(node,verbosity):
    """Get the real path of a file/directory"""
    if node == "": return ""
    have_subprocess=True
    if (int(verbosity) >= 2): startTime=time()
    try:
        import subprocess
    except ImportError:
        have_subprocess=False
    try:
        get_real_path = "realpath "+node
        if have_subprocess:
            p = subprocess.Popen(get_real_path,shell=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
            true_src = p.stdout.read()
            error_out = p.stderr.read()
        else:
            (stdin,stdout,stderr) = os.popen3(get_real_path,'r')
            true_src = stdout.read()
            error_out = stderr.read()
            stdin.close()
            stdout.close()
            stderr.close()
        if true_src == '(null)' or not true_src or re.search('No such file or directory$',true_src,re.M):
            print("Warning: real_path on " + node + " returned " + error_out)
            true_src = node
    except OSError:
        if (os.path.exists(node)):
            print("Warning: real_path does not exist or returned an error for "+src_file)
        true_src = node
    if (int(verbosity) >= 2): print("Info 2: getRealPath exec time: " + str( time() - startTime))
    return(true_src)

class LinkFile():
    """Structure for link file target information"""

    def __init__(self,link,target_host,target,link_only,verbosity=False):
        """Class constructor"""
        self.have_subprocess=True
        try:
            import subprocess
        except ImportError:
            self.have_subprocess=False
        self.link = link
        self.target_host = target_host
        self.target = target
        self.link_only = link_only
        self.verbosity = verbosity
        self.src = []
        self.host = []
        self._expandTarget()
        self._trueSources()
        self._setPrefixes()
        self._sourceTypes()

    def _expandTarget(self):
        """Complete target information through local or remote wildcard expansion"""
        import glob
        for i in range(0,len(self.target)):
            src_expanded = []
            try:
                hostname = self.target_host[i]
            except TypeError:
                hostname = None
            src_expanded = glob.glob(self.target[i])
            if len(src_expanded) < 1 and hostname:
                file_sep = '?'
                file_list = "ssh "+hostname+" \"python -c 'import glob; f=glob.glob(\\\""+self.target[i]+"\\\"); print(\\\""+file_sep+"\\\".join(f))'\""
                if self.have_subprocess:
                    import subprocess
                    p = subprocess.Popen(file_list,shell=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                    output = p.stdout.read()
                    error = p.stderr.read()
                else:
                    (stdin,stdout,stderr) = os.popen3(file_list,'r')
                    output = stdout.read()
                    error = stderr.read()
                src_expanded = [fname.rstrip('\n') for fname in re.split('\\'+file_sep,output.rstrip('\n')) if fname]
            if len(src_expanded) < 1:
                src_expanded = [self.target[i]]
            self.src.extend(src_expanded)
            self.host.extend([hostname for item in src_expanded])

    def _trueSources(self):
        """Get true source paths for entries"""
        self.true_src_file = [getRealPath(src_file,self.verbosity) for src_file in self.src]

    def _setPrefixes(self):
        """Set hosts and prefixes for entries"""
        self.src_file_prefix = [target_host and target_host+':' or '' for target_host in self.host]

    def _sourceTypes(self):
        """Determine the type of source"""
        self.remote_file_type = ['' for src in self.true_src_file]
        for host in set(self.host):
            if not host: continue
            idx = []
            for i in range(0,len(self.true_src_file)):
                if self.host[i] == host:
                    idx.append(i)
            check_file = "ssh "+host+" '"
            for i in idx:
                check_file += ' if [[ -d "'+self.true_src_file[i]+ \
                    '" ]] ; then echo 2 ; elif [[ -f "'+self.true_src_file[i]+ \
                    '" ]] ; then echo 1 ; else echo 0 ; fi;'
            check_file.rstrip(';')
            check_file += "'"
            if self.have_subprocess:
                import subprocess
                p = subprocess.Popen(check_file,shell=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                output = p.stdout.read().rstrip('\n')
                error = p.stderr.read()
            else:
                (stdin,stdout,stderr) = os.popen3(check_file,'r')
                output = stdout.read().rstrip('\n')
                error = stderr.read()
            if len(error) > 0:
                warnline = "Warning: STDERR returned from "+self.host[i]+" is "+error
                sys.stderr.write(warnline+'\n')
                if (self.verbosity): print(warnline)
            if len(output) > 0:
                output_list = output.split('\n')
                for i in range(0,len(idx)):
                    try:
                        ftype = int(output_list[i])
                    except:
                        if not self.link_only:
                            print("Warning: required file "+self.true_src_file[i]+" does not exist on host "+self.host[i])
                    if ftype == 1:
                        self.remote_file_type[idx[i]] = 'file'
                    elif ftype == 2:
                        self.remote_file_type[idx[i]] = 'directory'
            else:
                print("Warning: unable to login to target host "+self.host[i]+". See previous error statement for STDERR details.")

    def rephost(self):
        """Repeat host entry for all targets"""
        self.src.extend(self.target)
        self.host.extend([self.target_host[0] for item in self.target])

class Section(list):
    """Data and functions applicable to individual configuration sections"""

    # Class variables
    delimiter_exec = '`'
    delimiter_target = '(?<!<no)\n|\s+(?!value>)'
    verbosity = 0
    cleanup = False
    force = False

    def __init__(self,section,set=None,cfg=None,attrib={},varcache=None):
        """Class constructor"""
        self.section = section
        self.set = set
        self.cfg = cfg
        self.attrib = attrib
        self.varcacheFile = varcache
        no_loop = {'var':None,'steps':list('0')}
        self.loop = copy.copy(no_loop)
        if self._isType('loop'):
            try:
                self.loop['var'] = self.attrib['var']
                step = 'step' in self.attrib and int(self.attrib['step']) or 1
                format_string = "%0"+str(len(self.attrib['start']))+"d"
                self.loop['steps'] = [format_string % (i) for i in range(int(self.attrib['start']),int(self.attrib['end'])+1,step)]
            except KeyError:
                self.loop = no_loop
                warnline = "Warning: incomplete loop specification for <"+section+"> - no looping for this section"
                sys.stderr.write(warnline+'\n')
                if (self.verbosity): print(warnline)
            except ValueError:
                self.loop = no_loop
                warnline = "Warning: invalid loop specification for <"+section+"> - no looping for this section"
                sys.stderr.write(warnline+'\n')
                if (self.verbosity): print(warnline)

    def _isType(self,check_type):
        """Determine whether this section is of a specific type"""
        return('type' in self.attrib and self.attrib['type'] or None)

    def _splitHost(self,entry):
        """Split a set of strings into host:path form"""
        hostpath = {'host':[],'path':[]}
        for item in entry:
            try:
                (host,path) = re.split(':',item)
                host_noquote = re.sub('[\'\"]','',host)
            except ValueError:
                path = item
                host_noquote = None
            hostpath['host'].append(host_noquote)
            hostpath['path'].append(path)
        return(hostpath)

    def _sectionResolveKeywords(self,entry,internals=None):
        """Resolve special keywords in the entry"""
        return resolveKeywords(entry,delim_exec=self.delimiter_exec,set=self.set,verbose=self.verbosity,internals=internals)

    def _executeEmbedded(self,entry,internals={}):
        """Execute backtic embedded commands and substitute result"""
        have_subprocess=True
        try:
            import subprocess
        except ImportError:
            have_subprocess=False
        updated = [entry]
        delim = re.compile(self.delimiter_exec+'(.*?)'+self.delimiter_exec)
        shell_dot_config = (self.cfg) and '. '+self.cfg+' >/dev/null 2>&1; ' or 'true; '
        if (self.varcacheFile):
            shell_gen_cachefile = 'task_setup_cachegen '+self.cfg+' '+self.varcacheFile+' ; . '+self.varcacheFile+' ; '
            command_prefix = 'if [[ -s '+self.varcacheFile+' ]] ; then . '+self.varcacheFile+' >/dev/null 2>&1 ; else '+shell_gen_cachefile+'fi ; '
        else:
            command_prefix = shell_dot_config
        for command in delim.finditer(entry):
            for var in list(internals.keys()):
                command_prefix = command_prefix+str(var)+'='+str(internals[var])+'; '
            if have_subprocess:
                p = subprocess.Popen(command_prefix+command.group(1),shell=True,universal_newlines=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
                error_message = p.stderr.read().rstrip('\n')
                outbuf = p.stdout.read().rstrip('\n ')
            else:
                (stdin,stdout,stderr) = os.popen3(command_prefix+command.group(1),'r')
                error_message = stderr.read().rstrip('\n')
                outbuf = stdout.read().rstrip('\n ')
                stdin.close()
                stdout.close()
                stderr.close()
            elements = re.split(self.delimiter_target,outbuf)
            target_list = []
            for j in range(0,len(updated)):
                for i in range(0,len(elements)):
                    target_list.append(command.re.sub(elements[i],updated[j],count=1))
            updated = target_list
            if error_message:
                print("Warning: the embedded command "+self.delimiter_exec+command.string+self.delimiter_exec+ \
                      " in the configuration file returned an error: "+error_message)
        return(updated)

    def add(self,line,search_path):
        """Add data to the section"""
        data = re.split('\s+',re.sub('^(#)+',' ',line))
        entry = {}
        try:
            rawLink = data[1]
            rawTarget = ' '.join(data[2:]).rstrip()
        except IndexError:
            warnline = "Warning: ignoring malformed configuration line: "+line
            sys.stderr.write(warnline)
            if (self.verbosity): print(warnline)
            return(False)
        bin_path = None
        if search_path:
            try:
                bin_path = self.set['PATH']
            except:
                pass
        lastSlash = re.compile('/$',re.M)
        noval = re.compile('^\s*[\'\"]*<no\svalue>',re.M)
        comment = re.compile('^#',re.M)
        for step in self.loop['steps']:
            loopInternals={self.loop['var']:step}
            link = self._sectionResolveKeywords(lastSlash.sub('',rawLink),internals=loopInternals)
            link_split = self._splitHost([link['string']])
            entry["link_host"] = link_split["host"][0]
            entry["link"] = link_split["path"][0]
            target = self._sectionResolveKeywords(rawTarget,internals=loopInternals)
            target_executed = [str(item).replace("'","").rstrip() for item in self._executeEmbedded(target['string'],internals=loopInternals)]
            target_list = re.split(self.delimiter_target,' '.join(target_executed))
            target_split = self._splitHost(target_list)
            if ([True for target_string in target_split["path"] if noval.match(target_string)]):
                if any([True for target_string in target_split["path"] if not noval.match(target_string)]):
                    print("Info 1: some entries for "+link['string']+" will not be added because of special target value '<no value>'")
                else:
                    print("Info 1: will not create link for "+link['string']+" because of special target value '<no value>'")
                    continue
            if step != self.loop['steps'][0] and not (link['contains_internal'] or target['contains_internal']):
                continue
            entry["target"] = []
            entry["target_host"] = []
            for i in range(0,len(target_split["path"])):
                if comment.match(target_split["path"][i]):
                    break
                if not noval.match(target_split["path"][i]):
                    entry["target_host"].append(target_split["host"][i])
                    entry["target"].append(target_split["path"][i])
            entry["target_type"] = (lastSlash.search(rawLink) or len(entry["target"]) > 1) and 'directory' or 'file'
            if search_path:
                entry["target"] = [which(target,path=bin_path) for target in entry["target"]]
            entry["copy"] = False
            entry["cleanup"] = False
            entry["create_target"] = False
            entry["link_only"] = False
            if self.section == 'output':
                entry["create_target"] = True
                entry["link_only"] = True
            self.append(copy.deepcopy(entry))

class Config(dict):
    """Data and functions applicable to the task setup"""

    # Class variables
    configData = {}
    configFile = None
    taskdir = None
    basepath = None
    taskname = None
    verbosity = 0
    cleanup = False
    force = False
    error = 0
    ok = 1
    subdir_sectionMap = {'input':       'input',
                         'executables': 'bin',
                         'work':        'work',
                         'output':      'output',
                         'setup':       '.setup'}  #Tags in config files (keys) are mapped to these subdir names (values)
    force_sections = ['work','setup']              #Force the creation of these sections regardless of config file contents
    search_path_sections = ['executables','setup'] #These sections will search the PATH for non-fully-qualified targets
    ignore_sections = ['seq_scheduler']            #Ignore these sections in the configuration file
    varcache_name = 'task_setup_varcache.txt'      #Name for environment caching in embedded commands

    def __init__(self,file=None,taskdir=None,set=None,varcache=None):
        """Class constructor"""
        self.configFile = file
        self.taskdir = taskdir
        self.setFile = set
        self.set = None
        self.sectionList = []
        self.callFile = self._createTmpFile(sys.argv)
        self.envFile = self._createTmpFile(os.environ)
        self.varcacheFile = (varcache) and varcache or self._createTmpFile(None)
        self["file"] = file
        if set:
            self._readSetFile(set)
        if not self.configData:
            self._readConfigFile(self["file"])

    def __del__(self):
        """Class destructor"""
        os.unlink(self.callFile)
        os.unlink(self.envFile)

    def _createTmpFile(self,contents):
        """Create and fill a temporary file, returning the file name"""
        try:
            (fdunit,filename) = tempfile.mkstemp()
            if not contents: return(filename)
            fd = os.fdopen(fdunit,"w")
        except OSError:
            print("Warning: Unable to create temporary file for call statement")
            return(None)
        if isinstance(contents, dict):
            keys = list(contents.keys())
            keys.sort()
            for key in keys:
                fd.write(str(key)+'='+str(contents[key])+'\n')
        elif isinstance(contents, list):
            fd.write(' '.join(contents)+'\n')
        else:
            fd.write(str(contents)+'\n')
        fd.close()
        return(filename)

    def _readSetFile(self,file):
        """Read set file"""
        try:
            fd = open(file,"r")
            setData = fd.readlines()
        except IOError:
            print("Warning: unable to read set from "+file)
            self.set = None
            return()
        fd.close()
        sep = re.compile(r"(?<!\\)'")
        equal = re.compile("=")
        quote_count = 0
        self.set = {}; concat_line = ''
        for line in setData:
            quote_count += len(re.findall(sep,line))
            if quote_count%2 == 0:
                concat_line += line
                try:
                    (key,value) = equal.split(concat_line,maxsplit=1)
                except ValueError:
                    if quote_count == 0:
                        concat_line = ''
                    continue
                self.set[key] = value.rstrip('\n')
                concat_line = ''
            else:
                concat_line += line

    def _readConfigFile(self,file):
        """Read configuration file"""
        if not file:
            self.configData = ''
            return
        try:
            fd = open(file,"r")
            try:
                self.configData = fd.readlines()
            finally:
                fd.close()
        except IOError:
            warnline = "Warning: unable to read from configuration file "+file
            sys.stderr.write(warnline+'\n')
            if (self.verbosity): print(warnline)
            self.configData = []
            self.configFile = '/dev/null'

    def _map(self,section):
        """Map a section name to a task subdirectory name"""
        try:
            subdir = self.subdir_sectionMap[section]
        except KeyError:
            print("Warning: unknown section "+section+" encountered ... no mapping done")
            return(section)
        return(subdir)

    def _subdir_setup(self,subdir):
        """Set up the requested subdirectory for the task"""
        status = self.ok
        if not os.path.isdir(subdir):
            if (self.verbosity): print("Info 1: creating subdirectory "+subdir)
            try:
                mkdir_p(subdir)
            except OSError:
                print("Error: could not create "+subdir)
                status = self.error
        return(status)

    def _get_subdirs(self,dir,absolute=True):
        """Return a list of relative or absolute expected subdirectories"""
        subdirs = [self._map(section) for section in list(self["sections"].keys())]
        subdirs.sort()
        if absolute:
            return([os.path.join(dir,subdir) for subdir in subdirs])
        else:
            return(subdirs)

    def _taskdir_setup(self):
        """Set up task base directory"""
        status = self.ok
        if self.cleanup:
            if os.path.isdir(self.taskdir):
                contents = [entry for entry in os.listdir(self.taskdir) if os.path.isdir(os.path.join(self.taskdir,entry))]
                if len(contents) > 0:
                    if self.force:
                        for sub in contents:
                            try:
                                shutil.rmtree(os.path.join(self.taskdir,sub))
                            except:
                                print("Error: unable to force clean workspace subdirectory "+sub)
                                return(self.error)
                    else:
                        contents.sort()
                        if contents == self._get_subdirs(self.taskdir,absolute=False):
                            for sub in self._get_subdirs(self.taskdir,absolute=True):
                                try:
                                    shutil.rmtree(sub)
                                except:
                                    print("Error: unable to remove task subdirectory "+sub)
                                    return(self.error)
                        else:
                            print("Error: Invalid and/or changed subdirectory <-> section mapping in task_setup.py.")
                            print("   The requested task base directory "+self.taskdir+" contains a subdirectory that")
                            print("   is not recognized based on the configuration file "+self["file"]+"  If")
                            print("   this is a valid task base directory, please remove it manually and relaunch.")
                            print("       Task subdirectories: "+str(contents))
                            print("       Mapped config sections: "+str(self._get_subdirs(self.taskdir,absolute=False)))
                            return(self.error)
        if not os.path.isdir(self.taskdir):
            try:
                mkdir_p(self.taskdir)
            except OSError:
                print("Error: could not create task base directory "+self.taskdir)
                return(self.error)
        elif not os.access(self.taskdir,os.W_OK):
            print("Error: task directory "+self.taskdir+" is not writeable ... exiting")
            return(self.error)
        # Set task name and working path (needs to exist for `real_path` so it can't be done during construction)
        basedir = getRealPath(self.taskdir,self.verbosity)
        self.basepath = os.path.dirname(basedir)
        self.taskname = os.path.basename(basedir)
        return(status)

    def _append_meta(self,section,meta):
        """Append metadata to the specified section"""
        status = self.ok
        try:
            self["sections"][section].append(meta)
        except KeyError:
            status = self.error
        return(status)

    def _special_appends(self):
        """Add special values to sections"""
        self._append_meta("setup",{"link":"task_setup",
                                   "target":[sys.argv[0]],
                                   "target_type":'file',
                                   "target_host":[None],
                                   "copy":False,
                                   "cleanup":False,
                                   "create_target":False,
                                   "link_host":None,
                                   "link_only":False})
        if self["file"]:
            self._append_meta("setup",{"link":"task_setup.cfg",
                                       "target":[self.configFile],
                                       "target_type":'file',
                                       "target_host":[None],
                                       "copy":True,
                                       "cleanup":False,
                                       "create_target":False,
                                       "link_host":None,
                                       "link_only":False})
        if self.varcacheFile:
            self._append_meta("setup",{"link":"task_setup_varcache.txt",
                                       "target":[self.varcacheFile],
                                       "target_type":'file',
                                       "target_host":[None],
                                       "copy":True,
                                       "cleanup":True,
                                       "create_target":False,
                                       "link_host":None,
                                       "link_only":False})
        self._append_meta("setup",{"link":"task_setup_call.txt",
                                   "target":[self.callFile],
                                   "target_type":'file',
                                   "target_host":[None],
                                   "copy":True,
                                   "cleanup":False,
                                   "create_target":False,
                                   "link_host":None,
                                   "link_only":False})
        self._append_meta("setup",{"link":"task_setup_env.txt",
                                   "target":[self.envFile],
                                   "target_host":[None],
                                   "target_type":'file',
                                   "copy":True,
                                   "cleanup":False,
                                   "create_target":False,
                                   "link_host":None,
                                   "link_only":False})
        if self.setFile:
            self._append_meta("setup",{"link":"task_setup_set.txt",
                                       "target":[self.setFile],
                                       "target_type":'file',
                                       "target_host":[None],
                                       "copy":True,
                                       "cleanup":True,
                                       "create_target":False,
                                       "link_host":None,
                                       "link_only":False})
        cachegen = which('task_setup_cachegen',verbose=self.verbosity)
        if cachegen:
            self._append_meta("setup",{"link":"task_setup_cachegen",
                                       "target":[cachegen],
                                       "target_type":'file',
                                       "target_host":[None],
                                       "copy":False,
                                       "cleanup":False,
                                       "create_target":False,
                                       "link_host":None,
                                       "link_only":False})
        real_path=which('realpath',verbose=self.verbosity)
        if real_path:
            self._append_meta("setup",{"link":"task_setup_realpath",
                                       "target":[real_path],
                                       "target_type":'file',
                                       "target_host":[None],
                                       "copy":False,
                                       "cleanup":False,
                                       "create_target":False,
                                       "link_host":None,
                                       "link_only":False})
        return(self.ok)

    def _createTarget(self,entry,host,path):
        """Create target directory"""
        have_subprocess=True
        try:
            import subprocess
        except ImportError:
            have_subprocess=False
        status = self.ok
        if not entry["create_target"]: return(status)
        directory = (entry["target_type"] == 'directory') and path or os.path.split(path)[0]
        if not directory:
            print("Error: no directory specified target in request for "+entry["link"])
            status = self.error
            return(status)
        if host:
            make_dir = "echo \"s.mkdir_onebyone "+directory+"; if [[ -d "+directory+ \
                       " ]] ; then echo TASK_SETUP_SUCCESS ; else echo TASK_SETUP_FAILURE ; fi\" | ssh "+ \
                       host+" bash --login"
            if have_subprocess:
                p = subprocess.Popen(make_dir,shell=True,universal_newlines=True,stderr=subprocess.PIPE,stdout=subprocess.PIPE)
                error = p.stderr.read()
                output = p.stdout.read()
            else:
                (stdin,stdout,stderr) = os.popen3(make_dir,'r')
                error = stderr.read()
                output = stdout.read()
            if not re.search("TASK_SETUP_SUCCESS",output):
                status = self.error
                if re.search("TASK_SETUP_FAILURE",output):
                    print("Error: login to "+host+" successful but "+directory+" not created")
                else:
                    print("Error: unable to obtain directory status on "+host)
            if len(error) > 0:
                sys.stderr.write("task_setup.py::_createTarget() attempt to connect to "+host+" returned STDERR "+error+"\n")
        else:
            if not os.path.isdir(directory):
                try:
                    mkdir_p(directory)
                    if (self.verbosity): print("Info 1: created directory "+directory+" to complete target request")
                except:
                    print("Error: unable to create "+directory+" to complete target request")
                    status = self.error
        return(status)

    def _parseSectionHead(self,head):
        """Parse section header into individual attributes"""
        head = resolveKeywords(head,set=self.set,verbose=self.verbosity)
        try:
            att_string = re.split('\s+',head['string'],maxsplit=1)[1]
        except IndexError:
            return({})
        return(dict(token.split('=') for token in shlex.split(att_string)))

    def setOption(self,option,value):
        """Option handling dispatcher"""
        try:
            getattr(Section,option)
        except AttributeError:
            print("Error: attempt to change invalid setting "+option)
            return (self.error)
        setattr(Section,option,value)
        setattr(self,option,value)
        return(self.ok)

    def getSections(self):
        """Break input data into individual sections"""
        currentSection = None
        prefix='^\s*#\s*'
        validLine = re.compile(prefix+'[^#](.+)',re.M)
        sectionHead = re.compile(prefix+'<([^/]\S+)(.*)>',re.M)
        sectionFoot = re.compile(prefix+'</(.*)>',re.M)
        self["sections"] = {}
        for raw_line in self.configData:
            line = re.sub('^\s+','',raw_line,re.M)
            head = False
            valid = validLine.search(line)
            if (valid):
                foot = sectionFoot.search(line)
                if foot and currentSection:
                    if foot.group(1) != currentSection:
                        print("Warning: section head <"+currentSection+"> does not match the section foot </"+foot.group(1)+"> in "+self["file"])
                    currentSection = None
                else:
                    head = sectionHead.search(line)
                    if head:
                        if currentSection:
                            print("Error: found header for "+head.group(1)+" while still in open section for "+currentSection)
                            print("  Perhaps the configuration file "+self["file"]+" is missing an </"+currentSection+"> end section tag?")
                            sys.stderr.write("task_setup.py::getSections() failed parsing "+self["file"]+" at <"+head.group(1)+">\n")
                            self["sections"] = {}
                            return(self.error)
                        currentSection = head.group(1)
                        if currentSection in self.ignore_sections:
                            currentSection = None
                        else:
                            headAttrib = self._parseSectionHead(head.group(2))
                            self["sections"][currentSection] = Section(currentSection,set=self.set,cfg=self["file"],attrib=headAttrib,varcache=self.varcacheFile)
                            self.sectionList.append(currentSection)
                if (currentSection and not head):
                    self["sections"][currentSection].add(line,currentSection in self.search_path_sections)
        for force in self.force_sections:
            self["sections"][force] = Section(force)
        self._special_appends()
        return(self.ok)

    def write(self,fd):
        """Write the config file sections"""
        for section in self.sectionList:
            fd.write('#<'+section+'>\n')
            for entry in self["sections"][section]:
                append = (entry["target_type"] == 'directory') and '/' or ''
                target = ''
                for i in range(0,len(entry["target"])):
                    host = (entry["target_host"][i]) and entry["target_host"][i]+':' or ''
                    target += ' '+host+entry["target"][i]
                fd.write('# '+entry["link"]+append+' '+target+'\n')
            fd.write('#</'+section+'>\n')

    def link(self):
        """Perform subdirectory creation and linking operations"""
        have_subprocess=True
        try:
            import subprocess
        except ImportError:
            have_subprocess=False
        status = self.ok
        sub_status = self._taskdir_setup()
        if sub_status != self.ok: return(sub_status)
        for section in list(self["sections"].keys()):
            if (self.verbosity): print("  <"+section+">")
            abs_subdir = os.path.join(self.taskdir,self._map(section))
            sub_status = self._subdir_setup(abs_subdir)
            if sub_status != self.ok: return(sub_status)
            for entry in self["sections"][section]:
                if (int(self.verbosity) >= 2): startTime=time()
                line = LinkFile(entry["link"],entry["target_host"],entry["target"],entry["link_only"],verbosity=self.verbosity)
                if len(line.target) == 0:
                    print("Error: empty target for "+line.link+" ... skipping")
                    status = self.error
                    continue
                link_only = entry["link_only"]
                dest = os.path.join(abs_subdir,entry["link"])
                if not os.path.isdir(os.path.dirname(dest)):
                    mkdir_p(os.path.dirname(dest))
                if os.path.islink(dest): os.remove(dest)
                dest_is_dir = False
                if len(line.src) == 0:
                    line.rephost()
                elif entry["target_type"] == 'directory' and not link_only or len(line.src) > 1:
                    dest_is_dir = True
                    if not os.path.isdir(dest):
                        try:
                            mkdir_p(dest)
                        except OSError:
                            print("Error: could not create "+section+" subdirectory "+dest)
                            dest_is_dir = False
                            status = self.error

                # Process each file on the line separately
                for i in range(len(line.src)-1,-1,-1):

                    # Retrieve information about the source file
                    true_src_file = line.true_src_file[i].rstrip('\n')
                    src_file_prefix = line.src_file_prefix[i]

                    # Retrieve information about the destination
                    if dest_is_dir and not link_only:
                        dest_file = os.path.join(dest,os.path.basename(line.src[i]))
                    else:
                        dest_file = dest
                    dest_path_short = dest_file.replace(self.taskdir,'')

                    # Check that the source file information is valid
                    if not true_src_file:
                        print("Error: skipping entry because no source file given for "+dest_path_short)
                        status = self.error
                        continue

                    # Take care of creating directory links
                    if os.path.isdir(true_src_file) or line.remote_file_type[i] == 'directory':
                        if entry["target_type"] != 'directory':
                            if (self.verbosity): print("Warning: "+entry["target_type"]+" link "+entry["link"]+ \
                               " refers to a directory target "+str(entry["target"]))
                        if os.path.islink(dest_file):
                            print("Warning: updating directory link to "+dest_path_short+" => "+src_file_prefix+true_src_file+" (previous target was "+os.readlink(dest_file)+")")
                            os.remove(dest_file)
                        try:
                            os.symlink(path2host(line.host[i],true_src_file),dest_file)
                            if (self.verbosity): print("Info 1: linked directory "+dest_path_short+" => "+src_file_prefix+true_src_file)
                        except IOError:
                            print("Error: error creating symlink for directory "+dest_path_short+" => "+src_file_prefix+true_src_file)
                            status = self.error
                        except OSError:
                            status = self.error
                            if os.path.isdir(dest_file) and link_only:
                                print("Error: multiple entries for "+dest_path_short+" in the "+section+" section are not supported")
                            else:
                                raise

                    # Take care of creating file links or copies
                    else:
                        isfile = True
                        if line.remote_file_type[i] != 'file':
                            try:
                                fd = open(true_src_file,'r')
                            except IOError:
                                isfile = False

                        if isfile and entry["target_type"] != 'file' and len(line.src) == 1:
                            if (self.verbosity): print("Warning: "+entry["target_type"]+" link "+entry["link"]+ \
                               "/ refers to a file target "+str(entry["target"]))
                        if isfile or link_only:
                            try:
                                if entry["copy"] and not link_only:
                                    if entry["cleanup"]:
                                        shutil.move(true_src_file,dest_file)
                                        link_type = "moved"
                                    else:
                                        shutil.copyfile(true_src_file,dest_file)
                                        link_type = "copied"
                                else:
                                    if entry["create_target"]:
                                        status_create = self._createTarget(entry,line.host[i],true_src_file)
                                        if status == self.ok: status = status_create
                                        true_src_file = getRealPath(true_src_file,self.verbosity)
                                        if true_src_file == "":
                                           print("Error: attempting to create link to empty target string.")
                                           status = self.error
                                    if os.path.islink(dest_file):
                                        print("Warning: updating file link to "+dest_path_short+" => "+src_file_prefix+true_src_file+" (previous target was "+os.readlink(dest_file)+")")
                                        os.remove(dest_file)
                                    os.symlink(path2host(line.host[i],true_src_file),dest_file)
                                    link_type = "linked"
                                if (self.verbosity): print("Info 1: "+link_type+" file "+dest_path_short+" => "+src_file_prefix+true_src_file)
                            except OSError:
                                print("Error: error creating symlink for file "+dest_path_short+" => "+src_file_prefix+true_src_file)
                                raise
                                status = self.error
                        else:
                            print("Error: unable to link "+dest_path_short+" => "+src_file_prefix+true_src_file+" ... source file is unavailable")
                            status = self.error
                if (int(self.verbosity) >= 2): print(("Info 2: Link creation time: " + str( time() - startTime)))
            if (self.verbosity): print("  </"+section+">")
        return(status)

# Executable segment
if __name__ == "__main__":

    # Command line argument parsing
    usage = "For complete and up to date information on this command, see the man page by typing 'man task_setup.py'."
    parser = optparse.OptionParser(usage=usage)
    parser.add_option("-b","--base",dest="basedir",default='.',
                      help="task base DIRECTORY",metavar="DIRECTORY")
    parser.add_option("-v","--verbose",dest="verbose",action="count",
                      help="verbose runtime output",default=0)
    parser.add_option("-c","--clean",dest="clean",action="store_true",
                      help="clean task directory before setup",default=False)
    parser.add_option("-r","--force",dest="force",action="store_true",
                      help="force action (ignore warnings)",default=False)
    parser.add_option("-e","--environment",dest="environment",default=None,
                      help="text FILE containing the set namespace in which to run",metavar="FILE")
    parser.add_option("","--varcache",dest="varcache",default=None,
                      help="text FILE containing a 'sourceable' version of the set namespace",metavar="FILE")
    parser.add_option("-d","--dry-run",dest="dryrun",action="store_true",
                      help="handle configuration file without acting on it",default=False)
    (options,args) = parser.parse_args()

    # Ensure that the user has provided a configuration file
    try:
        cfgFile = args[0]
    except IndexError:
        cfgFile = None

    # Read, parse and act on configuration file for task setup
    undef_list = Store()
    cfg = Config(file=cfgFile,taskdir=options.basedir,set=options.environment,varcache=options.varcache)
    cfg.setOption('cleanup',options.clean)
    cfg.setOption('force',options.force)
    cfg.setOption('verbosity',options.verbose)
    if cfg.getSections():
        pass
    else:
        if cfg.verbosity: print(" *** Error: task_setup.py unable to continue *** ")
        sys.exit(1)
    if options.dryrun:
        cfg.write(sys.stdout)
        del cfg
        sys.exit(0)
    if cfg.link():
        del cfg
        sys.exit(0)
    else:
        if cfg.verbosity: print(" *** Error: problematic completion from task_setup.py *** ")
        del cfg
        sys.exit(1)
