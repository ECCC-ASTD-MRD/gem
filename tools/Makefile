################################### LICENCE BEGIN ###############################
# GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
# Copyright (C) 2007-2013 - Air Quality Research Division &
#                           National Prediction Operations division
#                           Environnement Canada
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
################################### LICENCE END #################################

SHELL = /bin/sh

INCLUDE_DIR=../lib/$(ARCH)
FC                               = s.compile 

FFLAGS_AIX-powerpc7              = -O 2 -openmp -includes $(INCLUDE_DIR)  -optf="-qmaxmem=-1"
FFLAGS_Linux_x86-64              = -O 2 -openmp -includes $(INCLUDE_DIR) -optf="-Minfo=all"

FFLAGS_DEBUG_AIX-powerpc7        = -debug -defines=\"-DDEBUG_OUTPUT\" -includes $(INCLUDE_DIR) -openmp -optf="-qcheck -qflttrap=zero:imp:ov:und:inv:en -qhalt=w -qlinedebug "
FFLAGS_DEBUG_Linux_x86-64        = -O 0 -debug  -includes $(INCLUDE_DIR) -openmp -optf="-C -Kieee -Ktrap=align -Ktrap=denorm -Ktrap=divz -Ktrap=fp -Ktrap=inv -Ktrap=ovf -Ktrap=unf"

NBJOBS_Linux_x86-64              = 8
NBJOBS_AIX-powerpc7              = 16

LFLAG_Linux_x86-64               = -libpath ../lib/$(ARCH) -libappl "rpn_comm_40007 gemmach" -librmn rmn_013
LFLAG_AIX-powerpc7               = -libpath ../lib/$(ARCH) -libappl "rpn_comm_40007 gemmach" -librmn rmn_013
MOD      = print_species_master_$(BASE_ARCH).abs

.PHONY: clean distclean

.SUFFIXES: $(SUFFIXES) .ftn90

.FFILES =	print_species_master.ftn90

.OFILES  = $(.FFILES:.ftn90=.o)

.OMODULES = $(.FMODULES:.ftn90=.o)

.ftn90.o:
	@echo "Compilation de $(<)"
	$(FC) $(FFLAGS_$(DEBUG)$(BASE_ARCH)) $(TENDENCIES) -src $<

all:${.OMODULES} $(.OFILES)
	@echo "Target all"
	r.build -obj ${.OMODULES} $(.OFILES) -o $(MOD) $(LFLAG_$(BASE_ARCH))
debug:
	@echo "Target debug"
	$(MAKE) DEBUG=DEBUG_ all

tendencies:
	$(MAKE) $(MAKEFLAGS) TENDENCIES='-defines="-DMACH_TENDENCIES"' all

distclean:
	rm -fr *.o *.f90 .fo gmon.out  $(MOD) $(MOD)*.c $(MOD)*.o *.mod *~

clean:
	rm -fr $(.OFILES) $(.FFILES:.ftn90=.f90) .fo gmon.out $(MOD)*.c $(MOD)*.o *.mod
