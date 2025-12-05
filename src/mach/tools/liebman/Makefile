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
#. r.load.dot rpn/libs/16.2

SHELL = /bin/sh

FC                               = s.compile

RMNLIB=rmn_016.2

LFLAG_Linux_x86-64               = -O 2 -librmn ${RMNLIB}
#LFLAG_AIX-powerpc7               = -O 2 -librmn ${RMNLIB}
LFLAG_AIX-powerpc7               = -O 2 -librmn
MOD                              = liebman_main_$(BASE_ARCH).Abs



.PHONY: clean distclean

.SUFFIXES: $(SUFFIXES) .ftn90

.FFILES =   liebman_main.ftn90 \
            lieb_limite.ftn90 \
            negmask.ftn90 \

.OFILES  = $(.FFILES:.ftn90=.o)

.ftn90.o:
	@echo "Compilation de $(<)"
	$(FC) -src $<

all:$(.OFILES)
	@echo "Target all"
	s.compile -obj $(.OFILES) -o $(MOD) $(LFLAG_$(BASE_ARCH))

arch:
	@echo "Target arch"
	@echo $(BASE_ARCH) $(RMNLIB)  $(ARCH) 

distclean: clean
	@echo "Target distclean"
	rm -fr $(MOD) 

clean:
	@echo "Target clean"
	rm -fr *~ *.o *.f90 .fo gmon.out .fo gmon.out $(MOD)  rii_files
