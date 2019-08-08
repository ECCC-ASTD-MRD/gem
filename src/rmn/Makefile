SHELL = /bin/bash

include rmnlib_version.inc

LIBNAME = librmn_$(RMNLIB_VERSION)

WORKDIR = ./WorkDir

LIBDIR = `pwd`

DEBUG = No

default: genlib

genlib: 
	mkdir -p $(WORKDIR)
	mkdir -p $(LIBDIR)/$(EC_ARCH)
	if [ $(DEBUG) == "yes" ] ; \
	then \
	echo 'Compiling with DEBUG option' ; \
	sleep 2 ; \
	./make_locallib_packages-d ; \
	rm -f $(LIBDIR)/$(EC_ARCH)/$(LIBNAME)_d.a ; \
	./merge_rmnlib_packages $(WORKDIR) $(LIBDIR) $(LIBNAME)_d ; \
	else \
	./make_locallib_packages ; \
	rm -f $(LIBDIR)/$(EC_ARCH)/$(LIBNAME).a ; \
	./merge_rmnlib_packages $(WORKDIR) $(LIBDIR) $(LIBNAME) ; \
	fi

clean:
	cd template_utils/gmm ; make veryclean
	find . -name '*.o' -exec rm {} \;

distclean: clean
	find . -name lib_local.a -exec rm {} \;
