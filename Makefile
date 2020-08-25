SHELL = /bin/bash

#Makefile for Environment Canada systems

GEM_COMPILER ?= "intel"

default: build

cmake-intel: ; $(MAKE) GEM_COMPILER=intel cmake
cmake-gnu: ; $(MAKE) GEM_COMPILER=gnu cmake
cmake:
	( cd build-${GEM_ARCH} && cmake -DCOMPILER=${GEM_COMPILER} ${GEM_GIT_DIR}/project )

.PHONY: build
build:
	( cd build-${GEM_ARCH} && make -j )

.PHONY: work
work: 
	( cd build-${GEM_ARCH} && make work )

#à compléter
buildclean:
	./clean.sh build

buildwork:
	./clean.sh work

cleanall:
	./clean.sh all
