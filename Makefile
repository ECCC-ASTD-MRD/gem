SHELL = /bin/bash

#Makefile for Environment Canada systems

COMPILER_name ?= "intel"

default: build

cmake-intel: ; $(MAKE) COMPILER_name=intel cmake
cmake-gnu: ; $(MAKE) COMPILER_name=gnu cmake
cmake:
	( cd build-${GEM_ARCH} && cmake -DCOMPILER=${COMPILER_name} ${GEM_GIT_DIR}/project )

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
