SHELL = /bin/bash

#Makefile for Environment Canada systems

COMPILER_suite ?= "intel"

default: build

cmake-intel: ; $(MAKE) COMPILER_suite=intel cmake
cmake-gnu: ; $(MAKE) COMPILER_suite=gnu cmake
cmake:
	( cd build-${GEM_ARCH} && cmake -DCOMPILER=${COMPILER_suite} ${GEM_GIT_DIR}/project )

.PHONY: build
build:
	( cd build-${GEM_ARCH} && make -j )

.PHONY: work
work: 
	( cd build-${GEM_ARCH} && make work )
