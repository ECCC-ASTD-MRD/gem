SHELL = /bin/bash

#Makefile for Environment Canada systems

COMPILER_suite ?= "intel"

default: build

cmake-intel: ; $(MAKE) COMPILER_suite=intel cmake
cmake-gnu: ; $(MAKE) COMPILER_suite=gnu cmake
cmake-gnu-debug: ; $(MAKE) COMPILER_suite=gnu gnu-debug
cmake-gnu-debug-plus: ; $(MAKE) COMPILER_suite=gnu gnu-debug-plus

cmake:
	( cd build-${GEM_ARCH} && cmake -DCOMPILER=${COMPILER_suite} ${GEM_GIT_DIR}/project )

gnu-debug:
	( cd build-${GEM_ARCH} && cmake -DCOMPILER=gnu -DCMAKE_BUILD_TYPE=Debug ${GEM_GIT_DIR}/project )

gnu-debug-plus:
	( cd build-${GEM_ARCH} && cmake -DCOMPILER=gnu -DCMAKE_BUILD_TYPE=Debug -DEXTRA_DEBUG=1 ${GEM_GIT_DIR}/project )

.PHONY: build
build:
	( cd build-${GEM_ARCH} && make -j )

.PHONY: work
work: 
	( cd build-${GEM_ARCH} && make work )
