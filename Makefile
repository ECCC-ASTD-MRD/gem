SHELL = /bin/bash

#Makefile for Environment Canada systems

COMPILER_SUITE ?= "intel"

default: build

cmake-intel: ; $(MAKE) COMPILER_SUITE=intel cmake
cmake-gnu: ; $(MAKE) COMPILER_SUITE=gnu cmake
cmake-gnu-debug: ; $(MAKE) COMPILER_SUITE=gnu gnu-debug
cmake-gnu-debug-plus: ; $(MAKE) COMPILER_SUITE=gnu gnu-debug-plus

cmake:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCOMPILER_SUITE=${COMPILER_SUITE} ${GEM_GIT_DIR} )

gnu-debug:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCOMPILER_SUITE=gnu -DCMAKE_BUILD_TYPE=Debug ${GEM_GIT_DIR} )

gnu-debug-plus:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCOMPILER_SUITE=gnu -DCMAKE_BUILD_TYPE=Debug -DEXTRA_DEBUG=1 ${GEM_GIT_DIR} )

.PHONY: build
build:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) )

.PHONY: work
work: 
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) work )
