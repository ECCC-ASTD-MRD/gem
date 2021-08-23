SHELL = /bin/bash

#Makefile for Environment Canada systems

default: build

cmake-gnu-debug: ; $(MAKE) COMPILER_SUITE=gnu gnu-debug
cmake-gnu-debug-plus: ; $(MAKE) COMPILER_SUITE=gnu gnu-debug-plus

cmake:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake ${GEM_GIT_DIR} )

gnu-debug:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCMAKE_BUILD_TYPE=Debug ${GEM_GIT_DIR} )

gnu-debug-plus:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCMAKE_BUILD_TYPE=Debug -DEXTRA_DEBUG=1 ${GEM_GIT_DIR} )

.PHONY: build
build:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) )

.PHONY: work
work: 
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) work )
