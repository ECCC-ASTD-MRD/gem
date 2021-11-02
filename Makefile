SHELL = /bin/bash

#Makefile for Environment Canada systems

default: build

debug: ; $(MAKE) cmake-debug
debug-plus: ; $(MAKE) cmake-debug-plus

cmake-with-system-rpn:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DWITH_SYSTEM_RPN=TRUE ${GEM_GIT_DIR} )

cmake:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DWITH_SYSTEM_RPN=FALSE ${GEM_GIT_DIR} )

cmake-debug:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCMAKE_BUILD_TYPE=Debug ${GEM_GIT_DIR} )

cmake-debug-plus:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCMAKE_BUILD_TYPE=Debug -DEXTRA_DEBUG=1 ${GEM_GIT_DIR} )

.PHONY: build
build:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) )

.PHONY: work
work: 
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) work )
