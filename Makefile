SHELL = /bin/bash

# Makefile for Environment Canada systems
# Make sure you update the appropriate git submodules, according to what you want to build

default: build

MAKEFLAGS += --no-print-directory

# Using installed RPN libraries (rmn, vgrid, rpncomm, tdpack)
cmake:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DWITH_SYSTEM_RPN=TRUE ${GEM_GIT_DIR} )

# Compiling everything: you need to update rpn-si libraries (rmn, vgrid, rpncomm, tdpack) submodules to do this
cmake-all:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DWITH_SYSTEM_RPN=FALSE ${GEM_GIT_DIR} )

# with CMAKE_BUILD_TYPE=Debug
cmake-debug:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCMAKE_BUILD_TYPE=Debug -DWITH_SYSTEM_RPN=TRUE $${GEM_GIT_DIR} )

# Extra debug (see extra checks defined in cmake_rpn compiler presets and in CMakeLists.txt)
cmake-debug-extra:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCMAKE_BUILD_TYPE=Debug -DWITH_SYSTEM_RPN=TRUE -DEXTRA_CHECKS=ON ${GEM_GIT_DIR} )

.PHONY: build
build:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) )

.PHONY: work
work: 
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) work )

package: 
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) package )

clean:
	. ./.initial_setup

distclean:
	. ./.clean_all


