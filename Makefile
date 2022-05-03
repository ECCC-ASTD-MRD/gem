SHELL = /bin/bash

# Makefile for Environment Canada systems
# Make sure you update the appropriate git submodules, according to what you want to build

default: build

debug: ; $(MAKE) cmake-debug
debug-plus: ; $(MAKE) cmake-debug-plus

# Using installed RPN libraries (rmn, vgrid, rpncomm, massv, tdpack)
cmake-with-system-rpn:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DWITH_SYSTEM_RPN=TRUE ${GEM_GIT_DIR} )

# Compiling everything
cmake:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DWITH_SYSTEM_RPN=FALSE ${GEM_GIT_DIR} )

# with CMAKE_BUILD_TYPE=Debug
cmake-debug:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCMAKE_BUILD_TYPE=Debug ${GEM_GIT_DIR} )

# Extra debug (see extra checks defined in cmake_rpn compiler presets and in CMakeLists.txt)
cmake-debug-plus:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCMAKE_BUILD_TYPE=Debug -DEXTRA_DEBUG=1 ${GEM_GIT_DIR} )

.PHONY: build
build:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) )

.PHONY: work
work: 
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) work )
