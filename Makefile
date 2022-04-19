SHELL = /bin/bash

# Makefile with some debug options
# .initial_setup must have been sourced beforehand

default: build

cmake:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake ${GEM_GIT_DIR} )

# Outside CMC: COMPILER_SUITE must be specified for Intel
cmake-intel:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCOMPILER_SUITE=intel ${GEM_GIT_DIR} )

# with CMAKE_BUILD_TYPE=Debug
cmake-debug:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCMAKE_BUILD_TYPE=Debug ${GEM_GIT_DIR} )

# Extra debug (see extra checks defined in cmake_rpn compiler presets and in CMakeLists.txt)
cmake-debug-extra:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCMAKE_BUILD_TYPE=Debug -DEXTRA_CHECKS=ON ${GEM_GIT_DIR} )

.PHONY: build
build:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) )

.PHONY: work
work: 
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) work )
