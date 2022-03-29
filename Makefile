SHELL = /bin/bash

# Makefile with some debug options
# .initial_setup must have been sourced beforehand

COMPILER_SUITE ?= "intel"

default: build

cmake-intel: ; $(MAKE) COMPILER_SUITE=intel cmake
cmake-gnu: ; $(MAKE) COMPILER_SUITE=gnu cmake
cmake-intel-debug: ; $(MAKE) COMPILER_SUITE=intel cmake-debug
cmake-intel-debug-extra: ; $(MAKE) COMPILER_SUITE=intel cmake-debug-extra
cmake-gnu-debug: ; $(MAKE) COMPILER_SUITE=gnu cmake-debug
cmake-gnu-debug-extra: ; $(MAKE) COMPILER_SUITE=gnu cmake-debug-extra

cmake:
	( cd build-${GEM_ARCH} && cmake -DCOMPILER_SUITE=${COMPILER_SUITE} ${GEM_GIT_DIR} )

# with CMAKE_BUILD_TYPE=Debug
cmake-debug:
	( cd build-${GEM_ARCH} && cmake -DCOMPILER_SUITE=${COMPILER_SUITE} -DCMAKE_BUILD_TYPE=Debug ${GEM_GIT_DIR} )

# Extra debug (see extra checks defined in cmake_rpn compiler presets and in CMakeLists.txt)
cmake-debug-extra:
	( cd build-${GEM_ARCH} && cmake -DCOMPILER_SUITE=${COMPILER_SUITE} -DCMAKE_BUILD_TYPE=Debug -DEXTRA_CHECKS=ON ${GEM_GIT_DIR} )

.PHONY: build
build:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) )

.PHONY: work
work: 
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) work )
