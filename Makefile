SHELL = /bin/bash

# Makefile with some debug options
# .initial_setup must have been sourced beforehand

default: build

MAKEFLAGS += --no-print-directory

cmake:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake ${gem_DIR} )

# Using static Intel libraries
cmake-static:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DSTATIC_INTEL=ON ${gem_DIR} )

# Outside CMC: COMPILER_SUITE may need to be specified for Intel
cmake-intel:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCOMPILER_SUITE=intel ${gem_DIR} )

# with CMAKE_BUILD_TYPE=Debug
cmake-debug:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCMAKE_BUILD_TYPE=Debug ${gem_DIR} )

# Extra debug (see extra checks defined in cmake_rpn compiler presets and in CMakeLists.txt)
cmake-debug-extra:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && cmake -DCMAKE_BUILD_TYPE=Debug -DEXTRA_CHECKS=ON ${gem_DIR} )

.PHONY: build
build:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) )

.PHONY: work
work: 
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) work )

# make clean in build directory, to remove compiler and linker generated files
clean:
	( cd build-${GEM_ARCH} && cd `/bin/pwd` && $(MAKE) clean )

# Delete the build and work directories
distclean:
	. ./.clean_all
