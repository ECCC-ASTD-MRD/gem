# Add the target distclean to delete everything that cmake created
add_custom_target (distclean @echo cleaning for source distribution)
add_custom_command(
    DEPENDS clean
    COMMENT "distribution clean"
    COMMAND find
    ARGS    ${CMAKE_CURRENT_BINARY_DIR} "\\(" -name CMakeCache.txt
            -o -name cmake_install.cmake
            -o -name Makefile
            -o -name CMakeFiles
            -o -name ${BUILD}
            -o -name bin-${ARCH}
            -o -name install_manifest.txt
            "\\)" -print | xargs rm -fr
    TARGET  distclean
)
