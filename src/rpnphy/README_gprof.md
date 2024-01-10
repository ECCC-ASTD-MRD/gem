
Compile with profiling options:

  Add -p/-pg to compiler/linker options

Run the executable:

  Make sure no cleanup is done afterward to keep the gmon.out files
  May want to add to the env:
    export GMON_OUT_PREFIX=gmon.out-

Sum all gmon.out files:

  gprof -s ${gem_DIR}/${GEM_WORK}/bin/maingemdm* $(find  ${gem_DIR}/${GEM_WORK}/RUNMOD/ -name 'gmon.out*')

Visualize:

  gprof -p ${gem_DIR}/${GEM_WORK}/bin/maingemdm* gmon.sum > gmon.sum-p
  gprof -q ${gem_DIR}/${GEM_WORK}/bin/maingemdm* gmon.sum > gmon.sum-q


Ref:
* https://www.nas.nasa.gov/hecc/support/kb/using-gprof-for-performance-analysis_671.html
* https://web.archive.org/web/20200226032601/http://shwina.github.io/2014/11/profiling-parallel