#include <hwloc.h>
#include <stdio.h>

int cpu_per_numa() {

  hwloc_topology_t topology;
  hwloc_topology_init(&topology);
  hwloc_topology_load(topology);

  int depth = hwloc_get_type_depth(topology, HWLOC_OBJ_CORE);
  int nbcpu_total = hwloc_get_nbobjs_by_depth(topology, depth);

  int depth_numa = hwloc_get_type_depth(topology, HWLOC_OBJ_NUMANODE);
  int nb_numa = hwloc_get_nbobjs_by_depth(topology, depth_numa);

  hwloc_topology_destroy(topology);

  return nbcpu_total / nb_numa;
}
