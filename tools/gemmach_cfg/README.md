# Optional Configuration for GEM-MACH Integration Test

`tools/gemmach_cfg` directory within the GEM-MACH repository is created to be a place-holder for local version of configuration for GEM-MACH integration test that can be utilized to overwrite the configuration provided in the common remote repository (i.e. GitLab).

If this directory contains any of the three configuration files (i.e. `gm_phy_intable`, `gem_settings.nml`, and/or `outcfg.out`) are saved in `tools/gemmach_cfg`, they will be used instead of the ones provided in the common remote repository.

