subroutine cpl_declare_grid_information()
use iso_fortran_env
use iris_mod
use phygridmap, only: drv_glb_ni, drv_glb_nj, drv_glb_hx, drv_glb_hy, drv_glb_i0, drv_glb_j0, drv_lcl_ni, drv_lcl_nj, phy_yinyang_S, phy_yinyang_L
use rmn_gmm
use cplocn_mod
implicit none

      integer :: ierr
      integer :: j0, j1
      integer :: i0, i1
      integer :: gnj

      ! Iris is inclusive with what it does with tiles, that is for a tile defined
      ! by (i0,j0), (i1,j1), the (i1,j1) corner is considered to be part of the
      ! tile (this is different from the C mentality.  If we have a tile that is
      ! 319x131 starting at (1,1), then the other corner is (1+319-1,1+131-1) = (319,131)
      ! and we see that the set {1, ..., 312} has 319 elements.
      i0 = drv_glb_i0
      i1 = i0 + drv_lcl_ni - 1
      j0 = drv_glb_j0
      j1 = j0 + drv_lcl_nj - 1
      gnj = drv_glb_nj
      if(phy_yinyang_L) then
        ! drv_glb_ni is the size of a Yin or Yan subgrid rather than the size of the
        ! total grid
        gnj = 2 * drv_glb_nj
        if(phy_yinyang_S == "YAN") then
          ! Working with Yin-Yan grid: Yin has an G_NIxG_NJ grid and Yan has a
          ! G_NIxG_NJ grid.  They need to be put together by Iris into one
          ! G_NIx(2*G_NJ) grid.  Let (ti,tj) be indices of tiles in a Yin or Yan
          ! subgrid and (TI,TJ) be indices of tiles in the combined grid For
          ! Yin, TI == ti, TJ == tj, but for Yan, we will have TI = ti, TJ = tj
          ! + tnj
          j0 = j0 + drv_glb_nj
          j1 = j1 + drv_glb_nj
        endif
      endif

      write (*, '("GEM:", a, "iris%model_grid(ATMOS, ", i4, ", ", i4, ", ", i4, ", ", i4, ", ", i4,", ", i4,", ",i4, ", ",i4,", ",i4, ", ",i4,", ",i4,")")') &
          phy_yinyang_S, drv_glb_ni, gnj, 1, drv_glb_hx, drv_glb_hy, i0, i1, j0, j1, 0, 0
      iris_grid = Iris%Model_Grid("ATMOS", drv_glb_ni, gnj, 1, drv_glb_hx, drv_glb_hy, i0, i1, j0, j1, 0, 0)

      if(.not. c_associated(iris_grid)) then
        write (*,*) "Error calling Iris%Model_Grid in cplocn_init"
        return
      endif

end
