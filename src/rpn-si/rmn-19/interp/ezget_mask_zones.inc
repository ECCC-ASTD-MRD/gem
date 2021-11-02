   subroutine qqq_ezget_mask_zones(mask_zones, x, y, ni_out, nj_out, mask_in, ni_in, nj_in)
   implicit none

   integer :: ni_out, nj_out, ni_in, nj_in
   integer :: mask_zones(ni_out, nj_out),mask_in(ni_in, nj_in)
   real    :: x(ni_out, nj_out), y(ni_out, nj_out)
   real    :: rx, ry
   integer :: i,j, k, l, ix, iy, nix, niy, nmissing

   integer :: codes(8)
   integer :: original_mask, all_pts_present, three_points_present, two_points_present, &
      one_point_present, one_point_missing, all_pts_missing, nearest_point_missing, outside_src_grid

   all_pts_missing           = 0
   one_point_present         = 1
   two_points_present        = 2
   three_points_present      = 3
   all_pts_present           = 4
   nearest_point_missing     = 5
   original_mask             = 6
   outside_src_grid          = 7


   do j=1,nj_out
      do i=1,ni_out
         ix = int(x(i,j))
         iy = int(y(i,j))
         nix = nint(x(i,j))
         niy = nint(y(i,j))
         if (ix<1.or.ix>ni_in.or.iy<1.or.iy>nj_in) then
            mask_zones(i,j) = outside_src_grid
         else
            if (mask_in(nix,niy) == 0) then
               mask_zones(i,j) = nearest_point_missing
            endif
            nmissing = 0
            do k=1,2
               do l=1,2
                  if (mask_in(ix+k-1,iy+l-1) == 0) then
                     nmissing = nmissing+1
                  endif
               enddo
            enddo
            select case (nmissing)
               case (0)
                  mask_zones(i,j) = all_pts_present
               case (1)
                  mask_zones(i,j) = three_points_present
               case (2)
                  mask_zones(i,j) = two_points_present
               case (3)
                  mask_zones(i,j) = one_point_present
               case (4)
                  mask_zones(i,j) = all_pts_missing
            end select
         endif
      enddo
   enddo


   end subroutine qqq_ezget_mask_zones

   subroutine qqq_ezsint_mask(mask_out, x, y, ni_out, nj_out, mask_in, ni_in, nj_in)
   implicit none

   integer :: ni_out, nj_out, ni_in, nj_in
   integer :: mask_out(ni_out, nj_out),mask_in(ni_in, nj_in)
   real    :: x(ni_out, nj_out), y(ni_out, nj_out)
   real    :: rx, ry
   integer :: i,j, k, l, ix, iy, nix, niy, nmissing, ier

   integer ezgetopt
   external ezgetopt
   character(len=32) :: value

   ier = ezgetopt('cloud_interp_alg', value)

   mask_out = 1

   do j=1,nj_out
      do i=1,ni_out
         ix = int(x(i,j))
         iy = int(y(i,j))
         nix = nint(x(i,j))
         niy = nint(y(i,j))
         if (ix<1.or.ix>ni_in.or.iy<1.or.iy>nj_in) then
            mask_out(i,j) = 0
         else if (mask_in(nix,niy) == 0) then
            mask_out(i,j) = 0
         endif
      enddo
   enddo

   if (value(1:6) == 'linear') then
      do j=1,nj_out-1
         do i=1,ni_out-1
            if (mask_out(i,j) == 1) then
               ix = int(x(i,j))
               iy = int(y(i,j))
               if (mask_in(ix+1,iy) == 0 .or.mask_in(ix,iy+1) == 0 .or. mask_in(ix+1,iy+1) == 0) then
                  mask_out(i,j) = 0
               endif
            endif
         enddo
      enddo
   endif

   end subroutine qqq_ezsint_mask
