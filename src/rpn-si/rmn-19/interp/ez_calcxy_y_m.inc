   subroutine ez_calcxy_y_m(wts, idxs, x, y, gdout_lat, gdout_lon, gdout_masque, gdin_lat, gdin_lon, gdin_masque, ni_src, nj_src, ni_dst, nj_dst, num_wts)
   implicit none

   integer ni_src, nj_src, ni_dst, nj_dst, num_wts
   real :: gdin_lon(ni_src,nj_src),gdin_lat(ni_src,nj_src), dist, total_wgt
   real :: gdout_lon(ni_dst, nj_dst), gdout_lat(ni_dst, nj_dst)
   real, dimension(:,:,:) :: wts(ni_dst, nj_dst, num_wts), x(ni_dst, nj_dst), y(ni_dst, nj_dst)
   real, dimension(:) :: tmpwgts(num_wts)
   integer :: idxs(ni_dst, nj_dst, num_wts), gdout_masque(ni_dst,nj_dst), gdin_masque(ni_src, nj_src)
   integer :: tmp_idxs(num_wts,2)
   real rx, ry
   integer i,j,k,ii,jj,iz
   real lon_min, lon_max, lat_min, lat_max
   integer valmax, locmax
   integer ier, i1, j1, k1, exception
   logical point_found

   integer ezgetval
   external ezgetval

   real dist_thresh

   lon_min = minval(gdin_lon)
   lon_max = maxval(gdin_lon)
   lat_min = minval(gdin_lat)
   lat_max = maxval(gdin_lat)
   gdout_masque = 1
   where (gdout_lon < lon_min) gdout_masque = 0
   where (gdout_lon > lon_max) gdout_masque = 0
   where (gdout_lat < lat_min) gdout_masque = 0
   where (gdout_lat > lat_max) gdout_masque = 0

   wts = 1.0e30
   idxs = 0
   ier = ezgetval('missing_distance_threshold', dist_thresh)

   do j=1,nj_dst
      do i=1,ni_dst
         k = i + (j-1) * ni_dst
         if (gdout_masque(i,j) /= 0) then
            tmpwgts = 1.0e30
            locmax=1
            tmp_idxs = 0
            point_found = .false.
            do jj=1,nj_src
               do ii=1,ni_src
                  if (gdin_masque(ii,jj) == 1) then
                     rx = (gdout_lon(i,j)- gdin_lon(ii,jj))
                     ry = (gdout_lat(i,j)- gdin_lat(ii,jj))
                     if (abs(rx) > 180.0) then
                        rx = abs(abs(rx) - 360.0)
                     endif
                     if (abs(ry) < dist_thresh .and. abs(rx) < dist_thresh) then
                        dist = rx*rx + ry*ry
                        point_found = .true.
                        if (dist < tmpwgts(locmax)) then
                           tmpwgts(locmax) = dist
                           idxs(i,j,locmax) = ii + (jj-1)*ni_src
                           tmp_idxs(locmax, 1) = ii
                           tmp_idxs(locmax, 2) = jj
                           locmax = maxloc(tmpwgts, 1)
                        endif
                     endif
                  endif
               enddo
            enddo
            if (.not.point_found) then
               gdout_masque(i,j) = 0
            endif
            if (gdout_masque(i,j) == 1) then
               call inside_or_outside(gdout_masque(i,j),x(i,j), y(i,j), gdout_lat(i,j),gdout_lon(i,j), gdin_lat, &
                                       gdin_lon, ni_src, nj_src, tmpwgts, tmp_idxs, num_wts)
            endif
            if (gdout_masque(i,j) == 1) then
               do iz=1,num_wts
                  wts(i,j,iz) = tmpwgts(iz)
               enddo
               if (wts(i,j,1) > 6371000.0) then
                  wts(i,j,1:num_wts) = 1.0e30
                  gdout_masque(i,j) = 0
               else
                  total_wgt = 0.0
                  do iz=1,num_wts
                     if (wts(i,j,iz) < 1.0e-10) then
                        print *, 'wts(i,j,iz) < 1.0e-10',i,j,iz,wts(i,j,iz)
                     endif
                     wts(i,j,iz) = max(wts(i,j,iz), 1.0e-10)
                     wts(i,j,iz) = 1.0 / wts(i,j,iz)
                  enddo
                  do iz=1,num_wts
                     total_wgt = total_wgt + wts(i,j,iz)
                  enddo
                  do iz=1,num_wts
                     wts(i,j,iz) = wts(i,j,iz) / total_wgt
                  enddo
               endif
            endif
         endif
      enddo
   enddo

   return
   end subroutine ez_calcxy_y_m

